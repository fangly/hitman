////////////////////////////////////////////////////////////////////////////////

// Cluster and annotate OTUs, make a rarefied, copy-corrected microbial profile and graph

// Input is:
// 1/ *.fna: FASTA reads, clean, trimmed to the same length, and named by sample
// 2/ *.fna: FASTA file of primers used (called 'fwd' and 'rev')

// Usage example: bpipe run pipeline_analyze.groovy sediments_QA.fna fwd_primer.fna

// Options:
// Number of threads
NOF_THREADS=12
// Flag to prevent singleton sequences to constitute an OTU
EXCLUDE_SINGLETONS=1
// OTU identity (from 0 to 1)
OTU_MIN_IDENTITY=0.97
// Flag to remove reference-based chimeras
EXCLUDE_CHIMERAS=1
// Database for chimera detection
CHIMERA_DB="/srv/db/gold/micrombiomeutils_r20110519_GOLD.fna"
// Parameters for taxonomic assignment
TAXO_DB="/srv/db/merged_gg_silva/gg201210_silva115/merged_gg_silva_99clust_seqs.fna"
TAXO_STR="/srv/db/merged_gg_silva/gg201210_silva115/merged_gg_silva_99clust_taxo.txt"
// Rule of thumb: >97% indicates the same species, >95% the same genus
TAXO_MIN_IDENTITY=0.90
//TAXO_MIN_IDENTITY=0.95
// Samples to merge
MERGE_SAMPLES=""
// OTUs to remove (case-insensitive)
REMOVE_OTUS="Eukaryota"
// Rarefaction depth (number of reads per sample)
RAREFACTION_DEPTH=1000
// CopyRighter gene copy number database
GCN_DB="/srv/db/copyrighter/ssu_img40_gg201210_merged.txt"

////////////////////////////////////////////////////////////////////////////////


about title: "16S analysis pipeline - OTU analysis" 


@Transform("otus")
dereplicate = {
   doc title: "Collapse duplicate reads, record their abundance and sort them by abundance"
   // http://drive5.com/usearch/manual/derep_fulllength.html
   // http://drive5.com/usearch/manual/sortbysize.html
   // usearch -derep_fulllength filtered.fasta -output derep.fasta -sizeout
   // usearch -sortbysize derep.fasta -output derep2.fasta
   exec """
      module load usearch/7.0.1001 &&
      TEMP_FILE=`mktemp tmp_dereplicate_XXXXXXXX.otus` &&
      usearch -derep_fulllength $input.fna -output \$TEMP_FILE -sizeout -threads $NOF_THREADS &&
      usearch -sortbysize \$TEMP_FILE -output $output.otus &&
      rm \$TEMP_FILE
   """
}


@Filter("rm_singletons")
rm_singletons = {
   doc title: "Remove singleton reads"
   if (EXCLUDE_SINGLETONS == 1) {
      // http://drive5.com/usearch/manual/sortbysize.html
      // usearch -sortbysize derep.fasta -output derep2.fasta -minsize 2
      exec """
         NOF_SINGLETONS=`grep '^>' $input.otus | grep -c 'size=1;'` &&
         echo "Removing \$NOF_SINGLETONS singletons!" &&
         module load usearch/7.0.1001 &&
         usearch -sortbysize $input.otus -output $output.otus -minsize 2
      """
   } else {
      exec """
         NOF_SINGLETONS=`grep '^>' $input.otus | grep -c 'size=1;'` &&
         echo "Not removing \$NOF_SINGLETONS singletons..." &&
         cp $input.otus $output.otus
      """
   }
}


@Filter("otu_clustering")
otu_clustering = {
   doc title: "Cluster reads into OTUs, remove de novo chimeras"
   // http://www.drive5.com/usearch/manual/cluster_otus.html
   // http://drive5.com/python/fasta_number_py.html
   // usearch -cluster_otus derep2.fasta -otus otus.fasta -otuid 0.97
   exec """
      module load usearch/7.0.1001 &&
      TEMP_FILE=`mktemp tmp_otu_clustering_XXXXXXXX.otus` &&
      usearch -cluster_otus $input.otus -otus \$TEMP_FILE -otuid $OTU_MIN_IDENTITY &&
      fasta_number.py \$TEMP_FILE OTU_ > $output.otus &&
      rm \$TEMP_FILE
   """
}


@Filter("rm_chimeras")
rm_chimeras = {
   doc title: "Chimera filtering using a reference database"
   if (EXCLUDE_CHIMERAS == 1) {
      // http://www.drive5.com/usearch/manual/uchime_ref.html
      // usearch -uchime_ref otus1.fa -db $d/gold.fa -strand plus -nonchimeras otus2.fa
      exec """
         echo "Removing chimeras!" &&
         module load usearch/7.0.1001 &&   
         usearch -uchime_ref $input.otus -db $CHIMERA_DB -threads $NOF_THREADS -strand plus -nonchimeras $output.otus
      """
   } else {
      exec """
         echo "Not removing chimeras..." &&
         cp $input.otus $output.otus
      """
   }
}


uparse_otu = segment {
   // Run the UPARSE OTU clustering pipeline
   // Input is amplicon reads in FASTA (*.fna)
   // Output is OTU reference sequences in FASTA (*.otus)
   dereplicate + rm_singletons + otu_clustering + rm_chimeras
}


@Transform("uc")
otu_mapping = {
   doc title: "Map each reads to an OTU (if possible)"
   // The uc2otutab.py script expects reads headers with a tag called barcodelabel=
   // http://drive5.com/usearch/manual/usearch_global.html
   // usearch -usearch_global reads.fa -db otus.fa -strand plus -id 0.97 -uc readmap.uc
   // find how many reads did not match (lines that start with an N, i.e. no-hit records)
   // grep -c "^N" readmap.uc
   exec """
      module load usearch/7.0.1001 &&
      TEMP_FILE=`mktemp tmp_otu_mapping_XXXXXXXX.fna` &&
      perl -p -e 's/^>(\\S+)_(\\S+)/>\$2;barcodelabel=\$1;/' $input.fna > \$TEMP_FILE &&
      usearch -usearch_global \$TEMP_FILE -db $input.otus -strand plus -id $OTU_MIN_IDENTITY -threads $NOF_THREADS -uc $output.uc &&
      rm \$TEMP_FILE &&
      NO_MATCH=`grep -c "^N" $output.uc` &&
      echo "Number of non-matching reads: \$NO_MATCH"
   """
}


@Transform("generic")
uclust2generic = {
   doc title: "Convert microbial profile from UCLUST format to generic site-by-species"
   exec """
      module load usearch/7.0.1001 &&
      uc2otutab.py $input.uc > $output.generic
   """
}


@Transform("biom")
generic2biom = {
   doc title: "Convert microbial profile from site-by-species generic format to biom"
   exec """
      module load bio-community &&
      bc_convert_files -if $input.generic -of biom -op $output.prefix
   """
}


@filter("desc2id")
desc2id = {
   doc title: "Give OTU the ID that is in their description field"
   exec """
      module load bio-community &&
      bc_convert_ids -if $input -ma desc -op $output.prefix
   """
}


make_otu_table = segment {
   // Convert microbial profile from UCLUST format to generic, and then to biom
   uclust2generic + generic2biom + desc2id
}


trim_db = {
   doc title: "Trim a database of sequences to their amplicon and save it globally"
   exec """
      module load bioperl &&
      TRIM_LEN=`extract_first_seqs --input $input.fna --number 1 | get_seq_length | cut -f 2` &&
      echo "Trim length: \$TRIM_LEN" &&
      FWD_PRIMER=`extract_seqs_by_name --input $input2.fna --string fwd | convert_seq_format --format raw` &&
      echo "Fwd primer: \$FWD_PRIMER" &&
      REV_PRIMER=`extract_seqs_by_name --input $input2.fna --string rev | convert_seq_format --format raw` &&
      echo "Rev primer: \$REV_PRIMER" &&
      TAXO_DB=${TAXO_DB} &&
      TRIMMED_TAXO_DB=\${TAXO_DB%.*}_trimmed_\${FWD_PRIMER}_\${TRIM_LEN}_bp.fna &&
      if [ ! -e \$TRIMMED_TAXO_DB ]; then
         echo "Trimming reference sequences in \$TRIMMED_TAXO_DB" &&
         module load emboss &&
         extract_amplicons -i $TAXO_DB -f \$FWD_PRIMER -r \$REV_PRIMER -t \$TRIM_LEN -e 1 -o \$TRIMMED_TAXO_DB;
      else
         echo "Re-using trimmed reference sequences in \$TRIMMED_TAXO_DB";
      fi &&
      ln -s \$TRIMMED_TAXO_DB $output.fna
   """
}


usearch_global = {
   doc title: "Assign taxonomy using USEARCH global alignments"
   exec """
      module load usearch/7.0.1001 &&
      usearch -usearch_global $input.otus -db $input.fna -id $TAXO_MIN_IDENTITY -threads $NOF_THREADS -strand plus -blast6out $output.blast
   """
   // -maxhits (default 0; i.e. ignored)
   //    This indicates the maximum number of hits written to the output files.
   //    This is not a termination condition; the search continues after the
   //    maximum number of hits is exceeded. When the search terminates, hits
   //    are sorted by decreasing identity or increasing E-value and up to
   //    maxhits are written. Default is 0, meaning that the option is ignored.
   // ‑maxaccepts (default 1) & -maxrejects (default 32)
   //    These options stop the search for a given query sequence if a given
   //    number of accepts or rejects have occurred
}


assign_taxo_global = segment {
   trim_db + usearch_global
}


@filter("id2tax")
id2tax = {
   doc title: "Replace IDs in microbial profiles by taxonomically-informative IDs"
   exec """
      module load bio-community &&
      bc_convert_ids -if $input.biom -bf $input.blast -op $output.prefix
   """
}


@filter("desc2tax")
desc2tax = {
   doc title: "Add taxonomic lineage to an microbial profiles"
   exec """
      module load bio-community &&
      bc_add_taxonomy -if $input.biom -tf $TAXO_STR -op $output.prefix
   """
}


add_taxonomy = segment {
   id2tax + desc2tax
}


//
// Remove samples
//
// Merge samples
//
// Remove euks
//


@filter("rarefy")
rarefy = {
   doc title: "Rarefy microbial profiles to 1,000 reads"
   exec """
      module load bio-community &&
      bc_rarefy -if $input -nr inf -ss $RAREFACTION_DEPTH -op $output.prefix
   """
}


// CopyRighter

// Heatmap

// alpha diversity

// beta diversity



success = {
   exec """
      echo "---> Success!" &&
      echo "$input1" &&
      echo "$input2"
   """
}


Bpipe.run {
   uparse_otu + otu_mapping + [make_otu_table, assign_taxo_global] + add_taxonomy + rarefy + success
}
