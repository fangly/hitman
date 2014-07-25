// Shared modules for Bpipe


// Block size to use for GNU Parallel
BLOCK="10M"


gunzip = {
   doc title: "Uncompress GZIP files",
       desc:  """This step is skipped if input is not compressed. Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   def ext = (input =~ /\.([^\.]*)$/)[0][1]
   def out = output.prefix
   if (ext == "gz") {
      out = out.replaceFirst(~/\.[^\.]+$/, '')
      produce(out) {
         exec """
            echo "Uncompressing GZIP file" &&
            pigz -d -p $threads -c $input.gz > $out
         """
      }
   } else {
      produce(out) {
         exec """
            echo "Skipping GZIP decompression" &&
            cp $input $out
         """
      }
      // TODO: When bug #90 is fixed, replace by: forward input
      // See code.google.com/p/bpipe/issues/detail?id=90
   }
}


gzip = {
   doc title: "Compress GZIP files",
       desc:  """This step is skipped if input is compressed. Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   def ext = (input =~ /\.([^\.]*)$/)[0][1]
   def out = output.prefix
   if (ext == "gz") {
      produce(out) {
         exec """
            echo "Skipping GZIP compression" &&
            cp $input $out
         """
      }
      // TODO: When bug #90 is fixed, replace by: forward input
      // See code.google.com/p/bpipe/issues/detail?id=90

   } else {
      out = out + ".gz"
      produce(out) {
         exec """
            echo "Compressing to GZIP file" &&
            pigz -p $threads -c $input > $out
         """
      }
   }
}


split_libraries = {
   doc title: "Split FASTQ libraries by MID using fastq-multx",
       desc:  """Uses a 2-column tabular mapping file. This step is skipped if
                 no mapping file is provided. Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   def skip = 1
   for ( file in inputs ) {
      if ( (file =~ /\.([^\.]*)$/)[0][1] == "mapping" ) {
         skip = 0
         break
      }
   }
   if (skip == 1) {
      def source   = new File(input.fastq).getName()
      def basename = source.replaceFirst(~/\..+$/, '')
      def output   = input.prefix+"."+basename+".fastq"
      produce( output ) {
         // Simply prepare the file for the rename_seqs step
         exec """
            echo "Skipping library split by MID" &&
            ln -f -s $source $output
         """
      }
   } else {
      //http://code.google.com/p/ea-utils/wiki/FastqMultx
      produce("${input.prefix}.*.fastq") {
         exec """
            echo "Splitting libraries by MID" &&
            module load ea_utils &&
            fastq-multx -B $input.mapping $input.fastq -o ${input.prefix}.%.fastq -b -m 0 &&
            rm ${input.prefix}.unmatched.fastq
         """
      }
   }
   // Alternatives:
   // - QIIME split_libraries_fastq.py, but does a bunch of extra stuff that's not needed
   //   https://github.com/qiime/qiime/blob/master/scripts/split_libraries_fastq.py
   // - FASTX_toolkit fastx_barcode_splitter, but no support for variable-length MIDs.
   //   http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage
   // - NGOpt splitBC, but does not support different MID lengths
   //   http://code.google.com/p/ngopt/source/browse/trunk/tools/splitbc
   // - barcode_splitter.py, but needs a weird "index" file
   //   https://gist.github.com/dgrtwo/3725741
   // - ngsutils barcode_split.py
   //   https://github.com/ngsutils/ngsutils/blob/master/ngsutils/fastq/barcode_split.py
   // - sabre
   //   https://github.com/najoshi/sabre
}


sff2fastq = {
   doc title: "Convert an SFF file to FASTQ",
       desc:  """Do nothing if input is in another format. Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // sff2fastq: https://github.com/indraniel/sff2fastq
   def fmt = null
   for ( file in inputs ) {
      def ext = (file =~ /\.([^\.]*)$/)[0][1]
      if (ext == "sff") {
         fmt = ext
         break
      }
   }
   if (fmt == "sff") {
      transform("fastq") {
         exec """
            echo "Converting SFF to FASTQ" &&
            module load sff2fastq &&
            sff2fastq -o $output.fastq $input.sff
         """
      }
   } else {
      exec """
         echo "Skipping conversion from SFF to FASTQ"
      """
      forward input
   }
}


@Transform("fna")
fastq2fasta = {
   doc title: "Convert a FASTQ file to FASTA",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load gnu_parallel &&
      module load fastx_toolkit &&
      cat $input.fastq | parallel -j $threads --block $BLOCK -L 4 -k --pipe "
         fastq_to_fasta -n -Q 33
      " > $output.fna
   """
}


@Transform("fastq")
fasta2fastq = {
   doc title: "Convert a FASTA file to FASTQ, with artificial quality scores",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // NEBC fasta_to_fastq.py http://goo.gl/DzmgcD
   // Usage: fasta_to_fastq.py input.fna
   // where the qual file must be named input.fna.qual
   exec """
      module load gnu_parallel &&
      module load biopython &&
      cat $input.fna | parallel -j $threads --block $BLOCK -L 4 -k --pipe "
         F=\$(mktemp -u tmp_convert_XXXXXXXX.{#}) &&
         cat > \\${F}.fna &&
         fake_qual \\${F}.fna 1> /dev/null &&
         mv -f \\${F}.fna.qual \\${F}.qual &&
         fasta_to_fastq.py \\${F}.fna &&
         cat \\${F}.fastq &&
         rm \\${F}.fna \\${F}.fastq \\${F}.qual
      " > $output.fastq
   """
}


merge_pairs = {
   doc title: "Merge pairs of FASTQ reads that overlap using PEAR",
       desc:  """Skip if not exactly two FASTQ files were given. Parameters:
                    'keep_fwd_singles', boolean to keep forward reads that could not be merged (default: 0)
                    'keep_rev_singles', boolean to keep reverse reads that could not be merged (default: 0)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   def count = 0
   for ( file in inputs ) {
      if ( (file =~ /\.([^\.]*)$/)[0][1] == 'fastq' ) {
         count = count + 1
      }
   }
   if (count != 2) {
      exec """
         echo "Skipping paired-end merging"
      """
      forward input
   } else {
      produce(input.prefix+".merge_pairs.assembled.fastq") {
         // http://www.exelixis-lab.org/web/software/pear
         exec """
            echo "Merging paired-end reads" &&
            module load pear &&
            pear -f $input1.fastq -r $input2.fastq -o ${input.prefix}.merge_pairs -p 0.01 -v 1 -b 33 -j $threads &&
            NONMERGED=`grep -c '^@' ${input.prefix}.merge_pairs.unassembled.forward.fastq` &&
            echo "Approx. $NONMERGED pairs of reads could not be merged" &&
            if [ "$keep_fwd_singles" = "1" ]; then
               echo "Adding forward singles";
               cat ${input.prefix}.merge_pairs.unassembled.forward.fastq >> ${input.prefix}.merge_pairs.assembled.fastq;
            fi &&
            if [ "$keep_rev_singles" = "1" ]; then
               echo "Adding forward singles";
               cat ${input.prefix}.merge_pairs.unassembled.forward.fastq >> ${input.prefix}.merge_pairs.assembled.fastq;
            fi
         """
      }
   }
   // Alternatives:
   //   Pandaseq, nice, but has more false-positives and different type of quality scores (lower)
   //   USEARCH fastq_mergepairs: http://www.drive5.com/usearch/manual/fastq_mergepairs.html
   //   ea-utils fastq-join (see below)
   //   FLASH, COPE, XORRO
}

@Filter("rename_seqs")
rename_seqs = {
   doc title: "Prefix FASTQ sequence names based on their filename",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // From a file like "Gasket76.split_libraries.3833.fastq", use the "3833" prefix
   // A sequence called "HN758CG01AIWK2" would thus become "3833_HN758CG01AIWK2"
   // In comparison, QIIME would call it "3833_1"
   // Note: rely on the sequences not being wrapped
   exec """
      PREFIX=`basename $input.fastq | perl -n -e 'print((split /\\./)[-2])'` &&
      cat $input.fastq | paste - - | perl -p -e "s/^@/\\@\${PREFIX}_/; s/^\\+\\S*/+/" | tr '\\t' '\\n' > $output.fastq
   """
   // Alternatively, use FQTRIM or `usearch -fastq_filter seqs.fq ‑fastqout seqs_new.fq -relabel sample_`
}


@Filter("cat_files")
cat_files = {
   doc title: "Concatenate all files into a single one",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      cat $inputs > $output
   """
}


final_report = {
   doc title: "Report the name of output files",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   println("Outputs are : "+inputs)
   //println(inputs.dump())
}


acacia = {
   doc title: "Denoise FASTQ reads with Acacia",
       desc:  """Parameters:
                    'skip', boolean to skip this step (default: 0)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   produce(input.prefix+".acacia.fna") {
      if (skip == 1) {
         exec """
            echo "Skipping Acacia denoising" &&
            module load gnu_parallel &&
            module load fastx_toolkit &&
            cat $input.fastq | parallel -j $threads --block $BLOCK -L 4 -k --pipe "
               fastq_to_fasta -n -Q 33
            " > $output.fna
         """
      } else {
         exec """
            echo "Denoising using Acacia"                                &&
            module load acacia/1.52                                      &&
            CONF_FILE=`mktemp tmp_acacia_XXXXXXXX.conf`                  &&
            echo "FASTA=FALSE"                             > \$CONF_FILE &&
            echo "FASTA_LOCATION=null"                    >> \$CONF_FILE &&
            echo "QUAL_LOCATION=null"                     >> \$CONF_FILE &&
            echo "FASTQ=TRUE"                             >> \$CONF_FILE &&
            echo "FASTQ_LOCATION=$input.fastq"            >> \$CONF_FILE &&
            echo "OUTPUT_DIR=./"                          >> \$CONF_FILE &&
            echo "OUTPUT_PREFIX=$output.prefix"           >> \$CONF_FILE &&
            echo "ANY_DIFF_SIGNIFICANT_FOR_TWO_SEQS=TRUE" >> \$CONF_FILE &&
            echo "AVG_QUALITY_CUTOFF=0"                   >> \$CONF_FILE &&
            echo "ERROR_MODEL=Balzer"                     >> \$CONF_FILE &&
            echo "FILTER_N_BEFORE_POS=0"                  >> \$CONF_FILE &&
            echo "FLOW_KEY=TCAG"                          >> \$CONF_FILE &&
            echo "MAXIMUM_MANHATTAN_DISTANCE=13"          >> \$CONF_FILE &&
            echo "MAX_RECURSE_DEPTH=2"                    >> \$CONF_FILE &&
            echo "MAX_STD_DEV_LENGTH=100"                 >> \$CONF_FILE &&
            echo "MID_FILE=data/ROCHE_5BASE_ACACIA.mids"  >> \$CONF_FILE &&
            echo "MID_OPTION=NO_MID"                      >> \$CONF_FILE &&
            echo "MIN_FLOW_TRUNCATION=150"                >> \$CONF_FILE &&
            echo "MIN_READ_REP_BEFORE_TRUNCATION=0.0"     >> \$CONF_FILE &&
            echo "REPRESENTATIVE_SEQUENCE=Mode"           >> \$CONF_FILE &&
            echo "SIGNIFICANCE_LEVEL=-9"                  >> \$CONF_FILE &&
            echo "SPLIT_ON_MID=FALSE"                     >> \$CONF_FILE &&
            echo "TRIM_TO_LENGTH="                        >> \$CONF_FILE &&
            echo "TRUNCATE_READ_TO_FLOW="                 >> \$CONF_FILE &&
            acacia -c \$CONF_FILE                                        &&
            rm \$CONF_FILE                                               &&
            mv ${output.prefix}_all_tags.seqOut $output.fna
         """
      }
   }
}


@Filter("rm_small_seqs")
rm_small_seqs = {
   doc title: "Remove small FASTQ reads",
       desc:  """Parameters:
                    'length', the minimum length to keep""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.usadellab.org/cms/?page=trimmomatic
   exec """
      module load trimmomatic &&
      trimmomatic SE -threads $threads -phred33 $input.fastq $output.fastq MINLEN:$length
   """
}


@Filter("trunc_seqs")
trunc_seqs = {
   doc title: "Truncate FASTQ sequences to the specified length",
       desc:  """Parameters:
                    'length', the trimming length""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.usadellab.org/cms/?page=trimmomatic
   exec """
      module load trimmomatic &&
      trimmomatic SE -threads $threads -phred33 $input.fastq $output.fastq CROP:$length
   """
}


@Filter("trunc_N")
trunc_N = {
   doc title: "Truncate sequence 3' end when encountering an ambiguity",
       desc:  """Going from 5' to 3', truncate the 3' end at the first base that
                 is ambiguous (N). Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // With 454, N gets the quality score 0, i.e. '!'
   // With Illumina 1.8+, N gets the quality score 2, i.e. '#'
   // To accomodate for both, we keep only bases with a quality >= 3
   // http://www.usadellab.org/cms/?page=trimmomatic
   exec """
      module load trimmomatic &&
      trimmomatic SE -threads $threads -phred33 $input.fastq $output.fastq SLIDINGWINDOW:1:3
   """
}


@Filter("qual_trunc_seqs")
qual_trunc_seqs = {
   doc title: "Quality-based 3' end FASTQ sequence truncation",
       desc:  """Going from 5' to 3', truncate the 3' end at the first base with
                 a quality score below the threshold. All the remaining bases in
                 will have a quality equal or highger than the threshold.
                 Parameters:
                    'qual', the threshold quality score""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // Trimmomatic is a good alternative to USEARCH fastq_truncqual
   // http://www.usadellab.org/cms/?page=trimmomatic
   // http://www.drive5.com/usearch/manual/fastq_filter.html
   exec """
      module load trimmomatic &&
      trimmomatic SE -threads $threads -phred33 $input.fastq $output.fastq SLIDINGWINDOW:1:$qual
   """
}


clip_adapters = {
   doc title: "Cut library adapters from the reads using TRIMMOMATIC",
       desc:  """Parameters:
                    'adapters', a FASTA file containing the adapters""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.usadellab.org/cms/?page=trimmomatic
   // EA-util's fastq-mcf could be an alternative
   adapters = adapters.toString() == "true" ? "" : adapters   // Groovy issue 93
   if (adapters) {
      filter("clip_adapters") {
         exec """
            module load trimmomatic &&
            trimmomatic SE -threads $threads -phred33 $input.fastq $output.fastq ILLUMINACLIP:$adapters:2:30:10
         """
      }
   } else {
      exec """
         echo "Skipping adapter clipping"
      """
      forward input
   }
}


@Filter("ee_filter")
ee_filter = {
   doc title: "Remove FASTQ reads with too many expected errors",
       desc:  """Parameters:
                    'ee', the maximum number of errors""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // The calculation of expected errors uses each individual quality score
   // http://www.drive5.com/usearch/manual/fastq_filter.html
   exec """
      module load gnu_parallel &&
      module load usearch/7.0.1001 &&
      cat $input.fastq | parallel -j $threads --block $BLOCK -L 4 -k --pipe "
         F=\$(mktemp -u tmp_ee_filter_XXXXXXXX.{#}) &&
         cat > \\${F}.in &&
         usearch -fastq_filter \\${F}.in -fastqout \\${F}.out -fastq_maxee $ee -quiet 1> /dev/null &&
         cat \\${F}.out &&
         rm \\${F}.*
      " > $output.fastq
   """
}


LOG_FILE=""

seq_stats = {
   doc title: "Log basic FASTQ sequence statistics",
       desc:  """Parameters:
                    'stage', name of this step,
                    'skip', boolean to skip this step (default: 0)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://code.google.com/p/ea-utils/
   if (skip == 1) {
      exec """
         echo "Skipping computation of sequence statistics"
      """
   } else {
      if (LOG_FILE == "") {
         LOG_FILE=output.stats
      }
      produce(LOG_FILE) {
         exec """
            if [ ! -e $LOG_FILE ]; then
               echo "Creating log file $LOG_FILE" &&
               echo -e "stage\\tnum_seq\\tmin_L\\tavg_L\\tmax_L\\tmin_Q\\tavg_Q\\tmax_Q\\tnum_N" > $LOG_FILE;
            fi &&
            echo "Appending stage '$stage' to log file $LOG_FILE" &&
            module load ea_utils &&
            TEMP_FILE=`mktemp tmp_stats_${stage}_XXXXXXXX.txt` &&
            fastq-stats -D $input.fastq > \$TEMP_FILE &&
            NUM_SEQS=`grep "^reads" \$TEMP_FILE | cut -f 2` &&
            MAX_LEN=`grep "^len" \$TEMP_FILE | egrep -v "mean|min|stdev" | cut -f 2` &&
            AVG_LEN=`grep "^len mean" \$TEMP_FILE | cut -f 2` &&
            AVG_LEN=`printf "%.1f" $AVG_LEN`
            MIN_LEN=`grep "^len min" \$TEMP_FILE | cut -f 2` &&
            MAX_Q=`grep "^qual max" \$TEMP_FILE | cut -f 2` &&
            AVG_Q=`grep "^qual mean" \$TEMP_FILE | cut -f 2` &&
            AVG_Q=`printf "%.1f" $AVG_Q`
            MIN_Q=`grep "^qual min" \$TEMP_FILE | cut -f 2` &&
            PERC_N=`grep "^%N" \$TEMP_FILE | cut -f 2` &&
            TOT_BP=`grep "^total bases" \$TEMP_FILE | cut -f 2` &&
            NUM_N=`echo "\$PERC_N * \$TOT_BP / 100" | bc` &&
            echo -e "$stage\\t\$NUM_SEQS\\t\$MIN_LEN\\t\$AVG_LEN\\t\$MAX_LEN\\t\$MIN_Q\\t\$AVG_Q\\t$MAX_Q\\t$NUM_N" >> $LOG_FILE &&
            rm \$TEMP_FILE
         """
      }
   }
   forward input
}


rm_singletons = {
   doc title: "Remove singleton reads",
       desc:  """Parameters:
                    'skip', boolean to skip this step (default: 0)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   if (skip == 1) {
      exec """
         NOF_SINGLETONS=`grep '^>' $input.otus | grep -c 'size=1;'` &&
         echo "Skipping removal of \$NOF_SINGLETONS singletons"
      """
      forward input
   } else {
      filter("rm_singletons") {
         // http://drive5.com/usearch/manual/sortbysize.html
         // usearch -sortbysize derep.fasta -output derep2.fasta -minsize 2
         exec """
            NOF_SINGLETONS=`grep '^>' $input.otus | grep -c 'size=1;'` &&
            echo "Removing \$NOF_SINGLETONS singletons" &&
            module load usearch/7.0.1001 &&
            usearch -sortbysize $input.otus -output $output.otus -minsize 2
         """
      }
   }
}


@Filter("otu_clustering")
otu_clustering = {
   doc title: "Cluster reads into OTUs, remove de novo chimeras",
       desc:  """Parameters:
                    'perc', the minimum identity % required (e.g. 97%)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.drive5.com/usearch/manual/cluster_otus.html
   // http://drive5.com/python/fasta_number_py.html
   // usearch -cluster_otus derep2.fasta -otus otus.fasta -otuid 0.97
   def otuid = perc / 100
   exec """
      module load usearch/7.0.1001 &&
      TEMP_FILE=`mktemp tmp_otu_clustering_XXXXXXXX.otus` &&
      usearch -cluster_otus $input.otus -otus \$TEMP_FILE -otuid $otuid &&
      fasta_number.py \$TEMP_FILE OTU_ > $output.otus &&
      rm \$TEMP_FILE
   """
}


@Transform("generic")
uclust2generic = {
   doc title: "Convert microbial profile from UCLUST format to generic site-by-species",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load usearch/7.0.1001 &&
      uc2otutab.py $input.uc > $output.generic
   """
}


@Transform("biom")
generic2biom = {
   doc title: "Convert microbial profile from site-by-species generic format to biom",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bio-community &&
      bc_convert_files -if $input.generic -of biom -op $output.prefix
   """
}


@filter("desc2id")
desc2id = {
   doc title: "Give OTUs the ID that is in their description field",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bio-community &&
      bc_convert_ids -if $input -ma desc -op $output.prefix
   """
}


@filter("id2tax")
id2tax = {
   doc title: "Replace IDs in microbial profiles by taxonomically-informative IDs",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bio-community &&
      bc_convert_ids -if $input.biom -bf $input.blast -op $output.prefix
   """
}


@filter("desc2tax")
desc2tax = {
   doc title: "Add taxonomic lineage to an microbial profiles",
       desc:  """Parameters:
                    'db', database file containing taxonomic lineages""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bio-community &&
      bc_add_taxonomy -if $input.biom -tf $db -op $output.prefix
   """
}


@Transform("otus")
dereplicate = {
   doc title: "Collapse duplicate reads, record their abundance and sort them by abundance",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://drive5.com/usearch/manual/derep_fulllength.html
   // http://drive5.com/usearch/manual/sortbysize.html
   // usearch -derep_fulllength filtered.fasta -output derep.fasta -sizeout
   // usearch -sortbysize derep.fasta -output derep2.fasta
   exec """
      module load usearch/7.0.1001 &&
      TEMP_FILE=`mktemp tmp_dereplicate_XXXXXXXX.otus` &&
      usearch -derep_fulllength $input.fna -output \$TEMP_FILE -sizeout -threads $threads &&
      usearch -sortbysize \$TEMP_FILE -output $output.otus &&
      rm \$TEMP_FILE
   """
}


rm_chimeras = {
   doc title: "Chimera filtering using a reference database",
       desc:  """Parameters:
                    'db', FASTA file of high-quality, chimera-free sequences,
                    'skip', boolean to skip this step (default: 0)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   if (skip == 1) {
      exec """
         echo "Skipping chimera removal"
      """
      forward input
   } else {
      filter("rm_chimeras") {
         // http://www.drive5.com/usearch/manual/uchime_ref.html
         // usearch -uchime_ref otus1.fa -db $d/gold.fa -strand plus -nonchimeras otus2.fa
         exec """
            echo "Removing chimeras" &&
            module load usearch/7.0.1001 &&   
            usearch -uchime_ref $input.otus -db $db -threads $threads -strand plus -nonchimeras $output.otus
         """
      }
   }
}


@Transform("uc")
otu_mapping = {
   doc title: "Map each reads to an OTU (if possible)",
       desc:  """Parameters:
                    'perc', the minimum identity % required (e.g. 97%)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // The uc2otutab.py script expects reads headers with a tag called barcodelabel=
   // http://drive5.com/usearch/manual/usearch_global.html
   // usearch -usearch_global reads.fa -db otus.fa -strand plus -id 0.97 -uc readmap.uc
   // find how many reads did not match (lines that start with an N, i.e. no-hit records)
   // grep -c "^N" readmap.uc
   def id = perc / 100
   exec """
      module load usearch/7.0.1001 &&
      TEMP_FILE=`mktemp tmp_otu_mapping_XXXXXXXX.fna` &&
      perl -p -e 's/^>(\\S+)_(\\S+)/>\$2;barcodelabel=\$1;/' $input.fna > \$TEMP_FILE &&
      usearch -usearch_global \$TEMP_FILE -db $input.otus -strand plus -id $id -threads $threads -uc $output.uc &&
      rm \$TEMP_FILE &&
      NO_MATCH=`grep -c "^N" $output.uc` &&
      echo "Number of non-matching reads: \$NO_MATCH"
   """
}


@Transform("blast")
usearch_global = {
   doc title: "Assign taxonomy using USEARCH global alignments",
       desc:  """Parameters:
                    'perc', the minimum identity % required (e.g. 97%)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   def id = perc / 100
   exec """
      module load usearch/7.0.1001 &&
      usearch -usearch_global $input.otus -db $input.fna -id $id -threads $threads -strand plus -blast6out $output &&
      if [ ! -s "$output" ]; then
         echo "WARNING: No similarities were found! There might be an issue..." 1>&2;
      fi
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


@Transform("primers")
orient_primers = {
   doc title: "Find out which of the given primers is forward and reverse",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // Check orientation of >1 read to not be affected by sequencing errors
   exec """
      module load bioperl &&
      TEMP_FILE=`mktemp tmp_orient_primers_XXXXXXXX.fna` &&
      extract_first_seqs --input $input.fna --number 10 --output \$TEMP_FILE &&
      PRIMER1=`extract_first_seqs --input $input.primers --number 1 | convert_seq_format --format raw` &&
      echo "Primer 1: \$PRIMER1" &&
      PRIMER2=`extract_last_seqs --input $input.primers --number 1 | convert_seq_format --format raw` &&
      echo "Primer 2: \$PRIMER2" &&
      TEST1=`extract_amplicons -i \$TEMP_FILE -f \$PRIMER1 -m bioperl 2> /dev/null` &&
      TEST2=`extract_amplicons -i \$TEMP_FILE -f \$PRIMER2 -m bioperl 2> /dev/null` &&
      FWD="Primer 1" &&
      if [[ ! -z \$TEST1 ]]; then
         : ;
      elif [[ ! -z \$TEST2 ]]; then
         FWD="Primer 2" &&
         read PRIMER1 PRIMER2 <<<"\${PRIMER2} \${PRIMER1}";
      else
         echo "Error: Primers do not match amplicons" 1>&2 &&
         exit 1;
      fi &&
      echo ${FWD}" has forward orientation" &&
      echo -e ">fwd\\n\${PRIMER1}\\n>rev\\n\${PRIMER2}" > $output.primers &&
      rm \$TEMP_FILE
   """
}


extract_amplicons = {
   doc title: "Given primers, extract amplicons from a file of sequences, trim them and save the results globally",
       desc:  """The first primer should be the primer from which sequencing starts. Parameters:
                    'db', file of database sequences""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bioperl &&
      module load emboss &&
      TRIM_LEN=`extract_first_seqs --input $input.fna --number 1 | get_seq_length | cut -f 2` &&
      echo "Trim length: \$TRIM_LEN" &&
      FWD_PRIMER=`extract_first_seqs --input $input.primers --number 1 | convert_seq_format --format raw` &&
      echo "Fwd primer: \$FWD_PRIMER" &&
      REV_PRIMER=`extract_last_seqs --input $input.primers --number 1 | convert_seq_format --format raw` &&
      echo "Rev primer: \$REV_PRIMER" &&
      TAXO_DB=$db &&
      TRIMMED_TAXO_DB=\${TAXO_DB%.*}_trimmed_\${FWD_PRIMER}_\${TRIM_LEN}_bp.fna &&
      if [ ! -e \$TRIMMED_TAXO_DB ]; then
         echo "Trimming reference sequences in \$TRIMMED_TAXO_DB" &&
         module load emboss &&
         extract_amplicons -i $TAXO_DB -f \$FWD_PRIMER -r \$REV_PRIMER -t \$TRIM_LEN -e 1 -o \$TRIMMED_TAXO_DB;
      else
         echo "Re-using trimmed reference sequences in \$TRIMMED_TAXO_DB";
      fi &&
      ln -f -s \$TRIMMED_TAXO_DB $output.fna
   """
}


@filter("rarefy")
rarefy = {
   doc title: "Rarefy microbial profiles",
       desc:  """Parameters:
                    'num', number of reads to rarefy to (e.g. 1000)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bio-community &&
      bc_rarefy -if $input -nr inf -ss $num -op $output.prefix
   """
}


rm_otus = {
   doc title: "Remove OTUs that match the given names",
       desc:  """Parameters:
                    'str', name of OTUs to remove (e.g. "*bacteria Eukaryota")""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   rm_otus = str.toString() == "true" ? "" : str   // Groovy issue 93
   if (rm_otus) {
      filter("rm_otus") {
         exec """
            module load bio-community &&
            bc_manage_members -if $input -en $rm_otus -op $output.prefix
         """
      }
   } else {
      exec """
         echo "Skipping removal of named OTUs"
      """
      forward input
   }
}


rm_samples = {
   doc title: "Remove samples that match the given names",
       desc:  """Parameters:
                    'str', name of samples to remove (e.g. "*pelagic benthic")""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   rm_samples = str.toString() == "true" ? "" : str   // Groovy issue 93
   if (rm_samples) {
      filter("rm_samples") {
         exec """
            module load bio-community &&
            bc_manage_samples -if $input -en $rm_samples -op $output.prefix
         """
      }
   } else {
      exec """
         echo "Skipping removal of named samples"
      """
      forward input
   }
}


cat_samples = {
   doc title: "Concatenate samples that match the given names",
       desc:  """Parameters:
                    'str', name of samples to concatenate (e.g. "-mn A B -mn C D E")""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   cat_samples = str.toString() == "true" ? "" : str   // Groovy issue 93
   if (cat_samples) {
      filter("cat_samples") {
         exec """
            module load bio-community &&
            bc_manage_samples -if $input -op $output.prefix -rm first $cat_samples
         """
      }
   } else {
      exec """
         echo "Skipping removal of named samples"
      """
      forward input
   }
}


@Filter("copyrighter")
copyrighter = {
   doc title: "Perform GCN bias correction using CopyRighter",
       desc:  """Parameters:
                    'db', gene copy number database (skip if db is '')""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   copyrighter = db.toString() == "true" ? "" : db   // Groovy issue 93
   if (copyrighter) {
      exec """
         echo "Running CopyRighter" &&
         module load copyrighter &&
         copyrighter -i $input -d $copyrighter -o $output -l id -v
      """
   } else {
      exec """
         echo "Skipping CopyRighter... simply converting to relative abundance" &&
         module load bio-community &&
         bc_summarize -if $input -cr 1 -md 0 -rl 0 -op $output.prefix
      """
   }
}


tax_summary = {
   doc title: "Make a taxonomic summary and group taxa <1% in all samples together",
       desc:  """Parameters:
                    'lvl', level of taxonomic summary, e.g. 6 for Greengenes genera""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   produce("${input.prefix}.summary*.*") {
      exec """
         module load bio-community &&
         bc_summarize -if $input -cr 0 -md 0 -rl 1 -tl $lvl -op ${input.prefix}.summary
      """
   }
}
