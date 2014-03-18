// Shared modules for Bpipe


gunzip = {
   doc title: "Uncompress files compressed with GZIP",
       desc:  """Parameters:
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
            gunzip -c $input.gz > $out
         """
      }
   } else {
      produce(out) {
         exec """
            echo "Skipping GZIP decompression..." &&
            cp $input $output
         """
      }
      // TODO: When bug #90 is fixed, replace by: forward input
      // See code.google.com/p/bpipe/issues/detail?id=90
   }
}


@Transform("mcf")
qmapping2mcf = {
   doc title: "Create a 454 MID mapping file from a QIIME mapping file",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      echo -e "custom_MIDs\n{" > $output.mcf &&
      cut -f 1,2 $input.qmapping | grep -v ^# | perl -p -e 's/^(.+)\t(.+)\$/mid = "\$1", "\$2", 2;/;' >> $output.mcf &&
      echo -e "}\n" >> $output.mcf
   """
}


@Transform("tmapping")
qmapping2tmapping = {
   doc title: "Create a tabular mapping file (2 columns) from a QIIME mapping file",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      cut -f 1-2 $input.qmapping > $output.tmapping
   """
}


split_sff_libraries = {
   doc title: "Use MIDs in a MCF mapping file to separate SFF records belonging to different samples using sfffile",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   produce("${input.prefix}*.sff") {
      exec """
         module load 454 &&
         sfffile -s custom_MIDs -mcf $input.mcf -o $input.prefix $input.sff
      """
   }
}


split_fastq_libraries = {
   doc title: "Use MIDs in a tabular mapping file to separate FASTQ records belonging to different samples using fastq-multx",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   //http://code.google.com/p/ea-utils/wiki/FastqMultx
   produce("${input.prefix}.*.fastq") {
      exec """
         module load ea_utils &&
         fastq-multx -B $input.tmapping $input.fastq -o ${input.prefix}.%.fastq -b -m 0 &&
         rm ${input.prefix}.unmatched.fastq
      """
   }
   // Alternatives:
   // - QIIME split_libraries_fastq.py, but does a bunch of extra stuff that's not needed
   //   https://github.com/qiime/qiime/blob/master/scripts/split_libraries_fastq.py
   // - FASTX_toolkit fastx_barcode_splitter, but does not support different MID lengths.
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
         echo "Skipping conversion from SFF to FASTQ..."
      """
      forward input
   }
}


merge_pairs = {
   doc title: "Merge pairs of FASTQ reads that overlap using Pandaseq",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   println("Inputs are: "+inputs)
   println("Number of inputs: "+inputs.size())
   if (inputs.size() != 2) {
      exec """
         echo "Skipping paired-end merging..."
      """
      forward input
   } else {
      filter("merge_pairs") {
         exec """
            echo "Merging paired-end reads"
            module load pandaseq &&
            SINGLES_FILE=`mktemp tmp_pandaseq_singles_XXXXXXXX.txt` &&
            STATS_FILE=`mktemp tmp_pandaseq_stats_XXXXXXXX.txt` &&
            pandaseq -f $input1.fastq -r $input2.fastq -T $threads -F -w $output -g $STATS_FILE -u $SINGLES_FILE
         """
         // &&
         //   NO_MATCH=`grep -c "^N" $output.uc` &&
         //   echo "Number of non-matching reads: \$NO_MATCH"
      }
   }



   // Alternatives:
   //   USEARCH fastq_mergepairs: http://www.drive5.com/usearch/manual/fastq_mergepairs.html
   //   ea-utils fastq-join: : http://code.google.com/p/ea-utils/wiki/FastqJoin
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
   //forward inputs
}


acacia = {
   doc title: "Denoise FASTQ reads with Acacia",
       desc:  """Parameters:
                    'skip', boolean to skip this step (default: 0)""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   produce(input.prefix+".acacia.fna") {
      if (skip==1) {
         exec """
            echo "Skipping Acacia..." &&
            module load fastx_toolkit &&
            fastq_to_fasta -i $input.fastq -n -Q 33 -o $output.fna
         """
      } else {
         exec """
            echo "Running Acacia!"                                       &&
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


@Filter("rm_small_fasta")
rm_small_fasta = {
   doc title: "Remove small FASTA reads",
       desc:  """Parameters:
                    'length', the minimum length to keep""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load fastx_toolkit &&
      fastx_clipper -i $input.fna -a ZZZZ -C -n -l $length -Q 33 -o $output.fna
   """
}


@Filter("rm_small_fastq")
rm_small_fastq = {
   doc title: "Remove small FASTQ reads",
       desc:  """Parameters:
                    'length', the minimum length to keep""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load fastx_toolkit &&
      fastx_clipper -i $input.fastq -a ZZZZ -C -n -l $length -Q 33 -o $output.fastq
   """
}


@Filter("trim_fasta")
trim_fasta = {
   doc title: "Trim FASTA reads",
       desc:  """Parameters:
                    'length', the minimum length to keep""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load fastx_toolkit &&
      fastx_trimmer -i $input.fna -f 1 -l $length -Q 33 -o $output.fna
   """
}


@Filter("trim_fastq")
trim_fastq = {
   doc title: "Trim FASTQ reads",
       desc:  """Parameters:
                    'length', the minimum length to keep""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load fastx_toolkit &&
      fastx_trimmer -i $input.fastq -f 1 -l $length -Q 33 -o $output.fastq
   """
}


@Filter("qual_trim_fastq")
qual_trim_fastq = {
   doc title: "Quality-based 3' end sequence trimming",
       desc:  """Parameters:
                    'qual', the quality score at which to start 3' end trimming""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.drive5.com/usearch/manual/fastq_filter.html
   exec """
      module load usearch/7.0.1001 &&
      usearch -fastq_filter $input.fastq -fastqout $output.fastq -fastq_truncqual $qual
   """
}


@Filter("ee_filter")
ee_filter = {
   doc title: "Remove reads with too many expected errors",
       desc:  """Parameters:
                    'ee', the maximum number of errors""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // The calculation of expected errors uses each individual quality score
   // http://www.drive5.com/usearch/manual/fastq_filter.html
   exec """
      module load usearch/7.0.1001 &&
      usearch -fastq_filter $input.fastq -fastqout $output.fastq -fastq_maxee $ee
   """
}


LOG_FILE=""

@Transform("stats")
fastq_stats = {
   doc title: "Log basic FASTQ sequence statistics",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.drive5.com/usearch/manual/fastq_stats.html
   if(LOG_FILE == "") {
      LOG_FILE=output.stats
   }
   produce(LOG_FILE) {
      exec """
         if [ ! -e $LOG_FILE ]; then
            echo "Creating log file $LOG_FILE" &&
            echo -e "stage\\tnum_seq\\tmin_L\\tmed_L\\tmax_L\\tmin_Q\\tmed_Q\\tmax_Q\\tnum_N" > $LOG_FILE;
         fi &&
         echo "Appending stage $stage to log file $LOG_FILE" &&
         module load usearch/7.0.1001 &&
         TEMP_FILE=`mktemp tmp_stats_${stage}_XXXXXXXX.txt` &&
         usearch -fastq_stats $input.fastq -log \$TEMP_FILE &&
         NUM_SEQS=`grep Recs \$TEMP_FILE | tail -n 1 | perl -n -e 'print((split /\\s+/)[1])'` &&
         MAX_LEN=`grep '^>=' \$TEMP_FILE | perl -n -e 'print((split /\\s+/)[1]);'` &&
         MED_LEN=`grep '[5-9].\\..%\$' \$TEMP_FILE | perl -n -e '@arr=split/\\s+/; if (scalar @arr == 5) { print \$arr[1]; last; }'` &&
         MIN_LEN=`grep '100\\.0%\$' \$TEMP_FILE | perl -e '\$len="?"; while (<>) { @arr=split /\\s+/; \$len=\$arr[1] if scalar(@arr)==5; }; print \$len'` &&
         MAX_Q=`grep -A 2 ^ASCII \$TEMP_FILE | tail -n 1 | perl -n -e 'print((split /\\s+/)[2])'` &&
         MED_Q=`grep '[5-9].\\..%\$' \$TEMP_FILE | perl -e '\$len=""; while (<>) { @arr=split /\\s+/; if (scalar(@arr)==7 && not(\$len)) { \$len=\$arr[2]; last; } }; print \$len;'` &&
         MIN_Q=`grep '100\\.0%\$' \$TEMP_FILE | tail -n 1 | perl -n -e 'print((split /\\s+/)[2])'` &&
         NUM_N=`perl -n -e 'print \$_ if (\$.%4 == 2);' $input.fastq | grep -o N | wc -l` &&
         echo -e "$stage\\t\$NUM_SEQS\\t\$MIN_LEN\\t\$MED_LEN\\t\$MAX_LEN\\t\$MIN_Q\\t\$MED_Q\\t$MAX_Q\\t$NUM_N" >> $LOG_FILE &&
         rm \$TEMP_FILE
      """
   }
   forward input.fastq
}


@Transform("stats")
fasta_stats = {
   doc title: "Log basic FASTA sequence statistics",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   // http://www.drive5.com/usearch/manual/fastq_stats.html
   if(LOG_FILE == "") {
      LOG_FILE=output.stats
   }
   produce(LOG_FILE) {
      exec """
         if [ ! -e $LOG_FILE ]; then
            echo "Creating log file $LOG_FILE" &&
            echo -e "stage\\tnum_seq\\tmin_L\\tmed_L\\tmax_L\\tmin_Q\\tmed_Q\\tmax_Q\\tnum_N" > $LOG_FILE;
         fi &&
         echo "Appending stage $stage to log file $LOG_FILE" &&
         fake_qual $input.fna &&
         TEMP_QUAL=${input.fna}.qual &&
         TEMP_FASTQ=`mktemp tmp_stats_${stage}_XXXXXXXX.fastq` &&
         module load usearch/7.0.1001 &&
         faqual2fastq.py $input.fna \$TEMP_QUAL > \$TEMP_FASTQ &&
         rm \$TEMP_QUAL &&
         TEMP_FILE=`mktemp tmp_stats_${stage}_XXXXXXXX.txt` &&
         usearch -fastq_stats \$TEMP_FASTQ -log \$TEMP_FILE &&
         rm \$TEMP_FASTQ &&
         NUM_SEQS=`grep Recs \$TEMP_FILE | tail -n 1 | perl -n -e 'print((split /\\s+/)[1])'` &&
         MAX_LEN=`grep '^>=' \$TEMP_FILE | perl -n -e 'print((split /\\s+/)[1]);'` &&
         MED_LEN=`grep '[5-9].\\..%\$' \$TEMP_FILE | perl -n -e '@arr=split/\\s+/; if (scalar @arr == 6) { print \$arr[1]; last; }'` &&
         MIN_LEN=`grep '100\\.0%\$' \$TEMP_FILE | perl -e '\$len="?"; while (<>) { @arr=split /\\s+/; \$len=\$arr[1] if scalar(@arr)==5; }; print \$len'` &&
         MAX_Q='-' &&
         MED_Q='-' &&
         MIN_Q='-' &&
         NUM_N=`grep -v '^>' $input.fna | grep -o N | wc -l` &&
         echo -e "$stage\\t\$NUM_SEQS\\t\$MIN_LEN\\t\$MED_LEN\\t\$MAX_LEN\\t\$MIN_Q\\t\$MED_Q\\t$MAX_Q\\t$NUM_N" >> $LOG_FILE &&
         rm \$TEMP_FILE
      """
   }
   forward input.fna
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
         echo "Skipping removal of \$NOF_SINGLETONS singletons..."
      """
      forward input
   } else {
      filter("rm_singletons") {
         // http://drive5.com/usearch/manual/sortbysize.html
         // usearch -sortbysize derep.fasta -output derep2.fasta -minsize 2
         exec """
            NOF_SINGLETONS=`grep '^>' $input.otus | grep -c 'size=1;'` &&
            echo "Removing \$NOF_SINGLETONS singletons!" &&
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
         echo "Skipping chimera removal..."
      """
      forward input
   } else {
      filter("rm_chimeras") {
         // http://www.drive5.com/usearch/manual/uchime_ref.html
         // usearch -uchime_ref otus1.fa -db $d/gold.fa -strand plus -nonchimeras otus2.fa
         exec """
            echo "Removing chimeras!" &&
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
      usearch -usearch_global $input.otus -db $input.fna -id $id -threads $threads -strand plus -blast6out $output
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


extract_amplicons = {
   doc title: "Given primers, extract amplicons from a file of sequences, trim them and save the results globally",
       desc:  """The first primer should be the primer from which sequencing starts. Parameters:
                    'db', file of database sequences""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load bioperl &&
      TRIM_LEN=`extract_first_seqs --input $input.fna --number 1 | get_seq_length | cut -f 2` &&
      echo "Trim length: \$TRIM_LEN" &&
      FWD_PRIMER=`extract_first_seqs --input $input2.fna --number 1 | convert_seq_format --format raw` &&
      echo "Fwd primer: \$FWD_PRIMER" &&
      REV_PRIMER=`extract_last_seqs --input $input2.fna --number 1 | convert_seq_format --format raw` &&
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


@Filter("copyrighter")
copyrighter = {
   doc title: "Perform GCN bias correction using CopyRighter",
       desc:  """Parameters:
                    'db', gene copy number database (skip if db is '')""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   if (db) {
      exec """
         echo "Running CopyRighter!" &&
         module load copyrighter &&
         copyrighter -i $input -d $db -o $output -l id -v
      """
   } else {
      exec """
         echo "Skipping CopyRighter... simply converting to relative abundance" &&
         module load bio-community &&
         bc_summarize -if $input -cr 1 -md 0 -rl 0 -op $output.prefix
      """
   }
}


