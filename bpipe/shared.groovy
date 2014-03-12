@Transform("mcf")
mapping2mcf = {
   doc title: "Create a 454 MID config file from a QIIME mapping file",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      echo -e "custom_MIDs\n{" > $output.mcf &&
      cut -f 1,2 $input.mapping | grep -v ^# | perl -p -e 's/^(.+)\t(.+)\$/mid = "\$1", "\$2", 2;/;' >> $output.mcf &&
      echo -e "}\n" >> $output.mcf
   """
}


split_sff_libraries = {
   doc title: "Use MIDs to separate SFF records belonging to different samples",
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


@Transform("fastq")
sff2fastq = {
   doc title: "Convert an SFF file to FASTQ",
       desc:  """Parameters:
                    none""",
       constraints: "",
       author: "Florent Angly (florent.angly@gmail.com)"
   exec """
      module load sff2fastq &&
      sff2fastq -o $output.fastq $input.sff
   """
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
      PREFIX=`echo $input.fastq | perl -n -e 'print((split /\\./)[-2])'` &&
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
   println("The outputs are : "+inputs)
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

