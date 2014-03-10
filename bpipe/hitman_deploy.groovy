////////////////////////////////////////////////////////////////////////////////

// Clean and trim amplicon reads
// Input FASTQ reads must have Sanger quality scores and be labelled by sample.
// Note: the code for calculating sequence statistics assumes that the sequences
// are NOT wrapped.
// Usage example: bpipe run pipeline_QA.groovy sediments.fastq

// Options:

// Quality score from which to start truncating the 3' end of a read
QUAL_TRUNC=13
// e.g. if you pick 16, all quality scores will be >=17 after truncation
// Q20: 1% , Q19: 1.3%, Q18: 1.6%, Q17: 2%, Q16: 2.5%, Q15: 3.2%, Q14: 4%, Q13: 5%
// Ambiguities such as N typically have Q0. A Q13 threshold will essentially
// eliminate Ns and other low quality bases from the reads.

// Filter out reads with over this % of expected errors (based on quality scores)
EE_PERC=3.0

// Flag to perform or skip the Acacia denoising step
DENOISE=1

// Read trimming/discard length
TRIM_LEN=250
// 250bp captures the whole V6-V7 region using the ACE 926F-1392R primers (V6-V8)

////////////////////////////////////////////////////////////////////////////////


//////////
// TODO
// * Use very clean reads (qual_trunc and ee_filt) to determine OTUs,
//   but try to classify even not so clean reads (ee_filt only) into OTUs.
//   This means implementing ee_filt before qual_trunc, saving READS into variable
//////////


about title: "hitman_deploy - Sequence quality assurance" 


LOG_FILE=""

acacia = {
   doc title: "Denoise FASTQ reads with Acacia"
   // Also filter N before position $TRIM_LEN
   produce(input.prefix+".acacia.fna") {
      if (DENOISE == 1) {
         exec """
            echo "Running Acacia!"                                       &&
            module load acacia/1.52                                      &&
            TRIM_TO_LENGTH=\$(($TRIM_LEN+20))                            &&
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
            echo "FILTER_N_BEFORE_POS=$TRIM_TO_LENGTH"    >> \$CONF_FILE &&
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
      } else {
         exec """
            echo "Skipping Acacia..." &&
            module load fastx_toolkit &&
            fastq_to_fasta -i $input.fastq -n -Q 33 -o $output.fna
         """
      }
   }
}


@Filter("discard_strict")
discard_strict = {
   doc title: "Remove small FASTA reads"
   exec """
      module load fastx_toolkit &&
      fastx_clipper -i $input.fna -a ZZZZ -C -n -l $TRIM_LEN -Q 33 -o $output.fna
   """
}


@Filter("discard_lax")
discard_lax = {
   doc title: "Remove small FASTQ reads"
   exec """
      module load fastx_toolkit &&
      fastx_clipper -i $input.fastq -a ZZZZ -C -n -l \$(($TRIM_LEN-10)) -Q 33 -o $output.fastq
   """
}


@Filter("trim_strict")
trim_strict = {
   doc title: "Quality-trim FASTQ reads"
   exec """
      module load fastx_toolkit &&
      fastx_trimmer -i $input.fna -f 1 -l $TRIM_LEN -Q 33 -o $output.fna
   """
}


@Filter("trim_lax")
trim_lax = {
   doc title: "Quality-trim FASTA reads"
   exec """
      module load fastx_toolkit &&
      fastx_trimmer -i $input.fastq -f 1 -l \$(($TRIM_LEN+10)) -Q 33 -o $output.fastq
   """
}


@Filter("quality_trunc")
quality_trunc = {
   doc title: "Truncate the 3' end of reads when their quality drops too low"
   // http://www.drive5.com/usearch/manual/fastq_filter.html
   exec """
      module load usearch/7.0.1001 &&
      usearch -fastq_filter $input.fastq -fastqout $output.fastq -fastq_truncqual $QUAL_TRUNC
   """
}


@Filter("ee_filter")
ee_filter = {
   doc title: "Remove reads with a percentage of expected errors too large"
   // The calculation of expected errors uses each individual quality score
   // http://www.drive5.com/usearch/manual/fastq_filter.html
   exec """
      module load usearch/7.0.1001 &&
      MAX_EE=`echo "${TRIM_LEN}*${EE_PERC}/100" | bc -l` &&
      usearch -fastq_filter $input.fastq -fastqout $output.fastq -fastq_maxee $MAX_EE
   """
}


@Transform("stats")
seq_stats = {
   doc title: "Record basic FASTQ sequence statistics (length, quality) into a log file"
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
seq_stats_fasta = {
   doc title: "Record basic FASTA sequence statistics (length, quality) into a log file"
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


Bpipe.run {
   seq_stats.using(stage:'initial') +
   discard_lax +
   seq_stats.using(stage:'discard') +
   quality_trunc +
   seq_stats.using(stage:'qtrunc') +
   discard_lax +
   seq_stats.using(stage:'qdiscrd') +
   trim_lax +
   seq_stats.using(stage:'trim') +
   ee_filter +
   seq_stats.using(stage:'eefilt') +
   acacia +
   seq_stats_fasta.using(stage:'acacia') +
   discard_strict +
   seq_stats_fasta.using(stage:'fdiscrd') +
   trim_strict +
   seq_stats_fasta.using(stage:'ftrim')
}
