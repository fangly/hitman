@Transform("mcf")
mapping2mcf = {
   doc title: "Create a 454 MID config file from a QIIME mapping file"
   exec """
      echo -e "custom_MIDs\n{" > $output.mcf &&
      cut -f 1,2 $input.mapping | grep -v ^# | perl -p -e 's/^(.+)\t(.+)\$/mid = "\$1", "\$2", 2;/;' >> $output.mcf &&
      echo -e "}\n" >> $output.mcf
   """
}


split_sff_libraries = {
   doc title: "Use MIDs to separate SFF records from different samples"
   produce("${input.prefix}*.sff") {
      exec """
         module load 454 &&
         sfffile -s custom_MIDs -mcf $input.mcf -o $input.prefix $input.sff
      """
   }
}


@Transform("fastq")
sff2fastq = {
   doc title: "Convert an SFF file to FASTQ"
   exec """
      module load sff2fastq &&
      sff2fastq -o $output.fastq $input.sff
   """
}


@Filter("rename_seqs")
rename_seqs = {
   doc title: "Prefix FASTQ sequence names based on their filename"
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
   doc title: "Concatenate all files into a single one"
   exec """
      cat $inputs > $output
   """
}


final_report = {
   println("The outputs are : "+inputs)
   //forward inputs
}


