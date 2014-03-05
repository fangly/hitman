////////////////////////////////////////////////////////////////////////////////

// Given SFF and QIIME mapping files, extract the desired samples into a single
// FASTQ file.
// Usage example: bpipe run pipeline_prep.groovy Gasket67.sff Gasket67.mapping Gasket76.sff Gasket76.mapping

////////////////////////////////////////////////////////////////////////////////


about title: "16S analysis pipeline - Input preparation" 


@Transform("mcf")
mapping2mcf = {
   doc title: "Create a 454 MID config file from a QIIME mapping file"
   exec """
      echo -e "custom_MIDs\n{" > $output.mcf && 
      cut -f 1,2 $input.mapping | grep -v ^# | perl -p -e 's/^(.+)\t(.+)\$/mid = "\$1", "\$2", 2;/;' >> $output.mcf &&
      echo -e "}\n" >> $output.mcf
   """
}


split_libraries = {
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
   // Alternatively, use FQTRIM or `usearch -fastq_filter seqs.fna ‑fastqout seqs_new.fna -relabel sample_`
}


//sff2fastq_multi = {
//   doc title: "Convert SFF files to the FASTQ format"
//   // Note: sequences are not wrapped
//   println("input: "+input)
//   println("input.prefix: "+input.prefix)
//   println("inputs: "+inputs)
//   produce("${input.prefix}*.fastq") {
//      exec """
//         module load parallel &&
//         module load sff2fastq &&
//         parallel "sff2fastq -o {.}.fastq {}" ::: $inputs
//      """
//   }
//   println("output: "+output)
//}


//@Filter("rename_seqs_multi")
//rename_seqs_multi = {
//   doc title: "Prefix sequence names based on their filename"
//   // From a file like "Gasket76.split_libraries.3833.fastq", use the "3833" prefix
//   // A sequence called "HN758CG01AIWK2" would thus become "3833_HN758CG01AIWK2"
//   // In comparison, QIIME would call it "3833_1"
//   // Note: rely on the sequences not being wrapped
//   exec """
//      module load parallel &&
//      ls ${input.prefix}*.fastq | rev | cut -d '.' -f 2 | rev | parallel "
//         cat ${input.prefix}*{}.fastq | paste - - | perl -p -e 's/^@/\\@{}_/; s/^\\+\\S*/+/' | tr '\\t' '\\n' > ${input.prefix}*{}_renamed.fastq
//      "
//   """
//}


@Filter("cat_files")
cat_files = {
   doc title: "Concatenate all files into a single one"
   exec """
      cat $inputs > $output
   """
}


report = {
   println("The outputs are : "+inputs)
   forward inputs
}


Bpipe.run {
   //"%.*" * [ mapping2mcf + split_libraries + sff2fastq_multi + rename_seqs_multi ]
   "%.*" * [ mapping2mcf + split_libraries +
             "%.sff" * [ sff2fastq + rename_seqs ] ] + cat_files + report
}
