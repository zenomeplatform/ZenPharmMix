#!/usr/bin/env nextflow

//haps1 = Channel.fromFilePairs(params.in_pgvcf, type: 'file') { file -> file.name.replaceAll(/.vcf.gz|.tbi$/,'') }

haps1 = Channel.fromPath(params.in_hap1, type: 'file')
haps2 = Channel.fromPath(params.in_hap2, type: 'file') 

//haps2 = Channel.fromFilePairs(params.in_haplo, type: 'file') { file -> file.name.replaceAll(/.vcf.gz|.tbi$/,'') }

haps1
    .combine(haps2)
    .set {data}
//    .println()


process make_diplo {
//   maxForks 10
//   errorStrategy 'ignore'
//   tag "${name}" 
//   label 'big_mem'
   publishDir params.out_dir, mode: 'copy', overwrite: 'true'

   input:
      file(vcfs) from data   

   output:	        
      file("*_diplo.vcf.gz") into diplotypes_ch
     
   script:
    """ 
    tabix ${vcfs.get(0)}
    tabix ${vcfs.get(1)}
    bcftools merge ${vcfs.get(0)} ${vcfs.get(1)} | bgzip -c > ${vcfs.get(0).name.replace(/.vcf.gz/, "")}_${vcfs.get(1).name.replace(/.vcf.gz/, "")}_diplo.vcf.gz
    """

}