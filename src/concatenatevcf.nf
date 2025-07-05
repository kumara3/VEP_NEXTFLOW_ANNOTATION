process concatenate_vcfs {
    label "action_vcf"
    cache "lenient"
    errorStrategy "retry"
    maxRetries 2
    cpus = 2
    memory = '4 GB'
    publishDir "concatenated_vcfs/", mode: 'copy'

    input:
    tuple val(chrom), path(vcf_files)

    output:
    //tuple val(chrom), val(name), path("${name}.${chrom}.vep.vcf.gz"), path("${name}.${chrom}.vep.vcf.gz.csi")
    tuple path("${chrom}.vep.vcf.gz"),  path("${chrom}.vep.vcf.gz.tbi")
	
	"""
    for f in ${vcf_files}; do bcftools index \${f};done
    for f in ${vcf_files}; do echo "\${f}"; done |sort -V > vcf_files.txt
    bcftools concat -f vcf_files.txt -O z -o ${chrom}.vep.vcf.gz
    bcftools index -t ${chrom}.vep.vcf.gz
	"""
}

