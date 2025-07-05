process splitfile {
    label 'action_vcf'
    cpu = 1
    memory = '2 GB'
    errorStrategy 'retry'
    maxRetries 2
    cache "lenient"
    publishDir "vcf_splits/", mode: 'copy'

    input:
    tuple path(vcf_file), path(vcf_index), val(shard_chr)
    

    output:
    tuple val(shard_chr),path("*.vcf.gz"), path("*.vcf.gz.tbi")

    script:
    def filename = shard_chr.replaceAll("[:\\-]", "_")
    """
    
    bcftools view -r $shard_chr -O z -o ${filename}.vcf.gz $vcf_file
    bcftools index -t ${filename}.vcf.gz
    """
}
