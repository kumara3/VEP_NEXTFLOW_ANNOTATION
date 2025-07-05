process check_vcf_split {
    label 'action_vcf'
    cpus = 1
    memory = '1 GB'
    errorStrategy 'retry'
    maxRetries 2
    cache "lenient"
    //publishDir "check_vcf_split/", mode: 'move'
    def user = System.getenv('USER')
    def tmpdir = "/scratch/${user}/tmp"
    def taskdir = "/scratch/${user}/nxf_work"


    input:
    tuple val(shard_chr), path(split_vcf), path(split_vcf_index)

    output:
    path "check_vcf_split"

    script:
    """
    echo "split vcf: $split_vcf"
    echo "vcf index: $split_vcf_index"
    mkdir -p check_vcf_split

    if [ -z "\$(bcftools view -H $split_vcf |head -n 1)"]; then
       mv ${split_vcf} check_vcf_split/
       mv ${split_vcf_index} check_vcf_split/
    else
        echo "Variants found, do not move"
    fi
    """
    
}