process filter_non_empty_vcfs {
    label 'action_vcf'
    cpus = 1
    memory = '1 GB'
    errorStrategy 'ignore'
    maxRetries 2
    cache "lenient"
    publishDir "check_vcf_splits/", mode: 'copy'
    
    input:
    tuple val(chr), path(vcf), path(vcf_idx)

    output:
    tuple val(chr), path(vcf), path(vcf_idx), path("valid_vcf_list.txt"),path("empty_vcf_list.txt"), path("status.txt")

    script:
    """
    if [ -n "\$(bcftools view -H $vcf | head -n 1)" ]; then
        echo "VCF $vcf has variants."
        echo "$vcf" >> valid_vcf_list.txt
        touch empty_vcf_list.txt # required to satisfy declared outputs
        echo "true" > status.txt # status of vcf file, vcf file is fully content   
    else
        echo "VCF $vcf is empty."
        echo "$vcf" >> empty_vcf_list.txt
        touch valid_vcf_list.txt
        echo "false" > status.txt # status of vcf file, vcf file is empty
    fi
    """
}