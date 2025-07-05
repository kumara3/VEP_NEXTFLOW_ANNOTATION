process runvep {
    label 'vep'
    cpus = 4
    memory = '16 GB'
    errorStrategy 'retry'
    maxRetries 2
    cache "lenient"
    publishDir "vep_annotation_vcfs/", mode: 'copy'
    def user = System.getenv('USER')
    def tmpdir = "/scratch/${user}/tmp"
    def taskdir = "/scratch/${user}/nxf_work"

    env = [
    'APPTAINERENV_TMPDIR'        : "/scratch/${user}/tmp",
    'APPTAINERENV_NXF_TASK_WORKDIR': "/scratch/${user}/nxf_work"
    ]


    input:
    tuple val(shard_chr), path(split_vcf), path(split_vcf_index), path(vep_script)
    val species
    path bind_path
    //path config_file


    output:
    tuple  val("${split_vcf.simpleName.tokenize('_')[0]}"), path("${split_vcf.simpleName}_vep.vcf.gz")

    script:
    def split_vcf_outprefix = "${split_vcf.simpleName}"
    
    """
    echo "vep script: $vep_script"
    echo "split vcf: $split_vcf"
    echo "vcf index: $split_vcf_index"
    echo "$split_vcf_outprefix"
    source "/home/s163196/UTSW-Bioinfo-Lab/Pipeline/tiger/ica-pipelines/src/configFile.txt"
    sh $vep_script "$split_vcf" ${split_vcf_outprefix} $species
    """
    
}
// helps
// https://community.seqera.io/t/getting-the-file-basename-depending-on-the-file-extension/622/3
// https://github.com/CERC-Genomic-Medicine/vep_pipeline/blob/master/Annotation.nf
// https://github.com/CERC-Genomic-Medicine/vep_pipeline
// #
// vep \
//     -i $split_vcf -o ${split_vcf_outprefix}.annotated.tab --pick_allele --species $species --fork 4 --assembly "GRCh38" \
//     --dir_cache $cache_dir \ 
//     --offline \
//     --cache --tab \
//     --force_overwrite \
//     --warning_file "${outprefix}.warning.txt" \
//     --fasta $reffasta --hgvsg --hgvs --shift_hgvs 1 --sift b --polyphen b --domains --regulatory --canonical \
//     --af --max_af --af_1kg --af_gnomad --af_gnomadg \
//     --plugin NearestGene \
//     --plugin TSSDistance \
//     --plugin REVEL,$REVEL \
//     --plugin CADD,snv=$cadd_snv,indel=$cadd_indel \
//     --custom file=$conservation_score,short_name=GERP,format=bigwig \
//     --plugin dbNSFP,$dbNSFP,LRT_score,GERP++_RS,FATHM_score,MutationTaster_score \
//     --plugin AlphaMissense,file=$AlphaMissense,transcript_match=1 \
//     --plugin SpliceAI,snv=$snv_spliceAI,indel=$indel_spliceAI \
//     --plugin DosageSensitivity,file=$DosageSensitivity \
//     --custom file=$clinvar,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
//     --custom file=$phyloP100way,short_name=phyloP100way,format=bigwig \
//     --custom file=$phastCons100way,short_name=phastCons100way,format=bigwig
   
// source $config_file
//     echo "$config_file"