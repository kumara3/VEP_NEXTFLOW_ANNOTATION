#!/bin/bash
########################################################################################################
##Usage##
Usage="Usage: sh `basename "$0"` (Followed_by_paramters_in_below_order)\n\
1. Inputvcf (/FullPath/to/input.vcf.gz)\n\
2. OutPrefix (/FullPath/perfix)\n\
3. Species\n\
\n\
Example: sh `basename "$0"` /path/to/data/input.vcf.gz /path/to/result/output species\n\
which gives a final output file: /path/to/result/output.vcf.gz\n\n\
"
########################################################################################################

ExpectedArguments=3
if [ $# -ne $ExpectedArguments ]
then
        echo -e $Usage
        echo -e "Number of arguments specifed:\t $#  \nNumber of arguments required:\t $ExpectedArguments"
        echo -e "\nArguments provided:\n$*"
        exit 1
fi


########################################################################################################
set -e
set -o verbose

## set parameters
invcf=$1
outprefix=$2
species=$3
# module load apptainer/1.1.6
# module load java/11.0.4
# module load nextflow/24.04.4.5917
    
# source config file. The config file has path to all the database and plugins 
#source /home/s163196/UTSW-Bioinfo-Lab/Pipeline/tiger/ica-pipelines/src/configFile.txt 
vep \
    -i $invcf -o ${outprefix}_vep.vcf.gz --pick_allele --fork 4 --species $species --assembly "GRCh38" \
    --dir_cache $cache_dir \
    --dir_plugins $plugin_dir --offline \
    --cache --vcf --compress_output bgzip \
    --force_overwrite \
    --warning_file "${outprefix}.warning.txt" \
    --fasta $reffasta --show_ref_allele --symbol --variant_class --canonical --check_existing --biotype \
    --plugin LoFtool,$LofTools \
    --custom file=$clinvar,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG \
    --custom file=$gwas,short_name=GWAS_pmid,format=bed \
    --plugin PrimateAI,$primateai \
    --plugin EVE,file=$eve \
    --plugin CADD,snv=$cadd_snv,indel=$cadd_indel \
    --custom file=$Gnocchi,short_name=Gnocchi_scores,format=bed \
    --custom file=$conservation_score,short_name=GERP_RS,format=bigwig \
    --custom file=$phastCons100way,short_name=phastCons100way,format=bigwig \
    --custom file=$gnomad_genome_flags,short_name=gnomAD_genome_flags,format=bed,type=exact \
    --custom file=$gnomad_exome_flags,short_name=gnomAD_exome_flags,format=bed,type=exact \
    --plugin REVEL,$REVEL \
    --plugin AlphaMissense,file=$AlphaMissense,transcript_match=1 \
    --plugin PolyPhen_SIFT,db=$polyphen_sift \
    --plugin SpliceAI,snv=$snv_spliceAI,indel=$indel_spliceAI \
    --plugin LoF,loftee_path:$loftee_path,conservation_file:$loftee_sql,human_ancestor_fa:$loftee_ref_ancestor,gerp_bigwig:$gerp_conservation_score \
    --af --max_af --af_1kg --af_gnomad --af_gnomadg \
    --fields Existing_variation,VARIANT_CLASS,SYMBOL,BIOTYPE,LoFtool,Consequence,IMPACT,Protein_position,Amino_acids,ClinVar_CLNSIG,GWAS_pmid,CADD_PHRED,PrimateAI,EVE_SCORE,Gnocchi_scores,GERP_RS,phastCons100way,REVEL,am_pathogenicity,PolyPhen_humdiv_score,PolyPhen_humvar_pred,SIFT_score,SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,LoF,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,gnomAD_exome_flags,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_MID_AF,gnomADe_NFE_AF,gnomADe_SAS_AF,gnomAD_genome_flags,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_SAS_AF
#help:  #https://github.com/konradjk/loftee/releases/tag/v1.0.4_GRCh38
#https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/
#https://github.com/konradjk/loftee/issues/111
##https://grch37.ensembl.org/info/docs/tools/vep/script/vep_example.html
#https://grch37.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#revel
#https://community.ukbiobank.ac.uk/hc/en-gb/community/posts/21846172273053-Missing-Loftee-annotations-when-running-with-Hail
#--plugin LoF,human_ancestor_fa:$human_ancestor_fa,loftee_path:$VEP_PLUGINS,conservation_file:$conservation_file,gerp_file:$gerp_file,skip_lof_info:true
#https://grch37.ensembl.org/info/docs/tools/vep/script/vep_custom.html
#https://kitzlerlab.github.io/WES-pipeline/variant-annotation.html
#https://samtools.github.io/bcftools/howtos/plugin.split-vep.html
#https://groups.google.com/a/soe.ucsc.edu/g/genome/c/yAJNS7Yuel8?pli=1
#https://useast.ensembl.org/info/genome/variation/prediction/protein_function.html
#https://discuss.gnomad.broadinstitute.org/t/differences-for-lc-plof-flags-between-v4-1-0-and-v2-1-1/405/2
#https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/
#rsync -Lav rsync://ftp.ebi.ac.uk/ensemblorg/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1
#https://training.nextflow.io/2.1.1/hello_nextflow/05_hello_containers/#234-inspect-how-nextflow-launched-the-containerized-task
#https://github.com/edgano/VEP-nf/blob/main/main.nf
#https://sateeshperi.github.io/nextflow_varcal/nextflow/nextflow_processes#processes
#https://samtools.github.io/bcftools/howtos/plugin.split-vep.html
#/https://sateeshperi.github.io/nextflow_varcal/nextflow/nextflow_operators