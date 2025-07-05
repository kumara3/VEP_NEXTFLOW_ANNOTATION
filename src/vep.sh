#!/bin/bash
########################################################################################################
##Usage
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
#--dir_cache /data
#--bind $plugin_dir:/opt/vep/.vep/ \
#apptainer exec /home/s163196/vep_data/ensembl-vep_latest.sif \
#--bind $cache_dir:/data \
#--bind "$(dirname $outprefix)" \
#--plugin DosageSensitivity,file=$DosageSensitivity \
singularity run --bind /scratch/:/scratch/ /home/s163196/vep_data/ensembl-vep_latest.sif \
  vep \
    -i $invcf -o $outprefix.annotated.tab --pick_allele --fork 4 --species $species --assembly "GRCh38" \
    --dir_cache $cache_dir \
    --dir_plugins $plugin_dir --offline \
    --cache --vcf \
    --force_overwrite \
    --warning_file "${outprefix}.warning.txt" \
    --fasta $reffasta --show_ref_allele --nearest gene --symbol --ccds --variant_class --hgvsg --hgvs --shift_hgvs 1 --sift b --polyphen b --domains --regulatory --canonical --check_existing --clin_sig_allele 1 --biotype --gene_phenotype \
    --af --max_af --af_1kg --af_gnomad --af_gnomadg \
    --plugin NMD \
    --plugin REVEL,$REVEL \
    --plugin CADD,snv=$cadd_snv,indel=$cadd_indel \
    --plugin LoFtool,$LofTools \
    --plugin GWAS,file=$gwas \
    --custom file=$conservation_score,short_name=GERP,format=bigwig \
    --plugin dbNSFP,$dbNSFP,GERP++_RS \
    --plugin AlphaMissense,file=$AlphaMissense,transcript_match=1 \
    --plugin SpliceAI,snv=$snv_spliceAI,indel=$indel_spliceAI \
    --plugin PrimateAI,$primateai \
    --plugin PolyPhen_SIFT,db=$polyphen_sift \
    --plugin LoF,loftee_path:$loftee_path,conservation_file:$loftee_sql,human_ancestor_fa:$loftee_ref_ancestor,gerp_bigwig:$gerp_conservation_score \
    --plugin EVE,file=$eve \
    --custom file=$clinvar,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    --custom file=$phyloP100way,short_name=phyloP100way,format=bigwig \
    --custom file=$phastCons100way,short_name=phastCons100way,format=bigwig \
    --custom file=$Gnocchi,short_name=Gnocchi_scores,format=bed \
    #--fields "#Uploaded_variation,Location,REF_ALLELE,Allele,Gene,Feature,Feature type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,VARIANT_CLASS,SYMBOL,NEAREST,GENE_PHENO,CANONICAL,CCDS,DOMAINS,GENE_PHENO,BIOTYPE,HGVSc,HGVSp,HGVSg,CLIN_SIG,SIFT,PolyPhen,phyloP100way,phastCons100way,Gnocchi_scores,GERP++_RS,am_class,am_pathogenicity,SpliceAI_pred,PrimateAI,PolyPhen_humdiv_pred,PolyPhen_humdiv_score,PolyPhen_humvar_pred,PolyPhen_humvar_score,SIFT_pred,SIFT_score,LoF,LoF_filter,LoF_flags,LoF_info,EVE_CLASS,EVE_SCORE,GERP,ClinVar,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,REVEL,CADD_PHRED,CADD_RAW,LoFtool"
    # --plugin GWAS,file=$gwas \ --takes a lot of time when running for first time
    #--plugin LoF,human_ancestor_fa:$human_ancestor_fa,loftee_path:$VEP_PLUGINS,conservation_file:$conservation_file,gerp_file:$gerp_file,skip_lof_info:true
   
    #--fields "#Uploaded_variation,Location,REF_ALLELE,Allele,Gene,Feature,Feature type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,VARIANT_CLASS,SYMBOL,HGVSc,HGVSp,HGVSg,NEAREST,SIFT,PolyPhen,CANONICAL,CCDS,DOMAINS,SV,CLIN_SIG,GENE_PHENO,BIOTYPE"

#tabix -h $invcf "shrad" | 
#--config "vep.ini" 
#--plugin LoF,loftee_path:$loftee,$human_ancestor_fa \
#--dir_cache /data \
# --fasta $reffasta \
#--check_existing
  #--compress_output bgzip \
  #\
########################################################################################################
## convert to tab
# which python2
# python2 $VEP_PLUGINS/src/tableize_vcf.2.py \
# --input $outprefix.vcf.gz \
# --output $outprefix.tsv \
# --info AC,AN,AF,VQSLOD --vep_info ALLELE_NUM,SYMBOL,HGVSc,HGVSp,HGVSg,Consequence,IMPACT,MAX_AF,MAX_AF_POPS,CLIN_SIG,CADD_PHRED,SIFT,PolyPhen,LoF --sample_info GT,DP --hom_samples --samples --original_chr_pos_ref_alt --include_id --existing_allele_only

# ## zip and index
# bgzip --force $outprefix.tsv
# tabix -p vcf --force --skip-lines 1 $outprefix.tsv.gz


# bcftools +split-vep chr21_1_34599028.vcf.annotated.annotated.vcf -f '%CHROM:%POS:%REF %CSQ[\t%GT]\n' -d -A tab
