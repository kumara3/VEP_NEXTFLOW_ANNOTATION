process convertToTab {
    label 'action_vcf'
    cpu = 1
    memory = '2 GB'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "annotations_per_chr/", mode: 'copy'

    input:
    tuple path(vcf_file), path(vcf_index)
    

    output:
    path "*.tab.gz"

    script:
    def filename = vcf_file.getSimpleName()
    //.replaceAll(".vcf.gz", "")
    """
    echo "${filename}"
    echo -e "CHROM\tPOS\tREF\tALT\tExisting_variation\tVARIANT_CLASS\tSYMBOL\tBIOTYPE\tLoFtool\tConsequence\tIMPACT\tProtein_position\tAmino_acids\tClinVar_CLNSIG\tGWAS_pmid\tCADD_PHRED\tPrimateAI\tEVE_SCORE\tGnocchi_scores\tGERP++_RS\tphastCons100way\tREVEL\tam_pathogenicity\tPolyPhen_humdiv_score\tPolyPhen_humvar_pred\tSIFT_score\tSpliceAI_pred_DP_AG\tSpliceAI_pred_DP_AL\tSpliceAI_pred_DP_DG\tSpliceAI_pred_DP_DL\tLoF\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tgnomAD_exome_flags\tgnomADe_AF\tgnomADe_AFR_AF\tgnomADe_AMR_AF\tgnomADe_ASJ_AF\tgnomADe_EAS_AF\tgnomADe_FIN_AF\tgnomADe_MID_AF\tgnomADe_NFE_AF\tgnomADe_SAS_AF\tgnomAD_genome_flags\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_SAS_AF\t\$(bcftools query -l $vcf_file |sed 's/\$/\\.GT/g'|tr '\n' '\t')" |sed 's/\\t\$//g' > header.txt
    cat header.txt <(bcftools +split-vep $vcf_file -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%GT]\n' -d -A tab) |bgzip -c > ${filename}.vep.tab.gz
    
    """
}