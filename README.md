# VEP variant annotation pipeline

## Prerequisites

The following software is required:
- Apptainer (tested with version 1.1.6)
- Nextflow (tested with version 23.04.0)
- Python 3 
- java (tested with version java/11.04)

## Workflow
> [!CAUTION]
> - Your input VCFs must be indexed and have corresponding `.tbi` files
> - Input VCFs will be split by shards per chromosome e.g., chr21:1-34599028;chr21:34599029-46709983;chr22:1-34599028;chr22:34599029-50818468.

## 1. Installation
This section describes how to set up VEP, configure all the necessary files.

### 1.1. Setting up VEP

1. The running required following modules: `apptainer`, `nextflow`, `java` modules:
   ```
   module load apptainer/1.1.6 java/11.0.4 nextflow/23.04.0.5857
   ```
2. Build SIF images with additional tools: `vep from ensemble`,`samtools`, `bcftools`, `tabix`.

   The SIF images have already been built for you. 
   
   ```
   vep: ~/vep_data/ensembl-vep_latest.sif
   bcftools: ~/vep_data/bcftools_1.13--h3a49de5_0.sif
   tabix: ~/vep_data/tabix_v1.9-11-deb_cv1.sif
   ```
### 1.2. Downloading data source files
   
1. Download VEP cache files. Theses files have already been downloaded for you. 
   ```
   cache_dir="/reference/Homo_sapiens_vep/data/cache"
   plugin_dir="/reference/Homo_sapiens_vep/plugins/pm"
   ```
2. Loftee - More detailed instructions on how to set up LoFtee are [here](https://github.com/konradjk/loftee).
   ```
   loftee_path="/reference/Homo_sapiens_vep/plugins/pm"
   loftee_sql="/reference/Homo_sapiens_vep/plugins/loftee/loftee.sql"
   loftee_ref_ancestor="/reference/Homo_sapiens_vep/plugins/loftee/human_ancestor.fa.gz"
   ```
3. Download all necessary databases based on human genome build GRCh38  (https://github.com/konradjk/loftee) These should include: GERP conservation scores (only for GRCh38), human_ancestor.fa files, SQL databases with PhyloCSF metrics (SQL files must be unzipped).

4. The custom datasource file have been configured.

5. Make sure to follow the exact directory struture
   ```
   Run the command : rsync -av /scratch/s163196/reference/Homo_sapiens_vep /scratch/<username>/reference/ 
   ├── bind.txt
   ├── data
   │   └── cache
   ├── human_ancestor_fa
   │   ├── gerp_conservation_scores.homo_sapiens.GRCh38.bw
   │   ├── human_ancestor.fa.gz
   │   ├── human_ancestor.fa.gz.fai
   │   └── loftee.sql.gz
   └── plugins
      ├── AlphaMissense
      ├── bench-console_1.0.tar
      ├── CADD
      ├── ClinVar
      ├── dbNSFP4
      ├── DosageSensitivity
      ├── EVE
      ├── fasta
      ├── GERP
      ├── Gnocchi
      ├── gnomAD
      ├── GWAS
      ├── loftee
      ├── LoFTools
      ├── phastCons100way
      ├── phyloP100way
      ├── pm
      ├── Polyphen_sift
      ├── PrimateAI
      ├── REVEL
      ├── spliceAI
      └── StructuralVariantOverlap
   ```

### 1.3. How to execute pipeline

   ```
   nextflow main.nf --help:
   Usage:
   options:
   ==========================================
     V E P  P I P E L I N E : Annotates varaint called VCF file using vep plugins (https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic)

     Usage:
    -------------------------------------------
     --inputvcf          : Input VCF file
     --shardfile         : Genomic coordinates on which to run VEP        
     --help              : See the usage
   ===========================================

   nextflow main.nf --inputvcf <test.vcf.gz>
   ```
  
### 1.4. If you want to add custom VCFs

In this section, I will explain how to integrate custom VCF files,  into your VEP command line/config file using this pipeline.

1. Download the custom data file. The expected format is gff, gtf,vcf and bde file. For example, if you have to include data source "Gnocchi". Download the bed file from : https://www.nature.com/articles/s41586-023-06045-0. 

2. Include the bed file in config file as shown below

 ```
   export Gnocchi="/reference/Homo_sapiens_vep/plugins/Gnocchi/Supplementary_Data_2.bed.gz"
 ```

3. To integrate the custom VCF into your VEP command, add the relevant `--custom` flag to vep.v1.sh. Below is an example configuration:
```
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

```
