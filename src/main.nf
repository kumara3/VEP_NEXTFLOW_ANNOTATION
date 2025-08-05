nextflow.enable.dsl=2
// Execution environment setup

params.help = ""

if(params.help){

    Usage()
    exit(0)
}


def Usage(){
    println """

    Usage:
    options:
     ===========================================
         V E P  P I P E L I N E : Annotates varaint called VCF file using vep plugins (https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic)
  
         Usage:
        -------------------------------------------
         --inputvcf          : Input VCF file
         --shardfile         : shard file (genomic coordinates). This is the file with chr name and chopped chr coordinates
                               input options:  shards.1_M.txt (chr1-chr25); shard.test.txt (chr 21, 22)       
         --help              : See the usage
        ===========================================
         """
         .stripIndent()

}

log.info """\
         ${params.manifest.name} v${params.manifest.version}
         ==========================
         input from   : ${params.inputvcf}
         output to    : ${workflow.launchDir}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         container    : ${workflow.container}
         """
         .stripIndent()

/* check the correctness of file (datatype of fastq files and then datatype of config files)*/

def expectedParams  = [
    'inputvcf'         : String
]

include { splitfile } from './split.nf'
include { runvep } from './vep.nf'
include { concatenate_vcfs } from './concatenatevcf.nf'
include { convertToTab } from './convertTab.nf'
include { filter_non_empty_vcfs } from './filter_non_empty_vcfs.nf'




if(!params.inputvcf) {
    error "Please provide the input vcf file"
}

if(!params.shardfile){
    error "shard file is missing. You need to be more vigilant"
}

params.shard_file = "${workflow.projectDir}/$params.shardfile"
//shards.1_M.txt
//shard.test.txt
//shards.chr2.txt
if(!params.config_file){
    params.config_file = "${workflow.projectDir}/configFile.sh"

}
workflow {
    log.info "Validating config data type"
    def errors = []
    expectedParams.each { key, expectedType ->
        def value = params.get(key)
        if(value == null){
            errors << "Missing parameter key :$key"
        }
        else if(!(value in expectedType)){

            errors << "Invalid type for params.${key}: Expected ${expectedType.simpleName}, Found ${value.getClass().simpleName}"
        }
    }


    if (errors) {
        log.error "Parameter validation failed:"
        errors.each { log.error "- ${it}" }
        System.exit(1)
    } else {
        log.info "All parameters types are valid."
    }
  

    vcf_ch = Channel.fromPath(params.inputvcf)
    vcf_index_ch = Channel.fromPath(params.inputvcf.replaceAll('.vcf.gz$','.vcf.gz.tbi'))
    bind_path = file(params.bind_path)
    //out_pre = Channel.from(params.outputprefix)
    shard_ch = Channel.fromPath(params.shard_file)
    each_chr = shard_ch.splitText().map{ it.trim() }
    // each_chr.view{ item -> "${item.getClass()}" } // print data type
    
    vep_script_ch = Channel.fromPath("${workflow.projectDir}/vep.v1.sh")
    //vcf_ch.view{ file -> "VCF file: ${file}" } : print the output
    
    // split files
    split_input = vcf_ch.combine(vcf_index_ch).combine(each_chr) 
    //split_input.view{ item -> "${item.getClass()}" } : check the data type
    //split_input.view{ tuple -> "${tuple[0]},${tuple[1]},${tuple[2]}" } : print split_input channel contents
    
    split_out = splitfile(split_input) // function call

    split_out.view{ it -> "created VCF:${it[1]}, Index:${it[2]},for chromosome region ${it[0]}" }  // print output of splitfile func.
    //split_out.view{ it.getClass() } // prints the data type of each item: arrayList
    
    split_out.map{ it -> tuple(*it) }.set{ split_out_tuple }
    // split_out_tuple.view{ it.getClass() }

    vcf_ch = filter_non_empty_vcfs(split_out_tuple)
    vcf_ch.view()
    //vcf_ch.view{ it.getClass() } // prints array list
    
    
    vcf_ch.map{ fields -> 
                    def (chr, vcf, vcf_idx, valid_log, empty_log, status_file) = fields 
                    def status = status_file.text.trim() == 'true'
                    return tuple(status,chr,vcf,vcf_idx,valid_log, empty_log) }
                    .branch {
                        success: it[0] == true
                        failed: it[0] == false
                    }.set { results }

    valid_vcf = results.success.map{ it -> tuple(it[1], it[2], it[3]) }
    vep_input = valid_vcf.combine(vep_script_ch)
    //vep_input.view()
    
    // summary of success and failed vcf files
    results.success.map{ it[2].toString() + '\n' } //vcf
        .collectFile(name: 'vep_success_vcfs.txt', storeDir: './report_logs')

    results.failed.map{ it[2].toString() + '\n' } //vcf
        .collectFile(name: 'vep_failed_vcfs.txt', storeDir: './report_logs')



//     vep_input = valid_vcf_splits.combine(vep_script_ch)
//     //vep_input = split_out_tuple.combine(vep_script_ch)
    vep_results = runvep(vep_input,params.species,params.bind_path) // run vep on each split file parallely
    
    //vep_results.view{ it.getClass() }
    //vep_results.view()
    
    // group vcf by chromosome
    group_vcf_by_chr = vep_results.groupTuple()
    group_vcf_by_chr.view{ "Grouped by chromosome: $it" }

    // concatenate vcf files per chromosome
    concatenated_vcfs = concatenate_vcfs(group_vcf_by_chr)
    concatenated_vcfs.view{ it.getClass() } // data type of the output. It is class java.util.ArrayList
    concatenated_vcfs.map{ it -> tuple(*it) }.set{ concatenated_vcfs_tuple }
    concatenated_vcfs_tuple.view{ it.getClass() }

    // convert to a tab/tsv format
    convertToTab(concatenated_vcfs_tuple)
}
