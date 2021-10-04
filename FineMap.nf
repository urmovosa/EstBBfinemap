def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     FineMappingPipeline v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run FineMap.nf \
        --gwaslist mygwases.txt\
        --genotypefolder genotypes\
        --imputationfile InfoFile.txt\
        --outdir output\
        -profile slurm\
        -resume

    Mandatory arguments:
    --gwaslist          Tab separated file with header and two columns. First column: path to gwas summary statistics file. Second column: path to sample ID list which was part of this GWAS.
    --genotypefolder    Folder containing bgzipped and tabixed .vcf files on which all those GWAS's were ran. File names have to contain string "chr[1-23]".
    --imputationfile    Separate file containing imputation INFO score for each SNP in the genotype data.
    --outdir            Folder where output is written.
    
    Optional arguments:
    --PvalThresh        GWAS P-value threshold for defining significant loci. Defaults to 5e-8.
    --Win               Genomic window to extract loci for finemapping. Defaults 1000000bp to either side of lead SNP.
    --MafThresh         MAF threshold to filter the input GWAS data. Defaults to 0.001.
    --InfoThresh        INFO score threshold to filter the input GWAS data. Defaults to 0.4.
    --MaxIter           Maximal number of iterations for SuSiE analysis. Defaults to 100.
    --MaxCausalSnps     Maximal number of causal SNPs tested in SuSiE analysis. Defaults to 10.
"""
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Define input channels
Channel.fromPath(params.gwaslist)
    .ifEmpty { error "Cannot find gwas list file in: ${params.gwaslist}" }
    .splitCsv(header: true, sep: ' ', strip: true)
    .map{row -> [ file(row.SumStat), file(row.SampleList) ]}
    .set { gwaslist_ch }

Channel
    .fromPath(params.genotypefolder)
    .map { gen -> [file("${gen}.vcf.gz"), file("${gen}.vcf.gz.tbi")] }
    .set { genotype_ch }

Channel
    .fromPath(params.imputationfile)
    .set { imputationfile_ch }

/* 
Channel
    .fromPath('bin/Report_template.Rmd')
    .set { report_ch } */

params.PvalThresh = 5e-8
params.Win = 1000000
params.MafThresh = 0.01
params.InfoThresh = 0.4
params.MaxIter = 100
params.MaxCausalSnps = 10

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
FineMappingPipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'FineMappingPipeline'
summary['Pipeline Version']         = workflow.manifest.version
summary['Run Name']                 = custom_runName ?: workflow.runName
summary['SS list']                  = params.gwaslist
summary['Genotype folder']          = params.genotypefolder
summary['Imputation QC file']       = params.imputationfile
summary['Max Memory']               = params.max_memory
summary['Max CPUs']                 = params.max_cpus
summary['Max Time']                 = params.max_time
summary['P threshold']              = params.PvalThresh
summary['Genomic window']           = params.Win
summary['MAF threshold']            = params.MafThresh
summary['INFO threshold']           = params.InfoThresh
summary['Max SuSiE iterations']     = params.MaxIter
summary['Max causal SNPs']          = params.MaxCausalSnps
summary['Output dir']               = params.outdir
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']            = params.awsregion
   summary['AWS Queue']             = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

gwaslist_ch = gwaslist_ch.combine(imputationfile_ch)


process FindRegions {

    tag {FindRegions}

    input:
      set file(gwas), file(samplefile), file(imputationfile) from gwaslist_ch
      val Maf from params.MafThresh
      val Pval from params.PvalThresh
      val Imp from params.InfoThresh
      val Win from params.Win

    output:
      set val(gwas.simpleName), file(gwas), file(samplefile), file(imputationfile) into ss_ch
      file("*regions.txt") into regions_ch

    """
    Rscript --vanilla ${baseDir}/bin/FindFinemapRegions.R  \
    -g ${gwas} \
    -i ${imputationfile} \
    -P ${Pval} \
    -W ${Win} \
    -M ${Maf} \
    -I ${Imp} \
    """
}

regions_ch
    .collectFile(name: 'regions.txt', keepHeader: true)
    .ifEmpty { error "No regions in the file" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.SS, row.Region ]}
    .set { gwaslist_ch }

gwaslist_ch = gwaslist_ch.combine(ss_ch, by: 0)

process PrepareSumstatRegions {

    tag {PrepareSumstatRegions}

    input:
        set val(gwas_id), val(region), file(ss_file), file(samplelist), file(imputationfile) from gwaslist_ch
        val Maf from params.MafThresh
        val Pval from params.PvalThresh
        val Imp from params.InfoThresh

    output:
        set val(gwas_id), val(region), file("*_region.txt"), file(samplelist), file("variants_filter.txt") into gwaslist_ch2

        """
        echo ${region}

        Rscript --vanilla ${baseDir}/bin/FilterSummaryStats.R  \
        -g ${ss_file} \
        -i ${imputationfile} \
        -r ${region} \
        -P ${Pval} \
        -M ${Maf} \
        -I ${Imp} \
        """
}

//gwaslist_ch2.view()


 process FilterVcf {
     
     tag {FilterVcf}

     input:
        set val(gwas_id), val(region), file(sumstats), file(samplelist), file(variants) from gwaslist_ch2
        path vcf_path from genotype_ch
     output:
        set val(gwas_id), val(region), file(sumstats), file(samplelist), file(variants), file("*_filtered.vcf.gz") into gwaslist_ch3

        """
        chr=\$(echo ${region} | sed -e "s/:.*//g")
        input_vcf=\$(ls *chr${chr}*.vcf.gz)

        bcftools view \
        --regions ${region} \
        --sample-file ${samplelist} \
        --force-samples \
        \${input_vcf} > ${vcf.simpleName}_filtered.vcf

        bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' ${vcf.simpleName}_filtered.vcf > ${vcf.simpleName}_filtered.vcf

        bcftools view \
        -e 'ID=@${variants}' \
        ${vcf.simpleName}_filtered.vcf > ${vcf.simpleName}_filtered.vcf

        bgzip ${vcf.simpleName}_filtered.vcf
        tabix -p ${vcf.simpleName}_filtered.vcf.gz
        """
 }


//  process ExtractLdMatrix {

//      tag {ExtractLdMatrix}

//      input:

//      output:

//         """
        
//         """
//  }

// process PrepareInputs {

//     tag {PrepareInputs}

//     input:

//     output:

//         """

//         """
// }

// process RunSuSiE {

//     tag {RunSuSiE}

//     input:

//     output:

//         """

//         """
// }
