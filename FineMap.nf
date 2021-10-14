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
        -profile slurm,singularity\
        -resume

    Mandatory arguments:
    --gwaslist          Tab separated file with header (column names: PhenoName, SumStat, SampleFile) and three columns. First column: phenotype name, second column: path to gwas summary statistics file, third column: path to the file which contains measurements for given phenotype (binary: 0, 1, NA; continuous: continuous numbers). Has to contain phenotype as a column name and column named "VKOOD" for sample IDs.
    --genotypefolder    Folder containing bgen files on which all those GWAS's were ran. File names have to contain the string "chr[1-23]".
    --imputationfile    Separate file containing imputation INFO score for each SNP in the genotype data.
    --outdir            Folder where output is written.
    
    Optional arguments:
    Filtering:
    --PvalThresh        GWAS P-value threshold for defining significant loci. Defaults to 5e-8.
    --Win               Genomic window to extract loci for finemapping. Defaults 1000000bp to either side of lead SNP.
    --MafThresh         MAF threshold to filter the input GWAS data. Defaults to 0.001.
    --InfoThresh        INFO score threshold to filter the input GWAS data. Defaults to 0.4.
    SuSiE settings:
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
    .map{row -> [ file(row.SumStat), file(row.SampleFile), row.PhenoName ]}
    .set { gwaslist_ch }

Channel
    .fromFilePairs(params.genotypefolder + '/*.{bgen,bgen.bgi,samples}', flat: true, size: -1)
    .set { genotype_ch }

Channel
    .fromPath(params.imputationfile)
    .set { imputationfile_ch }

Channel
    .fromPath('bin/Report_template.Rmd')
    .set { report_ch }

params.PvalThresh = 5e-8
params.Win = 1000000
params.MafThresh = 0.01
params.InfoThresh = 0.4
params.MaxIter = 100
params.MaxCausalSnps = 10


Channel
    .fromPath('bin/Report_template.Rmd')
    .set { report_ch }


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
      tuple file(gwas), file(samplefile), val(trait), file(imputationfile) from gwaslist_ch
      val Maf from params.MafThresh
      val Pval from params.PvalThresh
      val Imp from params.InfoThresh
      val Win from params.Win

    output:
      tuple val(gwas.simpleName), file(gwas), file(samplefile), val(trait), file(imputationfile) into ss_ch
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
        tuple val(gwas_id), val(region), file(ss_file), file(samplefile), val(trait), file(imputationfile) from gwaslist_ch
        val Maf from params.MafThresh
        val Pval from params.PvalThresh
        val Imp from params.InfoThresh

    output:
        tuple env(chr), val(gwas_id), val(region), file("*_region.txt"), file(samplefile), val(trait), file("variants_filter.txt") into gwaslist_ch2
        file("*.log") into region_log_ch
        """
        Rscript --vanilla ${baseDir}/bin/FilterSummaryStats.R  \
        -g ${ss_file} \
        -i ${imputationfile} \
        -r ${region} \
        -P ${Pval} \
        -M ${Maf} \
        -I ${Imp} \

        chr=\$(echo ${region} | sed -e "s/:.*//g")
        """
}

process GetGenChr {

    tag {GetGenChr}

    input:
        tuple fileId, file(bgen), file(bgi), file(sample) from genotype_ch

    output:
        tuple env(chr), file(bgen), file(bgi), file(sample) into bgen_ch

        """
        chr="\$(ls ${bgen} |\
        sed -e "s/.*chr//g" |\
        grep -oP "[0-9]{1,2}")"
        """
}

// Filter in only those regions where ther are at least 50 variants after all QC filters
gwaslist_ch2=gwaslist_ch2.filter{ it[6].countLines() > 50 }

gwaslist_ch3=gwaslist_ch2.combine(bgen_ch, by: 0)

process PrepareSampleList {
     
     tag {PrepareSampleList}

    input:
        tuple env(chr), val(gwas_id), val(region), file(ss_file), file(samplefile), val(trait), file(variants), file(bgen), file(bgi), file(sample) from gwaslist_ch3
    
    output:
        tuple env(chr), val(gwas_id), val(region), file(ss_file), file("*_PhenoList.txt"), val(trait), file(variants), file(bgen), file(bgi), file(sample) into gwaslist_ch4

        """
        Rscript --vanilla ${baseDir}/bin/PrepareSampleList.R \
        --samples_file ${samplefile} \
        --phenotype ${trait}
        """
}

process FilterBgen {
     
     tag {FilterBgen}

     input:
        tuple env(chr), val(gwas_id), val(region), file(ss_file), file(samplelist), val(trait), file(variants), file(bgen), file(bgi), file(sample) from gwaslist_ch4

     output:
        tuple val(gwas_id), val(region), file(ss_file), file(samplelist), val(trait), file(variants), file('*.bgen'), file('*.bgi'), file('*.sample') into gwaslist_ch5

        """
        # Parse region
        chr=\$(echo ${region} | sed -e "s/:.*//g")
        start=\$(echo ${region} | sed -e "s/.*://g" | sed -e "s/-.*//g")
        end=\$(echo ${region} | sed -e "s/.*://g" | sed -e "s/.*-//g")

        # Parse sample list
        awk '{print \$1}' ${samplelist} > samplelist.txt

        # NB!!! it seems that ref/alt are swapped in EstBB .bgen file
        plink2 \
        --bgen ${bgen} ref-last \
        --chr \${chr} \
        --from-bp \${start} \
        --to-bp \${end} \
        --keep samplelist.txt \
        --set-all-var-ids @:#\\\$a,\\\$r \
        --export bgen-1.3 \
        --new-id-max-allele-len 500 \
        --threads 4 \
        --memory 25600 \
        --out ${region}

        bgenix -index -g ${region}.bgen
        """
}

process FilterInfoBgen {
     
     tag {FilterInfoBgen}

    errorStrategy 'ignore' //TODO!!! fix this: bgen is done, yet it crashes for some cases (chr20)

     input:
        tuple val(gwas_id), val(region), file(ss_file), file(samplelist), val(trait), file(variants), file(bgen), file(bgi), file(sample) from gwaslist_ch5

     output:
        tuple val(gwas_id), val(region), file(ss_file), file(samplelist), val(trait), file(variants), file('*filtered.bgen'), file('*filtered.bgen.bgi'), file('*filtered.sample'), file('*.z') into gwaslist_ch6

        """

        plink2 \
        --bgen ${region}.bgen ref-last \
        --extract ${variants} \
        --threads 4 \
        --memory 25600 \
        --export bgen-1.3 \
        --out ${region}_filtered

        bgenix -index -g ${region}_filtered.bgen

        plink2 \
        --bgen ${region}_filtered.bgen \
        --make-just-bim \
        --out ${region}_filtered

        echo "rsid chromosome position allele1 allele2" > ${region}.z
        awk 'BEGIN { OFS = " "} { print \$2, \$1, \$4, \$5, \$6 }' ${region}_filtered.bim >> ${region}.z

        #mv ${region}_filtered.bgen ${region}.bgen
        #mv ${region}_filtered.bgen.bgi ${region}.bgen.bgi
        #mv ${region}_filtered.sample ${region}.sample

        """
}

process ExtractLdMatrix {

    tag {ExtractLdMatrix}

    input:
        tuple val(gwas_id), val(region), file(ss_file), file(samplelist), val(trait), file(variants), file(bgen), file(bgi), file(sample), file(zfile) from gwaslist_ch6
    
    output:
        tuple file(ss_file), val(region), file(samplelist), val(trait), file(variants), file('*.ld.gz') into prepare_inp_ch

        """
        # Prepare master file
        nrows=\$(wc -l ${samplelist} | awk '{print \$1}')
        nsamples=\$((\${nrows} - 1))

        echo "z;bgen;bgi;bcor;ld;n_samples;bdose" > master.txt

        echo "${region}.z;${region}_filtered.bgen;${region}_filtered.bgen.bgi;${region}.bcor;${region}.ld;\${nsamples};${region}.bdose" >> master.txt

        # Calculate correlation
        /usr/bin/ldstore \
        --in-files master.txt \
        --write-bcor \
        --read-only-bgen \
        --n-threads 4 

        /usr/bin/ldstore \
        --in-files master.txt \
        --rsids ${variants} \
        --bcor-to-text 

        gzip ${region}.ld
        """
}

process RunSuSiE {

    tag {RunSuSiE}

    publishDir path: "${params.outdir}/PipelineOutput", mode: 'copy', overwrite: true

    input:
        tuple file(ss_file), val(region), file(samplelist), val(trait), file(variants), file(ldmatrix) from prepare_inp_ch
        val MaxCausalSnps from params.MaxCausalSnps

    output:
        tuple file("*.susie.snp.gz"), file("*.susie.cred.gz"), file("*.susie.log"), file("*.rds") into output_ch,to_report_ch

        """
        # Calculate var_y (function adapted from FinnGen repo: https://github.com/FINNGEN/finemapping-pipeline/blob/37d75d0451a18d56e713fb5cd7a2907a5b328a7f/wdl/finemap_sub.wdl)
        var_y=\$(cat ${samplelist} | awk -v ph=Phenotype '
        BEGIN {
            FS = "\\t"
        }
        NR == 1 {
            for(i = 1; i <= NF; i++) {
                h[\$i] = i
            }
            exists=ph in h
            if (!exists) {
                print "Phenotype:"ph" not found in the given phenotype file." > "/dev/stderr"
                err = 1
                exit 1
            }
            cases = 0
            controls = 0
        }
        NR > 1 {
            vals[\$h[ph]] += 1
            if (\$h[ph] != "NA" && \$h[ph] != 0 && \$h[ph] != 1) {
                print "Phenotype:"ph" seems a quantitative trait. Setting var_y = 1." > "/dev/stderr"
                print 1.0
                err = 1
                exit 0
            }
        }
        END {
            if (!err) {
                phi = vals["1"] / (vals["1"]+vals["0"])
                var_y = phi*(1-phi)
                printf var_y
            }
        }')

        if [[ \$? -ne 0 ]]
        then
            echo "Error occurred while getting var_y from case control counts:" \$var_y
            exit 1
        fi

        # Prepare master file
        nrows=\$(wc -l ${samplelist} | awk '{print \$1}')
        nsamples=\$((\${nrows} - 1))

        # Run SuSiE
        Rscript --vanilla ${baseDir}/bin/RunSusieFinemappingSuffStat.R \
        --z ${ss_file} \
        --ld ${ldmatrix} \
        -n \${nsamples} \
        --L ${MaxCausalSnps} \
        --var-y \${var_y} \
        --snp ${trait}_${region}.susie.snp \
        --cred ${trait}_${region}.susie.cred \
        --log ${trait}_${region}.susie.log \
        --susie-obj ${trait}_${region}.susie.rds \
        --save-susie-obj \
        --write-alpha \
        --write-single-effect \
        --min-cs-corr 0.5
        """
}

report_ch = report_ch.collect()

process MakeReport {

    tag {MakeReport}

     publishDir path: "${params.outdir}/PipelineOutput", mode: 'copy', overwrite: true

    input:
    path report from report_ch
    tuple file("*.susie.snp.gz"), file("*.susie.cred.gz"), file("*.susie.log"), file("*.rds") from to_report_ch

    output:
        path "Report_FineMapping_*" into report_output

        """
        # Make report
        cp -L ${report} notebook.Rmd

        R -e 'library(rmarkdown);rmarkdown::render("notebook.Rmd", "html_document", 
        output_file = "Report_FineMapping.html"'
        """
}
