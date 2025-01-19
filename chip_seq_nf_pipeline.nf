

//modules

FASTQC='fastqc/0.11.7'
FASTP='fastp/0.22.0'
BWA='bwa/0.7.17'
SAMTOOLS='samtools/1.16'
BEDTOOLS='bedtools/2.30.0'
DEEPTOOLS='deeptools/3.5.1'
MACS2='macs2/2.1.1'
R_LAN='r/4.3.2'


// making a parameter for the user to input a reference genome

params.mouse_ref_mm10='/gpfs/home/rj931/projects/chip_seq_practice/nextflow_workstation/mouse_genome_grcm38/mm10.fa.gz'


// process for downloading the mouse reference genome and gtf

process mouse_ref_genome_gtf {
    cache true
    publishDir './ref_genome_dir', mode:'copy', pattern:'*.fa.gz'
    publishDir './gtf_genome_dir', mode:'copy', pattern:'Mus_musculus*.gtf.gz'

    input:

    output:
    path '*.fa.gz', emit: ref_genome
    path 'Mus_musculus*.gtf.gz', emit: gtf_file

    """
    #!/bin/env bash

    # This is the mouse reference genome
    curl -O https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz

    # this is GRCm38.p6 version of the mouse genome
    #curl -o Mus_musculus.GRCm38.p6.fasta https://www.ebi.ac.uk/ena/browser/api/embl/GCA_000001635.8?download=true&gzip=true && gzip -f Mus_musculus.GRCm38.p6.fasta

    # this is the mouse gtf
    curl -O https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz


    """

}


process fastp_se { 
    cache true
    //executor 'slurm'
    memory '30 GB'
    
    publishDir './fastp_filter_output', mode:'copy', pattern:'*'

    input:
    tuple val(control), path(control_reads)
    tuple val(treatment), path(treatment_reads)
    val(filenames_control)
    val(filenames_treatment)
    
    output:

    tuple path("${filenames_control}_filt.fastq"), path("${filenames_treatment}_filt.fastq"), emit: filtered_fastq_c_t

    // output the .html to a channel

    path "*.html", emit: fastp_html
    
    script:

    """
    #!/bin/env bash

    module load $FASTP
    module load $FASTQC


    # this is single end use
    fastp \
    -i ${control_reads} \
    -o $filenames_control'_filt.fastq' \
    -h $filenames_control'.html' \
    --qualified_quality_phred=15 \
    --dedup \
    --dup_calc_accuracy=4 

    fastp \
    -i ${treatment_reads} \
    -o $filenames_treatment'_filt.fastq' \
    -h $filenames_treatment'.html' \
    --qualified_quality_phred=15 \
    --dedup \
    --dup_calc_accuracy=4 

    """
}



/*
// not doing this right now
process fastqc {

}
*/


// process for making the ref index

process bwa_index {
    cache true
    //memory '100 GB'

    publishDir './ref_indices', mode:'copy', pattern:'*'

    input:

    path ref

    

    output:
    path '*', emit: ref_index_files

    """
    #!/bin/env bash

    module load $BWA

    bwa index \
    -p $ref \
    -a bwtsw \
    $ref
    """
}


//process for alignment of filtered reads to reference genome.
// alignment of the filtered reads to the reference genome

process bwa_alignment {

    cache true

    publishDir './Chip_sam_storage', mode:'copy', pattern:'*.sam'

    input:

    tuple path(control_reads), path(treatment_reads)
    path ref
    path index_files
    path gtf
    tuple val(control_name), val(treatment_name)

    output:

    tuple path("${sam_control_file}"), path("${sam_treatment_file}"), emit: sam_files_tuple_ch


    script:

    sai_control_output = "${control_name}.sai"
    sam_control_file = "${control_name}.sam"

    sai_treatment_output = "${treatment_name}.sai"
    sam_treatment_file = "${treatment_name}.sam"

    """
    #!/bin/env bash

    module load $BWA

    # aligning the control files

    bwa aln \
    "${ref}" \
    "${control_reads}" \
    > "${sai_control_output}"


    bwa samse \
    "${ref}" \
    "${sai_control_output}" \
    "${control_reads}" \
    > "${sam_control_file}"


    # now aligning the treatment files

    bwa aln \
    "${ref}" \
    "${treatment_reads}" \
    > "${sai_treatment_output}"


    bwa samse \
    "${ref}" \
    "${sai_treatment_output}" \
    "${treatment_reads}" \
    > "${sam_treatment_file}"

    """
}


// I should use samtools to order the coordinates and then output that as a bam file

process samtools_convert {

    publishDir './Sorted_bam_files', mode: 'copy', pattern: '*.bam'
    publishDir './Sorted_bam_files', mode: 'copy', pattern: '*.bai'
    publishDir './Sorted_bai_files', mode: 'copy', pattern: '*.bai'
    publishDir './Bigwig_files', mode: 'copy', pattern: '*.bigwig'


    input:

    tuple path(sam_control), path(sam_treatment) 
    tuple val(control_basename), val(treatment_basename)
    
    output:

    tuple path("${out_control_file}"), path("${out_treatment_file}"), emit: sorted_bam_tuple_ch
    path '*.bai', emit: sorted_bai_ch
    tuple path("${out_control_bigwig}"), path("${out_treatment_bigwig}"), emit: bigwig_normalized_tuple
    
    tuple path("${out_control_unnormalized_bw}"), path("${out_treatment_unnormalized_bw}"), emit: bigwig_c_t_tuple

    script:

    out_control_file = "${control_basename}_sorted.bam"
    out_treatment_file = "${treatment_basename}_sorted.bam"

    out_control_bigwig = "${control_basename}_normalized.bigwig"
    out_treatment_bigwig = "${treatment_basename}_normalized.bigwig"

    out_control_unnormalized_bw = "${control_basename}_unnormalized.bigwig"
    out_treatment_unnormalized_bw = "${treatment_basename}_unnormalized.bigwig"

    """
    #!/bin/env bash

    module unload $MACS2
    module load $SAMTOOLS
    module load $DEEPTOOLS

    samtools sort \
    -O bam \
    -o "${out_control_file}" \
    "${sam_control}" 

    samtools sort \
    -O bam \
    -o "${out_treatment_file}" \
    "${sam_treatment}" 

    # I need to index the bam files for use in bamCoverage

    samtools index \
    "${out_control_file}"

    samtools index \
    "${out_treatment_file}" 


    # If I need fast random access, I can use samtools index and use the coordinate sorted bam file as input
    # the output will be a .bai file.

    # maybe make a normalized bigwig file using rpm

    # using deeptools I can convert bed to bigwig and normalize by rpm or other methods

    # making the control bigwig file first
    # using the table that shows the effective genome size for a mouse genome grcm39
    bamCoverage \
    -b "${out_control_file}" \
    -o "${out_control_bigwig}" \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2654621783

    # now doing the same for treatment bam to bigwig 
    bamCoverage \
    -b "${out_treatment_file}" \
    -o "${out_treatment_bigwig}" \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2654621783

    

    # I want to get the bigwig unnormalized files 

    bamCoverage \
    -b "${out_control_file}" \
    -o "${out_control_unnormalized_bw}" 
    

    # now doing the same for treatment bam to bigwig 
    bamCoverage \
    -b "${out_treatment_file}" \
    -o "${out_treatment_unnormalized_bw}" 


    ######### NOT USING THIS SECTION #############################
    #scaling_factor_control=\$(samtools idxstats "{out_control_file}" | '{sum+=\$3}END{print sum/1000000}')
    
    #bedtools genomecov \
    #-ibam "{out_control_file}" \
    #-scale scaling_factor_control \
    #-bg \
    #> "{out_control_bedgraph}"
    

    # now converting the bedgraph to bigwig
    
    
    #scaling_factor_treatment=\$(samtools idxstats "{out_treatment_file}" | '{sum+=\$3}END{print sum/1000000}')


    #bedtools genomecov \
    #-ibam "{out_treatment_file}" \
    #-scale scaling_factor_treatment \
    #-bg \
    #> "{out_treatment_bedgraph}"


    """
}


// now to use macs2 to call the peaks

process peak_calling_mouse {

    publishDir './macs_outputs', mode: 'copy', pattern:'*'

    input:
    // i am doing it this way because i didnt make separate processes for each treatment and control groups
    tuple path(control_bam), path(treatment_bam)
    tuple val(control_name), val(treatment_name)

    // inputting the bigwig files to make the peak count table with deeptools
    //tuple path(control_bigwig), path(treatment_bigwig)
    //tuple val(control_bw_basename), val(treatment_bw_basename)
    
    output:

    path "*", emit: peak_files

    script:
    // macs2 does not allow for out file name
    experiment_name = "${control_name}_vs_${treatment_name}"

    //peak_count_table_out = "${control_bw_basename}_vs_${treatment_bw_basename}"

    """
    #!/bin/env bash
    # module unload \$DEEPTOOLS
    module unload $MACS2
    module unload python
    module load $MACS2
    

     # the -g flag is showing that the genome size is comming from the mouse
    # the -t flag takes the treatment group of experiments
    # the -c flag takes the control group of experiments
    # the -f flag specifies that I am using bam files if BAM is set
    # -B flag is whether or not to save extended fragment pileup
    # --trackline tells macs to include trackline with gedgraph files, not using this due to downstream issues in chipseeker
    # --SPMR saves signal per million reads for fragment pileup profiles
    # -q flag is the minimum FDR cutoff for peak detection
    # --cutoff-analysis will give a summary of number of total peaks that can be called by different thresholds. so we can make a decision
    # -n flag is the name of the experiment that will be used as prefix for output file names
    macs2 callpeak \
    -t "${treatment_bam}" \
    -c "${control_bam}" \
    -f BAM \
    -g 'mm' \
    -B \
    --SPMR \
    -q 0.05 \
    --cutoff-analysis \
    --outdir '.' \
    -n "${experiment_name}"


    ########### DOING THIS IN THE MAKE PEAK COUNT TABLE PROCESS ################
    # I need to unload macs2 and then load deeptools for it to work
    #module unload \$MACS2
    #module unload \$DEEPTOOLS
    #module load \$DEEPTOOLS

    # now I want to create a peak count table using deeptools command multiBigwigSummary BED-file

    #multiBigwigSummary BED-file \
    #-b "\${control_bigwig}" "\${treatment_bigwig}" \
    #-o "\${peak_count_table_out}" \
    #--BED "*.bed"
    
    """

}

process make_peak_count_table {

    publishDir './macs_outputs', mode: 'copy', pattern:'*'

    input:

    // inputting the bigwig files to make the peak count table with deeptools
    tuple path(control_bigwig), path(treatment_bigwig)
    tuple path(control_unnormalized_bw), path(treatment_unnormalized_bw)
    tuple val(control_bw_basename), val(treatment_bw_basename)
    tuple val(control_unnormalized_base), val(treatment_unnormalized_base)
    tuple file(lambda_bdg), file(analysis_txt), file(model_r), file(peak_narrowpeak), file(peaks_xls), file(summits_bed), file(pileup_bdg)

    output:

    path "${peak_count_table_out}*", emit: peak_count_table_ch
    
    path "${out_unnormalized_basename}*", emit: unnormalized_peak_count_ch

    script:

    peak_count_table_out = "${control_bw_basename}_vs_${treatment_bw_basename}"

    out_unnormalized_basename = "${control_unnormalized_base}_vs_${treatment_unnormalized_base}"

    bed_file = "${summits_bed}"

    """
    # I need to unload macs2 and then load deeptools for it to work
    module purge
    #module unload \$MACS2
    module load $DEEPTOOLS

    # now I want to create a peak count table using deeptools command multiBigwigSummary BED-file

    multiBigwigSummary BED-file \
    -b "${control_bigwig}" "${treatment_bigwig}" \
    -o "${peak_count_table_out}_peak_counts.npz" \
    --outRawCounts "${peak_count_table_out}_peak_counts.tsv" \
    --BED "${bed_file}"

    multiBigwigSummary BED-file \
    -b "${control_unnormalized_bw}" "${treatment_unnormalized_bw}" \
    -o "${out_unnormalized_basename}_peak_counts.npz" \
    --outRawCounts "${out_unnormalized_basename}_peak_counts.tsv" \
    --BED "${bed_file}"

    """

}



// I need to create the plots with ngsplot

process ngs_plot {

    publishDir './Visualizations', mode: 'copy', pattern: '*.pdf'

    input:

    path bam
    path bai

    output:
    
    path "*.pdf"

    script:


    """
    #!/bin/env bash

    module load $R_LAN
    # setting up the environment and tool to use ngsplot with the mouse genome

    # first downloading the mouse genome from their google drive
    # i got the file_id and put it equal to the link for google drive provided elsewhere

    #curl -o ngsplotdb_mm10_75_3.00.tar.gz \
    # "https://drive.google.com/uc?export=download&id=0B5hDZ2BucCI6NXNzNjZveXdadU0"
    


    #git clone  https://github.com/shenlab-sinai/ngsplot.git

    #ngsplotdb.py install -y ngsplotdb_mm10_75_3.00.tar.gz


    # This is for the treatment file, which I put second in the channel in this nf script
    # creating an enrichment plot of the transcription start site (tss) using ngsplot
     ngs.plot.r \
     -G mm10 \
     -R tss \
     -C "${bam[1]}" \
     -O klf4.tss \
     -T klf4 \
     -L 3000 \
     -FL 50 \
     -RR 50 \
     -SC global \
     -P 12 \
     -D ensembl \
     -MQ 20 \
     -SE 0 \
     -FC 0.07 \
     -LOW 0 \
     -GO total


    # This is for the control file which I put first in the nf script
    # creating the same graphs for the input gene

    ngs.plot.r \
     -G mm10 \
     -R tss \
     -C "${bam[0]}" \
     -O input.tss \
     -T Input \
     -L 3000 \
     -FL 50 \
     -CO white:blue \
     -CD 1 \
     -RR 50 \
     -SC global \
     -P 12 \
     -D ensembl \
     -MQ 20 \
     -SE 0 \
     -FC 0.07 \
     -LOW 0 \
     -GO total

    """

}

// making a process for chipseeker

process chipSeeker_mouse {
    
    cache true
    publishDir './annotated_peaks', mode: 'copy', pattern: '*.tsv'
    
    input:

    tuple file(lambda_bdg), file(analysis_txt), file(model_r), file(peak_narrowpeak), file(peaks_xls), file(summits_bed), file(pileup_bdg)

    output:

    path "*.tsv", emit: annotated_peaks_ch

    script:


    """
    #!/usr/bin/env Rscript --no-save

    # can comment this out if you have chipseeker already installed

    BiocManager::install("ChIPseeker")
    
    library(ChIPseeker)
    
    # place holder code
    
    peak_file = readPeakFile("${peak_narrowpeak}")
    
    # I think there are gtf/ knowngene packages in r under Txdb moniker
    
    BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene") # I will install the version for mm10
    
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    
    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
    
    ann_peaks_df = annotatePeak(peak_file, TxDb = txdb, tssRegion= c(-3000, 3000), level = "gene", assignGenomicAnnotation= TRUE, genomicAnnotationPriority= "Downstream", annoDb="org.Mm.eg.db")
    
    # added the annoDb parameter to include genename, ensembl, and symbol columns
    peak_ann_df = as.data.frame(ann_peaks_df)
    
    # saving the table to a file tsv
    write.table(peak_ann_df, 'ann_peaks_placeholder_name.tsv', sep ='\t')



    """
}

// process to use deseq2 to find differendial binding of the tf in the Chip-seq data


process deseq2_db {

    publishDir './R_outputs' , mode: 'copy', pattern: '*tsv'

    input:

    tuple path(peak_table_npz), path(peak_table_tsv)

    output:

    path 'peak_count_table_w_DB_peaks.tsv', emit: peak_table_db

    script:

    peak_counts = "${peak_table_tsv}"

    """
    #!/usr/bin/env Rscript --no-save

    # loading in deseq2 for differential binding analysis between the input and treatment

    library(DESeq2)

    # this will be using normalized counts taken from the bigwig files. I know that deseq2 does
    # not take normalized counts but I will have another process where I get the non-normalized counts from a bam file instead

    # first read in the tsv file

    peak_count_table = read.csv("${peak_counts}", sep='\t', header = TRUE)

    peak_count_table\$key = paste(peak_count_table[,1], peak_count_table[,2], peak_count_table[,3], sep = ':')

    input_df = data.frame(row.names = peak_count_table\$key, Input = peak_count_table[,4], Treatment = peak_count_table[,5])

    # now I need to create a colData df to store my experimental design

    colData = data.frame(row.names = colnames(input_df), conditions = colnames(input_df))

    colData\$condition <- factor(colData\$conditions)

    # just testing how I should handel and error in this case

    #if (all(row.names(colData) == colnames(input_df)) == False) {
    #    print('Your colData rows names do not match your count data column names')  }
    #else {
    #    print('colData row names and count data column names match') }

    # now inporting the count matrix and coldata into deseq2

    # since I do not have multiple biological replicates for each condition I need to make the design = ~ 1
    # this is recommended by michael love on bioconductor blog questions

    dds = DESeqDataSetFromMatrix(countData=input_df, colData= colData, design = ~ 1)

    dds = DESeq(dds)
    res = results(dds)

    # using normal instead of apeglm. resLFC is useful for visualization by shrinking the effect size
    resLFC = lfcShrink(dds, coef="Intercept", type="normal")

    # ordering the rows by the padj column

    resorder = res[order(res\$padj),]

    # keeping the rows that have a padj value less than or eaqual to 0.05
    
    keep = which(resorder\$padj <= 0.05)

    resorder_signif = resorder[keep,] # this results in only 113 of 2528 significant differential binding peaks

    # to check the amount of upregulated and downregulated DB peaks I can use summary
    # but I need to use the function to filter with alpha so the accurate threshold is recorded

    res05 = results(dds, alpha=0.05)
    summary(res05) # as can see this method gives back more peaks differentially expressed that the other method 165 instead of 113

    # I can remove and count how many differential binding peaks are in res05
    sum(res05\$padj < 0.05, na.rm= TRUE) # this gives back how many peaks are differentially bound by a tf = 165
    
    summary(resorder_signif) # gives back 113 DB peaks. this has less than the previous method.

    
    # finding which rows have the differential binding

    keep_diff = which(res05\$padj < 0.05, na.exclude)

    res05_diff = res05[keep_diff,]

    # now to add a new column to the original table to show which rows are differentially bound


    names_diff_keep = which( (peak_count_table\$key) %in% row.names(res05_diff)) # hrecording the row names that are differentially expressed with the ones in the full peak table

    new_column = rep(FALSE, nrow(peak_count_table)) # making a vector of false

    new_column[names_diff_keep] = TRUE # making sure the new column knows which rows are the differentially bound peaks


    peak_count_table\$Differentially_bound = new_column # now adding the column that shows which peaks are diffentially bound

    write.table(peak_count_table, './peak_count_table_w_DB_peaks.tsv', sep='\t')

    """

}



workflow {

    mm10_ref=Channel.fromPath(params.mouse_ref_mm10, checkIfExists: true)

    mouse_ref_genome_gtf()
    //getting the ref and gtf channels
    mouse_ref_genome = mouse_ref_genome_gtf.out.ref_genome

    mouse_gtf_genome = mouse_ref_genome_gtf.out.gtf_file

    mouse_ref_genome.view()
    mouse_gtf_genome.view()


    // this can be over ridden by setting the --se_reads parameter for single end reads
    // this is NOT control but igg. this experiment does not have the control just the treatment
    params.se_reads_control = '../chip_fastqs/control/chip_*.fastq'
    params.se_reads_treatment = '../chip_fastqs/treatment/chip_*.fastq'
    Channel.fromPath(params.se_reads_control, checkIfExists: true)
            .map{file -> tuple('control', file)}
            .set{SE_reads_control }
    SE_reads_control.view()
    
    Channel.fromPath(params.se_reads_treatment, checkIfExists: true)
            .map{file -> tuple('control', file)}
            .set{SE_reads_treatment}
    
    SE_reads_control
        .map { file -> file[1].baseName}
        .set {filenames_control}

    filenames_control.view()

    SE_reads_treatment
        .map { file -> file[1].baseName}
        .set {filenames_treatment}
    
    
    // I want to try getting the srr files using nextflow
    //params.sra_txt = '../SRR_Acc_List.txt'
    //SRR_txt_list = Channel.fromPath(params.sra_txt, checkIfExists: true).splitText().collect() //.view()
    //SRR_txt_list.view()
    //ids = SRR_txt_list
    //println ids
    // now that I have my SRR list, I can put it into the channel that works with sra

    //sra_acc = Channel.fromSRA(SRR_txt_list, apiKey:'a5b95701789e19a261c98a7fcdb920f46f08', cache:true)

    //sra_acc.view()
    //Channel.fromSRA()

    //control_fastq =  channel.fromPath('../chip_fastqs/goat_igg.fastq')

    fastp_se(SE_reads_control, SE_reads_treatment, filenames_control, filenames_treatment )

    filt_files_tuple = fastp_se.out.filtered_fastq_c_t
    filt_files_tuple.view()
    
    
    
    // get the base names of the filtered files
    filt_files_tuple
        .map {file -> file.baseName}
        .set {filt_basenames}
    filt_basenames.view()

    
    html_fastp_file = fastp_se.out.fastp_html.collect() 

    //filt_files_tuple.view()
    html_fastp_file.view()

    //bwa_index(mouse_ref_genome)
    
    // other version of mouse ref from above
    bwa_index(mm10_ref)

    // now taking those reference index files and giving them their own channel

    index_files = bwa_index.out.ref_index_files

    index_files.view()
        //.tap {amb_index, ann_index, bwt_index, pac_index, sa_index}

   


    //bwa_alignment(filt_files_tuple, mouse_ref_genome, index_files, mouse_gtf_genome, filt_basenames)

    // using the mm10 genome
    bwa_alignment(filt_files_tuple, mm10_ref, index_files, mouse_gtf_genome, filt_basenames)

    sam_output_files_ch = bwa_alignment.out.sam_files_tuple_ch

    sam_output_files_ch.view()

    

    // give the basename of the sam files.

    sam_output_files_ch
        .map{file -> file.baseName}
        .set{sam_basenames}

    samtools_convert(sam_output_files_ch, sam_basenames)

    sorted_bams = samtools_convert.out.sorted_bam_tuple_ch
    sorted_bams.view()

    big_wig_tuple_ch = samtools_convert.out.bigwig_normalized_tuple
    big_wig_unnormalized_ch = samtools_convert.out.bigwig_c_t_tuple

    bai_files = samtools_convert.out.sorted_bai_ch
    // making the basename again from the sorted bams channel

    sorted_bams
        .map{file -> file.baseName}
        .set{bams_basename}

    big_wig_tuple_ch
        .map {file -> file.baseName}
        .set {big_wig_basename}

    peak_calling_mouse(sorted_bams, bams_basename )

    called_peaks_ch = peak_calling_mouse.out.peak_files
    called_peaks_ch.view()

    big_wig_unnormalized_ch
        .map{file -> file.baseName}
        .set{unnormalized_basename_bw}


    make_peak_count_table(big_wig_tuple_ch, big_wig_unnormalized_ch, big_wig_basename, unnormalized_basename_bw, called_peaks_ch )

    peak_table_ch = make_peak_count_table.out.unnormalized_peak_count_ch
    peak_table_ch.view()


    ngs_plot(sorted_bams, bai_files)

    chipSeeker_mouse(called_peaks_ch)

    deseq2_db(peak_table_ch)

    
}