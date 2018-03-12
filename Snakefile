#shell.rm prefix("set -eo pipefail; ")

configfile: "config.yaml"

###################### TARGETS ################################

ALL_SAMPLES = config["samples"] + config["inputs"]
ALL_BAM     = expand("04aln/{sample}.unsorted.bam 04aln/{sample}.sorted.bam 04aln/{sample}.sorted.bam.bai 04aln/{sample}.sorted.bam.flagstat".split(), \
        sample = ALL_SAMPLES)

ALL_CLEAN_FASTQ = expand("03seqClean/{sample}_{num}_{pe}.fastq", sample = ALL_SAMPLES, num=[1,2], pe=["paired","unpaired"])
ALL_FASTQC  = expand("02fqc/{sample}_{num}_fastqc.zip", sample = ALL_SAMPLES, num=[1,2])

ALL_UNSORTED_BAM = expand("04aln/{sample}.unsorted.bam", sample = ALL_SAMPLES)
ALL_SORTED_BAM = expand("04aln/{sample}.sorted.bam", sample = ALL_SAMPLES)
ALL_INDEXED_BAM = expand("04aln/{sample}.sorted.bam.bai", sample = ALL_SAMPLES)

NARROW_PEAKS = expand("05peak/{sample}_peaks.narrowPeak", sample = config["samples"]) 
FOLD_ENRICH = expand("05peak/{sample}_FE.bdg 05peak/{sample}_FE.bw".split(), sample = config["samples"])
BG_SUBTRACT = expand("05peak/{sample}_BGsub.bdg 05peak/{sample}_BGsub.bdg".split(), sample = config["samples"])

BROAD_PEAKS = expand("06broad_peak/{sample}_peaks.broadPeak", sample = config["samples"])
FOLD_ENRICH_BP_BW = expand("06broad_peak/{sample}_FE.bdg 06broad_peak/{sample}_FE.bw".split(), sample = config["samples"])


# Prints out
NB_SAMPLES = len(ALL_SAMPLES)
#~ print("Number of samples: " + str(NB_SAMPLES), file=sys.stderr)

#~ for smp in ALL_SAMPLES:
  #~ print("Sample " + smp + " will be processed", file=sys.stderr)


#################### RULES ################################

rule all:
    input: FOLD_ENRICH + BG_SUBTRACT
    #input: ALL_FASTQC + ALL_CLEAN_FASTQ + ALL_BAM + NARROW_PEAKS + BROAD_PEAKS # + FOLD_ENRICH + FOLD_ENRICH_BP_BW

rule fastqc:
    input:  "data/{sample}.fastq.gz"
    output: "02fqc/{sample}_fastqc.zip", "02fqc/{sample}_fastqc.html"
    log:    "00log/{sample}_fastqc"
    threads: 2
    resources:
        mem  = 2 ,
        time = 20
    message: "fastqc {input}: {threads} / {resources.mem}"
    shell:
        """
        fastqc {input[0]} -o 02fqc 
        """

#~ ## use trimmomatic to trim low quality bases and adaptors
rule clean_fastq:
    input:
        r1 = "data/{sample}_1.fastq.gz", 
        r2 = "data/{sample}_2.fastq.gz"
    output:  
        out1p = "03seqClean/{sample}_1_paired.fastq",
        out1u = "03seqClean/{sample}_1_unpaired.fastq",
        out2p = "03seqClean/{sample}_2_paired.fastq",
        out2u = "03seqClean/{sample}_2_unpaired.fastq",
    log:     "00log/{sample}_clean_fastq"
    params: "-phred33"
    threads: 4
    resources:
        mem= 2,
        time= 30
    message: "clean_fastq {input}: {threads} threads / {resources.mem}"
    shell:
        """
        trimmomatic PE \
        {params} \
        {input.r1} {input.r2} \
        {output.out1p} {output.out1u} \
        {output.out2p} {output.out2u} \
        ILLUMINACLIP:/home/rybalko/miniconda3/share/trimmomatic-0.36-3/adapters/TruSeq3-PE-2.fa:2:30:10 \
        TRAILING:25 CROP:65 MINLEN:60 2> {log}
        """
# nohup java -jar /home/rybalko/tools/Trimmomatic-0.36/trimmomatic-0.36.jar

## ---- Read mapping with Bowtie2--------------------------------------
# The NGS reads are aligned with Bowtie2 against the reference genome 
# sequence (Langmead and Salzberg, 2012). 
# In ChIP-Seq experiments, it is usually more appropriate to eliminate 
# reads mapping to multiple locations. To achieve this, users 
# can remove the argument setting ‘-k 50 –non-deterministic’ in the params.

rule align_bowtie2:
    input:
        genome = config["idx_bt2"]+"mm9/mouse.mm9.reference.genome.allChr.longChrName.fa",
        r1 = "03seqClean/{sample}_1_paired.fastq", 
        r2 = "03seqClean/{sample}_2_paired.fastq"
    output: temp("04aln/{sample}.sam")
    threads: 10
    params: bowtie = " -p 10 --very-sensitive"
    resources:
        mem    = 10,
        time   = 45
    message: "aligning {input}: {threads} threads / {resources.mem}"
    log: "00log/{sample}.align"
    shell:
        """
        bowtie2 {params.bowtie} \
        -x {input.genome} -1 {input.r1} -2 {input.r2} -S {output} 2> {log}
        """

rule sam2bam:
    input: "04aln/{sample}.sam"
    output: "04aln/{sample}.unsorted.bam"
    threads: 10
    resources:
        mem    = 10,
        time   = 45
    message: "sam --> bam {input}: {threads} threads / {resources.mem}"
    shell:
        """
        samtools view -Sb {input} >  {output}
        """

rule sort_bam:
    input: "04aln/{sample}.unsorted.bam"
    output: "04aln/{sample}.sorted.bam"
    log:    "00log/{sample}.sort_bam"
    threads: 10
    params: " -m 1G"
    resources:
        mem    = 10,
        time   = 45
    message: "bam --> sotred bam {input}: {threads} threads / {resources.mem}"
    shell:
        """
        samtools sort {params} -@ {threads} {input} -o {output} 2> {log}
        """

rule index_bam:
    input:  "04aln/{sample}.sorted.bam"
    output: "04aln/{sample}.sorted.bam.bai"
    log:    "00log/{sample}.index_bam"
    threads: 10
    resources:
        mem   = 500,
        time  = 10
    message: "index_bam {input}: {threads} threads / {resources.mem}"
    shell:
        """
        samtools index {input} {output} 2> {log}
        """

rule flagstat_bam:
    input:  "04aln/{sample}.sorted.bam"
    output: "04aln/{sample}.sorted.bam.flagstat"
    log:    "00log/{sample}.flagstat_bam"
    threads: 1
    params:
        mem   = "500M",
        time  = "10"
    message: "flagstat_bam {input}: {threads} threads / {params.mem}"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule merge_inputs:
    input:   expand("04aln/{sample}.sorted.bam", sample = config["inputs"])
    output:  "04aln/control.sorted.bam"
    log:     "00log/merge_controls"
    threads: 2
    params:
        mem = "1G",
        time = "15"
    message: "merge_controls {input}: {threads} threads / {params.mem}"
    shell:
        """
        samtools merge -r -@{threads} {output} {input}
        """

#Important to produce model in R and have  {sample}_model.r and
#{sample}_model.pdf as output file. This doesn't work with --nomodel
####################################################
rule find_narrow_peaks:
    input:  
        smpl = "04aln/{sample}.sorted.bam", 
        inp = "04aln/control.sorted.bam"
    output: "05peak/{sample}_model.r", "05peak/{sample}_peaks.narrowPeak",
            "05peak/{sample}_peaks.xls", "05peak/{sample}_summits.bed",
            #"05peak/{sample}_model.pdf"
    log:    "00log/{sample}.find_narrow_peaks"
    threads: 2
    params:
        macs = "--gsize mm --bdg",  
        mem  = "2G",
        time = "15"
    message: "find_narrow_peaks {input}: {threads} threads / {params.mem}"
    shell:
        """
        macs2 callpeak -t {input.smpl} \
                -c {input.inp} -f BAM {params.macs} \
                --outdir 05peak  -n {wildcards.sample} -q 0.001 &> {log}
        cd 05peak && Rscript {wildcards.sample}_model.r
        """

#~ Run MACS2 bdgcmp to generate fold-enrichment or logLR track
 #~ macs2 bdgcmp -t ${SAMPLE}_treat_pileup.bdg -c ${SAMPLE}_control_lambda.bdg -o ${SAMPLE}_logLR.bdg -m logLR -p 0.00001
 #~ * -m FE means to calculate fold enrichment. Other options can be logLR for log likelihood, subtract for subtracting noise from treatment sample.
 #~ * -p sets pseudocount. This number will be added to 'pileup per million reads' value. You don't need it while generating fold enrichment track because control lambda will always >0. But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.

rule fold_enrich:
    input:
        smpl = "05peak/{sample}_treat_pileup.bdg",
        inp =  "05peak/{sample}_control_lambda.bdg"
    output: "05peak/{sample}_FE.bdg"
    log: "00log/{sample}.fold_enrich"
    threads: 2
    message: "generate fold-enrichment {input}: {threads}"
    shell:
        """
        macs2 bdgcmp -t {input.smpl} -c {input.inp} -o {output} -m FE &> {log}
        """


rule bg_subtract:
    input:
        smpl = "05peak/{sample}_treat_pileup.bdg",
        inp =  "05peak/{sample}_control_lambda.bdg"
    output: "05peak/{sample}_BGsub.bdg"
    log: "00log/{sample}_merged.bg_subtract"
    threads: 2
    message: "generate subtract {input}: {threads}"
    shell:
        """
        macs2 bdgcmp -t {input.smpl} -c {input.inp} -o {output} -m subtract &> {log}
        """

rule bedGraph2bigWig:
    input:
        ref = config["idx_bt2"]+"mm9/mm9.len",
        smpl_fe = "05peak/{sample}_FE.bdg",
        smpl_bg = "05peak/{sample}_BGsub.bdg"
    output: "05peak/{sample}_FE.bw","05peak/{sample}_BGsub.bw" 
    log:  "00log/{sample}.bdg2bw"
    threads: 2
    message: "generate bigWig from bedGraph {input}: {threads}"
    shell:
        """
        ./bdg2bw {input.smpl_fe} {input.ref} > {output[0]} 
        ./bdg2bw {input.smpl_bg} {input.ref} > {output[1]} 
        """


rule find_broad_peaks:
    input:  
        smpl = "04aln/{sample}.sorted.bam", 
        inp = "04aln/control.sorted.bam"
    output: "06broad_peak/{sample}_peaks.xls", "06broad_peak/{sample}_peaks.broadPeak"
    log:    "00log/{sample}.find_broad_peaks"
    threads: 2
    params:
        macs = "--gsize mm --bdg",
        mem  = "2G",
        time = "15"
    message: "find_broad_peaks {input}: {threads} threads / {params.mem}"
    shell:
        """
        macs2 callpeak -t {input.smpl} \
                -c {input.inp} -f BAM {params.macs} \
                --broad --broad-cutoff 0.1 --nomodel --extsize 150 \
                --outdir 06broad_peak -n {wildcards.sample} -q 0.001 &> {log}
        """


rule fold_enrichBP:
    input:
        smpl = "06broad_peak/{sample}_treat_pileup.bdg",
        inp =  "06broad_peak/{sample}_control_lambda.bdg"
    output: "06broad_peak/{sample}_FE.bdg"
    log: "00log/{sample}.fold_enrichBP"
    threads: 2
    message: "generate fold-enrichment {input}: {threads}"
    shell:
        """
        macs2 bdgcmp -t {input.smpl} -c {input.inp} -o {output} -m FE &> {log}
        """

rule bedGraph2bigWigBP:
    input:
        ref = config["idx_bt2"]+"mm9/mm9.len",
        smpl = "06broad_peak/{sample}_FE.bdg"
    output: "06broad_peak/{sample}_FE.bw"
    log:  "00log/{sample}.bdg2bwBP"
    threads: 2
    message: "generate bigWig from bedGraph {input}: {threads}"
    shell:
        """
        bdg2bw {input.smpl} {input.ref} > {output}  
        """


### Visualization
#~ To visualize the workflow, one can use the option --dag. This creates a representation of the DAG in the graphviz dot language which has to be postprocessed by the graphviz tool dot. E.g. to visualize the DAG that would be executed, you can issue:
#~ $ snakemake --dag | dot | display
#~ For saving this to a file, you can specify the desired format:
#~ $ snakemake --dag | dot -Tsvg > dag.svg
#~ To visualize the whole DAG regardless of the eventual presence of files, the forceall option can be used:
#~ $ snakemake --forceall --dag | dot -Tpdf > dag.pdf
#~ Of course the visual appearance can be modified by providing further command line arguments to dot.
# Rule graph
# snakemake --rulegraph | dot -Tpng > rules.png





###############################################################



######################################## EXAMPLES OF RULES ###########################
######################################################################################

#~ localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

## list all samples
#~ CASES = config["samples"]
#~ CONTROLS = config["inputs"]


## list BAM files
#~ CONTROL_BAM = expand("04aln/{sample}.sorted.bam", sample=CONTROLS)
#~ CASE_BAM = expand("04aln/{sample}.sorted.bam", sample=CASES)

## create target for peak-calling: will call the rule call_peaks in order to generate bed files
## note: the "zip" function allow the correct pairing of the BAM files
#ALL_PEAKS   = expand("05peak/{case}_vs_{control}_peaks.bed", zip, case=CASES, control=CONTROLS)
#ALL_FASTQ   = expand("data/{sample}_{num}.fastq", sample = ALL_SAMPLES, num=[1,2])




## for each sample, there are multiple fastq.gz from different lanes, merge them
## because the number of the fastq.gz files in the folder is not fixed, need to be
## determined by a function

#~ import glob
#~ def get_fastqs(wildcards):
    #~ return glob.glob("rawfastqs/"+ wildcards.sample+"/"+ wildcards.sample + "_L00[0-9].fastq.gz")

#~ rule merge_fastqs:
    #~ input: get_fastqs
    #~ output: "01seq/{sample}.fastq"
    #~ log: "00log/{sample}_unzip"
    #~ threads: 1
    #~ message: "merging fastqs gunzip -c {input} > {output}"
    #~ shell: "gunzip -c {input} > {output} 2> {log}"




#~ rule align:
    #~ input:  "03seqClean/{sample}_clean.fastq"
    #~ output: temp("04aln/{sample}.unsorted.bam")
    #~ threads: 10
    #~ params: bowtie = "--chunkmbs 320 --best -p 10 "
    #~ resources:
        #~ mem    = 10,
        #~ time   = 45
    #~ message: "aligning {input}: {threads} threads / {resources.mem}"
    #~ log: "00log/{sample}.align"
    #~ shell:
        #~ """
        #~ module load bowtie/1.1.1 samtools/1.2
        #~ bowtie {params.bowtie} --threads={threads} {config[idx_bt1]} {input} 2> {log} \
            #~ | samtools view -Sb - \
            #~ >  {output}
        #~ """

#~ rule sort_bam:
    #~ input:  "04aln/{sample}.unsorted.bam"
    #~ output: "04aln/{sample}.sorted.bam"
    #~ log:    "00log/{sample}.sort_bam"
    #~ threads: 10
    #~ resources:
        #~ mem  = 12,
        #~ time = 15
    #~ message: "sort_bam {input}: {threads} threads / {resources.mem}"
    #~ shell:
        #~ """
        #~ module load samtools/1.2
        #~ samtools sort -m 1G -@ {threads} -O bam -T {output}.tmp {input} > {output} 2> {log}
        #~ """

#~ rule index_bam:
    #~ input:  "04aln/{sample}.sorted.bam"
    #~ output: "04aln/{sample}.sorted.bam.bai"
    #~ log:    "00log/{sample}.index_bam"
    #~ threads: 1
    #~ resources:
        #~ mem   = 500,
        #~ time  = 10
    #~ message: "index_bam {input}: {threads} threads / {resources.mem}"
    #~ shell:
        #~ """
        #~ module load samtools/1.2
        #~ samtools index {input} 2> {log}
        #~ """

#~ rule flagstat_bam:
    #~ input:  "04aln/{sample}.sorted.bam"
    #~ output: "04aln/{sample}.sorted.bam.flagstat"
    #~ log:    "00log/{sample}.flagstat_bam"
    #~ threads: 1
    #~ resources:
        #~ mem   = 500,
        #~ time  = 10
    #~ message: "flagstat_bam {input}: {threads} threads / {resources.mem}"
    #~ shell:
        #~ """
        #~ module load samtools/1.2
        #~ samtools flagstat {input} > {output} 2> {log}
        #~ """

#~ rule call_peaks:
    #~ input: control = "04aln/{control_id}.sorted.bam", case="04aln/{case_id}.sorted.bam"
    #~ output: bed="05peak/{case_id}_vs_{control_id}_peaks.bed"
    #~ log: "00log/{case_id}_vs_{control_id}_call_peaks.log"
    #~ params:
        #~ name = "{case_id}_vs_{control_id}"
    #~ resources:
        #~ mem  = 4,
        #~ time = 30

    #~ message: "call_peaks macs14 {input}: {threads} threads / {resources.mem}"
    #~ shell:
        #~ """
        #~ module load macs14 R
        #~ macs14 -t {input.case} \
            #~ -c {input.control} -f BAM -g {config[macs_g]} \
            #~ --outdir 04peak -n {params.name} -p 10e-5 &> {log}
            #~ cd 05peak && Rscript {params.name}_model.r
        #~ """

#~ rule merge_controls:
    #~ input:   expand("03aln/{sample}.sorted.bam", sample = config["controls"])
    #~ output:  "03aln/control.sorted.bam"
    #~ log:     "00log/merge_controls"
    #~ threads: 2
    #~ params:
        #~ mem = "1G",
        #~ time = "15"
    #~ message: "merge_controls {input}: {threads} threads / {params.mem}"
    #~ shell:
        #~ """
        #~ inbam=( {input} )
        #~ if [[ ${{#inbam[@]}} -eq 1 ]]; then 
            #~ ln -s $(cd $(dirname {input}) && pwd)/$(basename {input}) {output}
        #~ else
            #~ module load samtools/1.2
            #~ samtools merge -r -@{threads} {output} {input}
        #~ fi
        #~ """

#~ rule flagstat_bam:
    #~ input:  "03aln/{sample}.sorted.bam"
    #~ output: "03aln/{sample}.sorted.bam.flagstat"
    #~ log:    "00log/{sample}.flagstat_bam"
    #~ threads: 1
    #~ params:
        #~ mem   = "500M",
        #~ time  = "10"
    #~ message: "flagstat_bam {input}: {threads} threads / {params.mem}"
    #~ shell:
        #~ """
        #~ module load samtools/1.2
        #~ samtools flagstat {input} > {output}
        #~ """

#~ rule quantify_genes:
    #~ input:
        #~ genome = 'genome.fa',
        #~ r1 = 'fastq/{sample}.R1.fastq.gz',
        #~ r2 = 'fastq/{sample}.R2.fastq.gz'
    #~ output:
        #~ '{sample}.txt'
    #~ shell:
        #~ 'echo {input.genome} {input.r1} {input.r2} > {output}'



#~ rule find_narrow_peaks:
    #~ input:  "03aln/{sample}.sorted.bam", "03aln/control.sorted.bam"
    #~ output: "04peak/{sample}_model.r", "04peak/{sample}_peaks.narrowPeak",
            #~ "04peak/{sample}_peaks.xls", "04peak/{sample}_summits.bed",
            #~ "04peak/{sample}_model.pdf"
    #~ log:    "00log/{sample}.find_narrow_peaks"
    #~ threads: 2
    #~ params:
        #~ mem  = "2G",
        #~ time = "15"
    #~ message: "find_narrow_peaks {input}: {threads} threads / {params.mem}"
    #~ shell:
        #~ """
        #~ module load macs/2.1.0.20150420 R
        #~ macs2 callpeak -t {input[0]} \
                #~ -c {input[1]} -f BAM -g {config[macs_g]} \
                #~ --outdir 04peak -n {wildcards.sample} -q 0.001 &> {log}
        #~ cd 04peak && Rscript {wildcards.sample}_model.r
        #~ """
