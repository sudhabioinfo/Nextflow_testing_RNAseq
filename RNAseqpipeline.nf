params.dir = "${projectDir}"  // Project directory
params.fastq = "${projectDir}/data/*.fastq"  // Support multiple FASTQ files
params.QC_reports = "${projectDir}/results/QC_reports"  // Directory for FastQC reports
params.trim = "${projectDir}/results/trimmed"  // Directory for trimmed FASTQ files
params.trim_QC = "${projectDir}/results/trimmed_QC_reports"  // Directory for trimmed FastQC reports
params.trimfile = "${projectDir}/results/trimmed/trimmed_{name}.fastq"  // Template for trimmed FASTQ files
params.Hisat = "${projectDir}/HISAT2/grch38/genome"  // Path to HISAT2 genome index
params.Hisataligned = "${projectDir}/results/hisat2_aligned"  // outputs
params.FC = "${projectDir}/results/feature_counts"
params.gtf = "${projectDir}/Homo_sapiens.GRCh38.106.gtf.gz"

process QC{

    publishDir("${params.QC_reports}", mode: 'copy')

    input:
    path fastq

    output:
    path "*"

    script:
    """
    fastqc $fastq
    """
    }

process trim {

    publishDir params.trim, mode: 'copy'

    input:
    path fastq

    output:
    path "*"

    script:
    """
    java -jar /Users/sushant/Desktop/Github/nextflow/nextflowprojects/RNAseq_pipelines/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 $fastq trimmed.fastq TRAILING:10 -phred33
    echo "Trimmomatic Finished running!"
    """
    }

process trimQC{

    publishDir("${params.trim_QC}", mode: 'copy')

    input:
    path trimfile

    output:
    path "*"

    script:
    """
    fastqc $trimfile 
    """
    }

process hisat {

    publishDir("${params.Hisataligned}", mode: 'copy')

    input:
    path trimmed_fastq

    output:
    path "trimmed.bam"
    path "hisat2.log"

    script:
    """
    hisat2 -q --rna-strandness R -x $params.Hisat -U $trimmed_fastq 2> >(tee hisat2.log >&2) | \
    samtools sort -o trimmed.bam
    """
}

process featurecounts{

    publishDir("${params.FC}", mode: 'copy')

    input:
    path bam
    path gtf

    output:
    path "*"

    script:
    """
    # Decompress GTF if it is gzipped
    if [[ $gtf == *.gz ]]; then
        gunzip -c $gtf > annotation.gtf
    else
        cp $gtf annotation.gtf
    fi
    featureCounts -S 2 -a annotation.gtf -o feature.counts.txt -T 4 $bam 

    """

    }

workflow {

    fastq_ch=Channel.fromPath(params.fastq)
    QC(fastq_ch)
    QC.out.view()

    trim(fastq_ch)
    trim.out.view()

    trimqc_ch = Channel.fromPath(params.trimfile)
    trimQC(trimqc_ch)
    trimQC.out.view()

    hisat(trim.out)
    featurecounts(hisat.out[0], params.gtf)

    }


