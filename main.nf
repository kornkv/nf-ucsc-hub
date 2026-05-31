#!/usr/bin/env nextflow

/*
 * nf-ucsc-hub: Build UCSC track hubs from genome assemblies
 */

nextflow.enable.dsl=2

params.genome   = null
params.gtf      = null
params.trf      = null
params.assembly = null
params.outdir   = './hub'
params.email    = 'your@email.cz'

if (!params.genome || !params.gtf || !params.trf) {
    error """
    Required parameters missing. Usage:
      nextflow run main.nf --genome genome.fa --gtf genes.gtf --trf trf.dat --outdir ./hub
    """.stripIndent()
}

def assembly = params.assembly ?: file(params.genome).baseName

log.info """
    ========================
    nf-ucsc-hub
    ========================
    genome   : ${params.genome}
    gtf      : ${params.gtf}
    trf      : ${params.trf}
    assembly : ${assembly}
    outdir   : ${params.outdir}
    ========================
""".stripIndent()

// ---- Processes ----

process MAKE_2BIT {
    input:
    path genome

    output:
    path "${assembly}.2bit", emit: twobit

    script:
    """
    faToTwoBit ${genome} ${assembly}.2bit
    """
}

process CHROM_SIZES {
    input:
    path twobit

    output:
    path "chrom.sizes", emit: sizes

    script:
    """
    twoBitInfo ${twobit} stdout | sort -k2rn > chrom.sizes
    """
}

process TRF_TO_BED {
    input:
    path dat

    output:
    path "trf.sorted.bed", emit: bed

    script:
    """
    awk '/^Sequence:/{chr=\$2} NF>=15 && \$1~/^[0-9]+\$/{print chr, \$1-1, \$2}' OFS='\\t' ${dat} \
        | sort -k1,1 -k2,2n > trf.sorted.bed
    """
}

process TRF_TO_BIGBED {
    input:
    path bed
    path sizes

    output:
    path "trf.bb", emit: bb

    script:
    """
    bedToBigBed ${bed} ${sizes} trf.bb
    """
}

process GTF_TO_BED {
    input:
    path gtf

    output:
    path "genes.sorted.bed", emit: bed

    script:
    """
    gtfToGenePred -ignoreGroupsWithoutExons ${gtf} genes.genePred
    genePredToBed genes.genePred genes.bed
    sort -k1,1 -k2,2n genes.bed > genes.sorted.bed
    """
}

process GENES_TO_BIGBED {
    input:
    path bed
    path sizes

    output:
    path "genes.bb", emit: bb

    script:
    """
    bedToBigBed -type=bed12 ${bed} ${sizes} genes.bb
    """
}

process SEARCH_INDEX {
    input:
    path bed

    output:
    path "genes.ix",  emit: ix
    path "genes.ixx", emit: ixx

    script:
    """
    awk -F'\\t' '{print \$4"\\t"\$4}' ${bed} | sort -u > names.txt
    ixIxx names.txt genes.ix genes.ixx
    """
}

process ASSEMBLE_HUB {
    publishDir params.outdir, mode: 'copy'

    input:
    path twobit
    path sizes
    path trf_bb
    path genes_bb
    path ix
    path ixx

    output:
    path "**"

    script:
    """
    mkdir -p ${assembly}

    DEFPOS=\$(head -1 ${sizes} | awk '{print \$1":1-"(\$2<100000?\$2:100000)}')

    cat > hub.txt << 'HUBEOF'
hub ${assembly}
shortLabel ${assembly}
longLabel ${assembly} genome hub
genomesFile genomes.txt
email ${params.email}
HUBEOF

    cat > genomes.txt << GEOF
genome ${assembly}
trackDb ${assembly}/trackDb.txt
twoBitPath ${assembly}/${assembly}.2bit
organism ${assembly}
defaultPos \${DEFPOS}
GEOF

    cat > ${assembly}/trackDb.txt << 'TEOF'
track genes
bigDataUrl genes.bb
shortLabel Genes
longLabel Gene annotations
type bigBed 12
searchIndex name
searchTrix genes.ix
visibility pack

track trf
bigDataUrl trf.bb
shortLabel TRF
longLabel Tandem Repeats
type bigBed 3
visibility dense
TEOF

    cp ${twobit} ${assembly}/
    cp ${trf_bb} ${assembly}/
    cp ${genes_bb} ${assembly}/
    cp ${ix} ${assembly}/
    cp ${ixx} ${assembly}/
    """
}

// ---- Workflow ----

workflow {
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
    gtf_ch    = Channel.fromPath(params.gtf, checkIfExists: true)
    trf_ch    = Channel.fromPath(params.trf, checkIfExists: true)

    MAKE_2BIT(genome_ch)
    CHROM_SIZES(MAKE_2BIT.out.twobit)

    TRF_TO_BED(trf_ch)
    TRF_TO_BIGBED(TRF_TO_BED.out.bed, CHROM_SIZES.out.sizes)

    GTF_TO_BED(gtf_ch)
    GENES_TO_BIGBED(GTF_TO_BED.out.bed, CHROM_SIZES.out.sizes)
    SEARCH_INDEX(GTF_TO_BED.out.bed)

    ASSEMBLE_HUB(
        MAKE_2BIT.out.twobit,
        CHROM_SIZES.out.sizes,
        TRF_TO_BIGBED.out.bb,
        GENES_TO_BIGBED.out.bb,
        SEARCH_INDEX.out.ix,
        SEARCH_INDEX.out.ixx
    )
}
