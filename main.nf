#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { plink2 } from './modules/plink.nf'

Channel
    .fromPath(params.fam)
    .set{fam}

Channel
    .fromPath(params.ids_to_include)
    .set{ids_to_include}

Channel
    .fromPath(params.phenotypes_to_include)
    .set{phenotypes_to_include}

Channel
    .fromPath(params.genotypes)
    .map { file -> tuple(file.baseName, file) }
    .groupTuple(by:0)
    .combine(fam)
    .combine(ids_to_include)
    .dump()
    .set{genotypes_to_filter}

workflow {
    plink2(params)
}