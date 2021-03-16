process bolt_lmm {
    label "bolt"
    tag "${prefix}"

    input:
    output:
    script:
        """
        bolt --numThreads "${task.cpus}" --bfile "${prefix}" \
            --fam "${fam}" --exclude LMAO
        """
}