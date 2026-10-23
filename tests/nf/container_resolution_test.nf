// Test harness: print the container each process label resolves to. Driven by
// tests/test_sif_dir_resolution.sh with the apptainer engine switched off; the tasks need not
// succeed (on macOS they do not), the driver reads the container from the echoed command.
process C_FUNANNOTATE {
    label 'funannotate'

    output:
    stdout

    script:
    """
    echo funannotate=${task.container}
    """
}

process C_GENOME_CLEAN {
    label 'genome_clean'

    output:
    stdout

    script:
    """
    echo genome_clean=${task.container}
    """
}

process C_ANTISMASH {
    label 'antismash'

    output:
    stdout

    script:
    """
    echo antismash=${task.container}
    """
}

process C_SIGNALP {
    label 'signalp'

    output:
    stdout

    script:
    """
    echo signalp=${task.container}
    """
}

process C_DEEPTMHMM {
    label 'deeptmhmm'

    output:
    stdout

    script:
    """
    echo deeptmhmm=${task.container}
    """
}

workflow {
    C_FUNANNOTATE(); C_GENOME_CLEAN(); C_ANTISMASH(); C_SIGNALP(); C_DEEPTMHMM()
    // Where tasks run (Linux CI) the echoed line is the output; where they fail (macOS) the error
    // report quotes the command, which the driver also reads.
    C_FUNANNOTATE.out.mix(C_GENOME_CLEAN.out, C_ANTISMASH.out, C_SIGNALP.out, C_DEEPTMHMM.out)
        .view { it.trim() }
}
