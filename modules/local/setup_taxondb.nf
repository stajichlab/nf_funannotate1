process SETUP_TAXONDB {
    label 'setup'
    label 'process_single'

    storeDir params.taxondb

    output:
    path "names.dmp",    emit: ready
    path "nodes.dmp"
    path "merged.dmp"
    path "delnodes.dmp"
    path "division.dmp"
    path "gencode.dmp"
    path "citations.dmp"

    script:
    """
    set -euo pipefail
    wget --no-verbose https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar zxf taxdump.tar.gz
    rm taxdump.tar.gz
    """

    stub:
    """
    for f in names.dmp nodes.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp citations.dmp; do
        touch \$f
    done
    """
}
