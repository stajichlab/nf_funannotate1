#!/usr/bin/env bash
# k8s_avoid_nodes (conf/site_nrp.config) adds a hostname NotIn term to every task pod's node
# affinity, and only when set. Usage: bash tests/test_k8s_avoid_nodes.sh   (needs nextflow)
set -uo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
command -v nextflow >/dev/null || { echo "SKIP: nextflow not found"; exit 0; }
# The nrp profile needs Nextflow >= 26 (k8s/README.md): 25.10 applies profiles in file order and then fails
# with "Cannot enable more than one container engine".
major=$(nextflow -version 2>/dev/null | sed -n 's/.*version \([0-9][0-9]*\)\..*/\1/p' | head -1)
[ "${major:-0}" -ge 26 ] || { echo "SKIP: needs Nextflow >= 26 (found ${major:-unknown})"; exit 0; }
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
pass=0; fail=0
check() { if eval "$2"; then pass=$((pass+1)); echo "ok   $1"; else fail=$((fail+1)); echo "FAIL $1"; fi; }
pod() { ( cd "$T" && nextflow -c "$PROJECT_DIR/nextflow.config" run "$PROJECT_DIR/tests/nf/print_task_pod.nf" \
          -profile annotate,nrp --k8s_namespace x -ansi-log false "$@" 2>&1 ) | grep '^POD=' ; }
printf 'k8s_avoid_nodes: [bad-node.example.org, other-node]\n' > "$T/p.yaml"

default=$(pod)
check "default: GPU-node exclusion and opportunistic priority are kept" 'echo "$default" | grep -q "pci-10de.present" && echo "$default" | grep -q "priorityClassName:opportunistic"'
check "default: no hostname term"                                       '! echo "$default" | grep -q "kubernetes.io/hostname"'

withnodes=$(pod -params-file "$T/p.yaml")
check "params file: hostname NotIn term lists both nodes"               'echo "$withnodes" | grep -q "kubernetes.io/hostname, operator:NotIn, values:\[bad-node.example.org, other-node\]"'
check "params file: GPU-node exclusion still there, in the same term"   'echo "$withnodes" | grep -q "pci-10de.present, operator:NotIn, values:\[true\]\], \[key:kubernetes.io/hostname"'

echo "passed=$pass failed=$fail"; [ "$fail" -eq 0 ]
