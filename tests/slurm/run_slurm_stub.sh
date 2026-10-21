#!/usr/bin/env bash
# Stub-run the pipeline against a throwaway single-node SLURM (stock config, no
# site partitions) to check that `-profile test,slurm` submits and completes jobs
# on a plain cluster. Needs docker. Usage: bash tests/slurm/run_slurm_stub.sh
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
NAME=nf_slurm_test
NXF_VER="${NXF_VER:-25.10.0}"

trap 'docker rm -f "$NAME" >/dev/null 2>&1 || true' EXIT
docker rm -f "$NAME" >/dev/null 2>&1 || true

# Repo is mounted at the same path inside so the work dir is "shared" with the
# (same-host) compute node.
docker run --privileged --cgroupns=host --name "$NAME" -v "$REPO:$REPO" -w "$REPO" ubuntu:24.04 bash -euxc "
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq
  apt-get install -y -qq slurm-wlm munge default-jre-headless curl procps pigz >/dev/null

  mkdir -p /var/spool/slurmd /var/spool/slurmctld /var/log/slurm /run/munge
  chown slurm: /var/spool/slurmd /var/spool/slurmctld /var/log/slurm
  # Run munged/slurm as root-owned single-node setup; SlurmUser=root keeps it simple.
  dd if=/dev/urandom bs=1 count=1024 of=/etc/munge/munge.key 2>/dev/null
  chown munge: /etc/munge/munge.key && chmod 400 /etc/munge/munge.key
  chown -R munge: /run/munge /var/log/munge /var/lib/munge 2>/dev/null || true
  su -s /bin/sh munge -c 'munged'

  H=\$(hostname)
  cat > /etc/slurm/slurm.conf <<CONF
ClusterName=ci
SlurmctldHost=\$H
SlurmUser=root
StateSaveLocation=/var/spool/slurmctld
SlurmdSpoolDir=/var/spool/slurmd
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdLogFile=/var/log/slurm/slurmd.log
ProctrackType=proctrack/linuxproc
TaskPlugin=task/none
SelectType=select/cons_tres
SelectTypeParameters=CR_Core
NodeName=\$H CPUs=4 RealMemory=4000 State=UNKNOWN
SlurmdParameters=config_overrides
PartitionName=main Nodes=\$H Default=YES MaxTime=INFINITE State=UP
CONF
  printf 'CgroupPlugin=cgroup/v2\\nIgnoreSystemd=yes\\n' > /etc/slurm/cgroup.conf
  # slurmd's cgroup/v2 plugin wants <parent of own cgroup>/system.slice; absent in some containers.
  CG=/sys/fs/cgroup\$(cut -d: -f3 /proc/self/cgroup); mkdir -p \"\$(dirname \$CG)/system.slice\" || true
  trap 'tail -n 30 /var/log/slurm/*.log' EXIT
  slurmctld
  slurmd
  sleep 5
  sinfo

  (cd /usr/local/bin && curl -fsSL https://get.nextflow.io | NXF_VER=$NXF_VER bash)
  nextflow -c tests/slurm/stock_slurm.config run . -profile test,slurm -stub-run
  grep -c 'executor > slurm' .nextflow.log || true
  sacct -X --format=JobID,JobName%30,State 2>/dev/null | head -20 || squeue -h | head
"
