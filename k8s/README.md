# Running nf_funannotate1 on Kubernetes (NRP Nautilus)

`-profile annotate,nrp` runs every task as its own pod through Nextflow's
built-in k8s executor. No SLURM layer is involved. Nextflow itself runs as a
**Job** inside the cluster, one per run: the Job's command is the Nextflow run
and it exits when the pipeline does (NRP prohibits idle `sleep infinity` pods;
see "NRP fair use" below). The run Job and every task pod mount one shared
ReadWriteMany PVC at `/data`, which holds the source checkout, launch dirs,
`work/`, reference DBs and results.

This repository holds only site-neutral pieces. Your namespace, NRP project
and bucket names belong in your own (ideally private) deployment repo: a
kustomize overlay for the cluster objects and a small Nextflow config for the
run. `k8s/overlays/example/` is the template.

Config pieces:

| File | Role |
|---|---|
| `conf/provision_singularity.config` | label → image map (reused) |
| `conf/executor_k8s.config` | turns that map into pod-runnable OCI images; no nested apptainer (GeneMark / prodigal / MariaDB setup run inside their image); pods run as root; node-local predict |
| `conf/site_nrp.config` | what is true for any NRP namespace: 16 core / 32 GB task caps, `opportunistic` priority + no GPU nodes for task pods, PVC and `/data/refdb/*` layout, eggNOG mount, Ceph S3 endpoint |
| your site config (`-c`) | `params.k8s_namespace` (required), `params.nrp_project`, `executor.queueSize` (concurrent task pods), anything else site-specific |
| `k8s/base/` | kustomize base: PVC, ServiceAccount + Role (pods, configmaps, jobs) |
| `k8s/run/` | the Nextflow run Job (one per run directory) |
| `k8s/tools/` | finite Jobs `s3-sync`, `build-rust-tools`; `shell`, a 1 h setup/inspection pod |
| `k8s/overlays/example/` | fill-in-the-blank overlay (namespace, NRP project annotation), `runs/my-run/`, `tools/*`, `site.config` |
| `k8s/params_nrp_test.yaml` | first smoke-test toggles |

Always list the pipeline profile first (`annotate,nrp`), and use Nextflow
**26.x**. Its config parser applies profiles in command-line order, which lets
the k8s overrides win over `annotate`'s defaults. 25.10 applies them in
file-definition order instead, and the run would then try to use apptainer.

## One-time setup

1. **Overlay.** Copy `k8s/overlays/example/` into your deployment repo and fill
   in `namespace` and the NRP project. To track this repo instead of copying
   the base, point `resources` at it by URL, pinned to a tag or commit:
   `https://github.com/stajichlab/nf_funannotate1//k8s/base?ref=<tag-or-sha>`.
2. **S3 credentials** (optional, for `s3://` inputs/outputs and the staging
   pod), created directly in the cluster, never committed:
   ```bash
   kubectl create secret generic nrp-s3-creds -n <namespace> \
       --from-literal=AWS_ACCESS_KEY_ID=... --from-literal=AWS_SECRET_ACCESS_KEY=...
   ```
3. **Site config** for the run: copy `k8s/overlays/example/site.config` next
   to your overlay and fill in the same namespace and NRP project. It sets the
   `k8s {}` scope directly: `-c` files are read after the profiles, so params
   set there alone don't reach `conf/site_nrp.config`.
4. **Apply** (the RBAC part needs a namespace admin):
   ```bash
   kubectl apply -k <overlay dir>
   ```

**Setting up the PVC** uses the short-lived shell (`kubectl apply -k <overlay
dir>/tools/shell`, then `kubectl exec -it -n <namespace> nf-shell -- bash`;
delete it when done). Setup and inspection only, never computation:

```bash
# 1. Source checkout on the PVC (task pods read bin/ and assets/ from projectDir)
mkdir -p /data/src /data/refdb /data/runs
git clone https://github.com/stajichlab/nf_funannotate1 /data/src/nf_funannotate1
# ...and copy your site.config onto the PVC as /data/src/site.config (from your
#    workstation: kubectl exec -i -n <ns> nf-shell -- sh -c 'cat > /data/src/site.config' < site.config)

# 2. Reference files the pipeline can't download itself (e.g. from your S3
#    bucket with the s3-sync Job, or `kubectl exec -i ... -- sh -c 'cat > f' < f`):
#      /data/refdb/swissprot_fungi.faa   (params.proteins)
#      /data/refdb/template.sbt          (params.sbt_template)
#      /data/refdb/busco_lineages/       (odb10 lineages, for BUSCO steps)
#      /data/refdb/eggnog_db/            (only for --run_annotate)
#    funannotate_db and taxondb are built by SETUP_FUNANNOTATE_DB /
#    SETUP_TAXONDB on the first run (storeDir-cached under /data/refdb); see
#    "Reference DB archive" below to skip the ~2 h funannotate_db build.

# 3. A run directory: samples.csv (+ genomes) and a params file
mkdir -p /data/runs/nrp_test && cd /data/runs/nrp_test
cp -r /data/src/nf_funannotate1/samples.csv /data/src/nf_funannotate1/test_run .
cp /data/src/nf_funannotate1/k8s/params_nrp_test.yaml params.yaml
```

**Running** is a Job per run directory. Copy `overlays/example/runs/my-run/`,
set its name suffix, `RUN_DIR` and params file, then from your workstation:

```bash
kubectl apply  -k <overlay dir>/runs/nrp_test       # start, or resume after a stop
kubectl logs -f -n <namespace> job/nextflow-nrp_test  # follow (also nf_run.log in RUN_DIR)
kubectl delete -k <overlay dir>/runs/nrp_test       # stop; apply again to -resume
```

The Job always runs with `-resume` and keeps the previous `nf_run.log`
(timestamped). Watch task pods with `kubectl get pods -n <namespace> -w`.
Failed task pods are kept, so `kubectl logs <pod>` shows them.

**Resume cache.** Nextflow's cache (`.nextflow`: the LevelDB task cache and
run history) runs on the pod's local disk, not the PVC, with LevelDB's
memory-mapped writes turned off (`-Dleveldb.mmap=false`). `k8s/run/nf-run.sh`
copies it to `<RUN_DIR>/.nextflow-snapshots/snap-<time>` every
`CACHE_SYNC_SECONDS` (default 300) and once more when Nextflow exits, keeping
the newest `CACHE_SNAPSHOTS_KEEP` (3). Each snapshot is written to a `.tmp-`
directory, fsynced and renamed, so a `snap-` directory is always complete. The
next Job restores the newest one. Losing the node mid-run costs at most the
tasks finished since the last snapshot. With the cache on CephFS, a lost node
left the database corrupt and every retry failed. A snapshot that won't open is
renamed `bad-*`, and the Job's retry falls back to the one before. After
`kubectl delete`, the old pod can still be stopping when the new Job starts, so
the new pod waits for the old one's heartbeat (`RUN_DIR/.nf-run.heartbeat`) to
go stale before restoring. To stop a run and be sure its pod is gone, use
`kubectl delete -k <run> --cascade=foreground`. A run
started before this change has its `RUN_DIR/.nextflow` migrated on the first
launch and renamed to `.nextflow.migrated-<time>`.

**Helper Jobs** (`k8s/tools/`) go through the overlay too: `s3-sync` (set
`SRC` / `DST`, runs `aws s3 sync`, exits) and `build-rust-tools` (builds the
SRA Rust helpers into the checkout's `tools/bin`; only needed for host-tool
profiles or `sra_tools` older than 1.4.0).

## Turning features back on

In roughly this order, one per run (`-resume` keeps what already passed):

GeneMark is already on in the smoke test. GENEMARK_RUN runs in
`teambraker/braker3`, which needs no license key (`--gm_key` exists only for
a key-gated `--container_genemark` image). With RNA-seq it tries ET mode and
falls back to ES when ET can't train (intron-poor genomes).

1. `--run_sra_fetch true`: RNA-seq fetch plus Trinity/PASA training. PASA's
   MariaDB starts inside the funannotate image (`mariadbd` is bundled). If that
   misbehaves, `--pasa_mysql false` falls back to SQLite.
2. `--run_annotate true`: needs `/data/refdb/eggnog_db`, mounted into the pod
   at `/opt/databases/eggnog_db`.

## Not available on k8s yet

- **SignalP, DeepTMHMM, antiSMASH**: the images are locally built `.sif` files
  (SignalP and DeepTMHMM are licensed). They need pushing to a registry the
  pods can pull from. For the private ones, that means a private registry
  (on NRP, its GitLab registry) plus a `pod imagePullSecret`. GPU runs would
  use the `accelerator` directive in place of the SLURM `--gres` options.
- **FCS-GX purge** (GENOME_CLEAN without `--skip_fcs`): needs the ~470 GB
  database in `/dev/shm`.
- **prodigal** (`--run_prodigal`): host mode needs prodigal and python in one
  image. The funannotate image is the default and has not been checked for
  prodigal.

## NRP fair use and limits

From NRP's docs ([Jobs](https://nrp.ai/documentation/userdocs/running/jobs/),
[CPU-only](https://nrp.ai/documentation/userdocs/running/cpu-only/)):

- **No idle / interactive pods.** "Running in interactive mode (`sleep infinity`
  command and manual start of computation) or any command that doesn't end by
  itself is prohibited, and user can be banned." Hence the run Job and finite
  helper Jobs; `tools/shell` is only for brief setup and ends by itself in 1 h.
- **No fair queue.** "If you submit 1000 jobs, you block all other users."
  Nextflow's `executor.queueSize` is the bounded queue: set it in your site
  config (the profile default is 100; a few dozen is considerate).
- **CPU-only work** should be preemptible and stay off GPU nodes:
  `conf/site_nrp.config` gives every task pod `priorityClassName: opportunistic`
  and a node anti-affinity on `feature.node.kubernetes.io/pci-10de.present`. A
  preempted task fails and Nextflow retries it.
- **Pods without a controller** (the k8s executor's default) are capped at
  **16 cores / 32 GB** and a **6 h lifetime**; `conf/site_nrp.config` caps tasks
  at 16 / 32 GB. For bigger or longer tasks run them as Jobs: the base Role
  allows `batch/jobs`; set `k8s.computeResourceType = 'Job'` and raise
  `k8s_max_cpus` / `k8s_max_memory` in your site config.

## Reference DB archive

A fully built `/data/refdb/funannotate_db` is ~33 GB, ~9.5 GB as `.tar.gz`.
Archive it once to your bucket and restore it into new PVCs to skip the ~2 h
SETUP_FUNANNOTATE_DB build (odb10 lineages for predict's BUSCO training go in
it too):

```bash
# from a pod with the PVC at /data, S3 creds and tar (e.g. amazon/aws-cli + `yum install -y tar gzip pigz`)
EP=--endpoint-url=https://s3-west.nrp-nautilus.io
tar -cf - -C /data/refdb funannotate_db | pigz | aws $EP s3 cp - s3://<bucket>/refdb/funannotate_db.tar.gz --expected-size <bytes>
aws $EP s3 cp s3://<bucket>/refdb/funannotate_db.tar.gz - | tar -xzf - -C /data/refdb
```

## Known risks / untested

- Pods get the namespace default ephemeral storage for `/tmp` (50 Gi on
  NRP). Predict and Trinity scratch live there; very large genomes may need
  Jobs with an explicit `disk` request.
- `k8s_max_cpus` / `k8s_max_memory` cap every task (16 / 32 GB on NRP). Raise
  them only together with Jobs.

## GeneMark-only run

`genemark_sidecar.nf` runs GeneMark-ES (or -ET when `--rnaseq_reads_dir` is set) once per genome and
stops. On NRP use it through the same run Job, with these changes to the Job's environment:

| Variable | Value |
|---|---|
| `PROFILE` | `genemark_sidecar,nrp` (pipeline profile first) |
| `NF_EXTRA_ARGS` | `-main-script genemark_sidecar.nf` |
| `PARAMS_FILE` | a YAML with `samples`, `target`, `publish_per_genome: true`, `genemark_force_container: false` |

- `genemark_force_container: false` is required on k8s. The default (`true`) forces an `apptainer exec`
  that does not exist inside a pod.
- `publish_per_genome: true` writes `<target>/<out>/<out>.genemark.{gtf,mod}` and an empty
  `<out>.other.gff3`. The default writes the flat layout that `--genemark_sidecar_dir` reads.
- An empty `.gtf` means GeneMark skipped the genome (pre-flight or "too small"). There is no `.mod` then.
- First-attempt task request is 4 CPUs / 16 GB (8 / 32 GB on retry), from `conf/site_nrp.config`.
