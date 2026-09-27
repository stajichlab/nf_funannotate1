# Running nf_funannotate1 on Kubernetes (NRP Nautilus)

`-profile annotate,nrp` runs every task as its own pod through Nextflow's
built-in k8s executor. No SLURM layer is involved. Nextflow itself runs in a
long-lived **head pod** (a Deployment) inside the cluster. The head pod and every task pod
mount one shared ReadWriteMany PVC at `/data`, which holds the source checkout,
launch dirs, `work/`, reference DBs and results.

This repository holds only site-neutral pieces. Your namespace, NRP project
and bucket names belong in your own (ideally private) deployment repo: a
kustomize overlay for the cluster objects and a small Nextflow config for the
run. `k8s/overlays/example/` is the template.

Config pieces:

| File | Role |
|---|---|
| `conf/provision_singularity.config` | label → image map (reused) |
| `conf/executor_k8s.config` | turns that map into pod-runnable OCI images; no nested apptainer (GeneMark / prodigal / MariaDB setup run inside their image); pods run as root; node-local predict |
| `conf/site_nrp.config` | what is true for any NRP namespace: 16 core / 32 GB task caps, PVC and `/data/refdb/*` layout, eggNOG mount, Ceph S3 endpoint |
| your site config (`-c`) | `params.k8s_namespace` (required), `params.nrp_project`, anything else site-specific |
| `k8s/base/` | kustomize base: PVC, ServiceAccount + Role (pods, configmaps, jobs), head Deployment |
| `k8s/tools/` | one-off pods: `build-rust-tools`, `s3-stage` |
| `k8s/overlays/example/` | fill-in-the-blank overlay (namespace, NRP project annotation) |
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
   kubectl exec -it -n <namespace> deploy/funannotate-nextflow-head -- bash
   ```

Inside the head pod:

```bash
# 1. Source checkout on the PVC (task pods read bin/ and assets/ from projectDir)
mkdir -p /data/src /data/refdb /data/runs
git clone https://github.com/stajichlab/nf_funannotate1 /data/src/nf_funannotate1
# ...and copy your site.config somewhere on the PVC, e.g. /data/src/site.config

# 2. Reference files the pipeline can't download itself (e.g. from your S3
#    bucket with the s3-stage pod, or `kubectl exec -i ... -- sh -c 'cat > f' < f`):
#      /data/refdb/swissprot_fungi.faa   (params.proteins)
#      /data/refdb/template.sbt          (params.sbt_template)
#      /data/refdb/busco_lineages/       (odb10 lineages, for BUSCO steps)
#      /data/refdb/eggnog_db/            (only for --run_annotate)
#    funannotate_db and taxondb are built by SETUP_FUNANNOTATE_DB /
#    SETUP_TAXONDB on the first run (storeDir-cached under /data/refdb); see
#    "Reference DB archive" below to skip the ~2 h funannotate_db build.

# 3. Smoke test
mkdir -p /data/runs/nrp_test && cd /data/runs/nrp_test
cp -r /data/src/nf_funannotate1/samples.csv /data/src/nf_funannotate1/test_run .
nextflow run /data/src/nf_funannotate1 -profile annotate,nrp -c /data/src/site.config \
    -params-file /data/src/nf_funannotate1/k8s/params_nrp_test.yaml
```

Launch long runs detached (`setsid nohup nextflow run ... &`, or tmux) so a
dropped `kubectl exec` session doesn't end them. Watch pods with
`kubectl get pods -n <namespace> -w`. Failed task pods are kept, so
`kubectl logs <pod>` shows them. Successful ones are cleaned up.

**One-off pods** (`k8s/tools/`) go through the overlay too, so they get your
namespace. The example has `tools/build-rust-tools` and `tools/s3-stage`:
```bash
kubectl apply -k <overlay dir>/tools/s3-stage
kubectl delete pod -n <namespace> funannotate-s3-stage   # when done
```
`build-rust-tools` builds the SRA Rust helpers into `tools/bin` of the PVC
checkout. It's only needed with `sra_tools` images older than 1.4.0; the
container profiles use 1.4.0, which ships them.

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

## NRP limits on pods

NRP's admission policy treats pods without a controller (which is what the
k8s executor creates by default) specially: at most **16 cores / 32 GB**, and
a **6 h lifetime** (activeDeadline). `conf/site_nrp.config` caps tasks at
16 / 32 GB accordingly, and the head runs as a Deployment so it isn't killed
at 6 h. For bigger or longer tasks (Trinity, large-genome predict), run tasks
as Jobs instead: the base Role already allows `batch/jobs`; set
`k8s.computeResourceType = 'Job'` and raise `k8s_max_cpus` / `k8s_max_memory`
in your site config.

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
