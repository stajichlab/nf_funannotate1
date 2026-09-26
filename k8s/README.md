# Running nf_funannotate1 on Kubernetes (NRP Nautilus)

`-profile annotate,nrp` runs every task as its own pod through Nextflow's
built-in k8s executor. No SLURM layer is involved. Nextflow itself runs in a
long-lived **head pod** inside the cluster. The head pod and every task pod
mount one shared ReadWriteMany PVC at `/data`, which holds the source checkout,
launch dirs, `work/`, reference DBs and results.

Config pieces:

| File | Role |
|---|---|
| `conf/provision_singularity.config` | label → image map (reused) |
| `conf/executor_k8s.config` | turns that map into pod-runnable OCI images; no nested apptainer (GeneMark / prodigal / MariaDB setup run inside their image); pods run as root |
| `conf/site_nrp.config` | namespace, PVC, service account, `/data/refdb/*` paths, eggNOG mount, Ceph S3 endpoint |
| `k8s/params_nrp_test.yaml` | first smoke-test toggles |

Always list the pipeline profile first (`annotate,nrp`), and use Nextflow
**26.x**. Its config parser applies profiles in command-line order, which lets
the k8s overrides win over `annotate`'s defaults. 25.10 applies them in
file-definition order instead, and the run would then try to use apptainer.

## One-time setup

The `ucr-stajichlab` namespace already has the `nextflow-runner`
ServiceAccount/RBAC and the `nrp-s3-creds` Secret from the BFD project
(`../BFD/k8s/rbac.yaml`), so they're reused here.

```bash
kubectl apply -f k8s/pvc.yaml
kubectl apply -f k8s/head-pod.yaml
kubectl exec -it -n ucr-stajichlab funannotate-nextflow-head -- bash
```

Inside the head pod:

```bash
# 1. Source checkout on the PVC (task pods read bin/ and assets/ from projectDir)
mkdir -p /data/src /data/refdb /data/runs
git clone https://github.com/stajichlab/nf_funannotate1 /data/src/nf_funannotate1

# 2. Reference files the pipeline can't download itself (from UCR HPCC,
#    e.g. via s3://stajichlab/... with the aws{} endpoint in site_nrp.config,
#    or `kubectl cp` from your workstation):
#      /data/refdb/swissprot_fungi.faa   (params.proteins)
#      /data/refdb/template.sbt          (params.sbt_template)
#      /data/refdb/busco_lineages/       (BUSCO v10 lineage tree, for BUSCO steps)
#      /data/refdb/eggnog_db/            (only for --run_annotate)
#    funannotate_db and taxondb are built by SETUP_FUNANNOTATE_DB /
#    SETUP_TAXONDB on the first run (storeDir-cached under /data/refdb).

# 3. Smoke test
mkdir -p /data/runs/nrp_test && cd /data/runs/nrp_test
cp -r /data/src/nf_funannotate1/samples.csv /data/src/nf_funannotate1/test_run .
nextflow run /data/src/nf_funannotate1 -profile annotate,nrp \
    -params-file /data/src/nf_funannotate1/k8s/params_nrp_test.yaml
```

Watch pods from your workstation with
`kubectl get pods -n ucr-stajichlab -w`. Failed task pods are kept, so
`kubectl logs <pod>` shows them. Successful ones are cleaned up.

## Turning features back on

In roughly this order, one per run (`-resume` keeps what already passed):

GeneMark is already on in the smoke test. GENEMARK_RUN runs in
`teambraker/braker3`, which needs no license key (`--gm_key` exists only for
a key-gated `--container_genemark` image).

1. `--run_sra_fetch true`: RNA-seq fetch plus Trinity/PASA training. PASA's
   MariaDB starts inside the funannotate image (`mariadbd` is bundled). If that
   misbehaves, `--pasa_mysql false` falls back to SQLite.
2. `--run_annotate true`: needs `/data/refdb/eggnog_db`, mounted into the pod
   at `/opt/databases/eggnog_db`.

## Not available on k8s yet

- **SignalP, DeepTMHMM, antiSMASH**: the images are locally built `.sif` files
  (SignalP and DeepTMHMM are licensed). They need pushing to a registry the
  pods can pull from. For the private ones, that means NRP's GitLab registry
  plus `pod imagePullSecret: 'gitlab-registry-cred'` (that secret already
  exists in the namespace). GPU runs would use the `accelerator` directive in
  place of the SLURM `--gres` options.
- **FCS-GX purge** (GENOME_CLEAN without `--skip_fcs`): needs the ~470 GB
  database in `/dev/shm`.
- **prodigal** (`--run_prodigal`): host mode needs prodigal and python in one
  image. The funannotate image is the default and has not been checked for
  prodigal.

## Known risks / untested

- Pods get the namespace default of 50 Gi ephemeral storage for `/tmp`.
  Most tools write into the task work dir on the PVC, but a tool that puts
  large temp files in `/tmp` would be evicted.
- `k8s_max_cpus` / `k8s_max_memory` (32 / 128 GB) cap every task. Raise them
  if NRP schedules bigger pods for you in reasonable time.
- The pod annotation `nrp-nautilus.io/project` reuses the BFD project value.
