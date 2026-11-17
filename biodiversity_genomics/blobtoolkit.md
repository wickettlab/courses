# Exercise: Screening the Assembly for Contamination with BlobToolKit

Biodiversity Genomics

Norman Wickett and Daria Shipilina

17 November, 2026

## 1. Introduction

BlobToolKit checks an assembly for contamination. Every scaffold is placed by its GC content and read coverage, and coloured by its best taxonomic hit (BLAST/DIAMOND). BUSCO completeness is calculated at the same time. We run it with the Nextflow pipeline ```sanger-tol/blobtoolkit```.

Documentation: https://pipelines.tol.sanger.ac.uk/blobtoolkit



## 2. Working Directory

```text
cd /lisc/data/scratch/course/2026w300485/username
mkdir BlobTK
cd BlobTK
```

## 3. Submission Script

```bash
#!/bin/bash
#SBATCH --job-name=blobtk
#SBATCH --output=%j.%x.out
#SBATCH --error=%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --time=16:00:00
#SBATCH --no-requeue
#SBATCH --mail-user=your_email@univie.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL

set -euxo pipefail

module purge
module load Nextflow/26.04.4

### CHANGE THIS ###
### Delete your working folder, create a new one (can be the same name)
BASE="/lisc/data/scratch/course/2026w300485/dasha/BlobTK"  #create a new folder where analysis will run
FASTA="/lisc/data/scratch/course/2026w300485/dasha/Dserotinus_scaffolds_final.fa"  #assembly to evaluate
TAXON="Dianthus serotinus"
TAXON_ID="Dserotinus" # no spaces!
READS="/path/to/Dserotinus_hifi.fastq.gz"   # path to your raw HiFi reads (not Hi-C)

### !COPY  BELOW AGAIN! ###
SAMPLESHEET="${BASE}/samplesheet_blobtools.csv"
OUTDIR="${BASE}/results"
WORKDIR="${BASE}/work"
TMPDIR="${BASE}/tmp"

mkdir -p "${OUTDIR}" "${WORKDIR}" "${TMPDIR}"

# Databases (shared databases from the lab)
#BUSCO_DB="/lisc/data/work/botany/diospyrosProj/Amin_tmp/blob_tool_kit/busco_db/v5/data/"
BUSCO_DB="/lisc/data/scratch/course/2026w300485/databases/data"
#TAXDUMP="/lisc/data/work/botany/diospyrosProj/Amin_tmp/blob_tool_kit/taxdump_2026_04"
TAXDUMP="/lisc/data/scratch/course/2026w300485/databases/taxdump_2026_04"
#BLASTP="/lisc/data/work/botany/diospyrosProj/Amin_tmp/blob_tool_kit/blastp_blastx_db/reference_proteomes.dmnd"
BLASTP="/lisc/data/scratch/course/2026w300485/databases/reference_proteomes.dmnd"
#BLASTX="/lisc/data/work/botany/diospyrosProj/Amin_tmp/blob_tool_kit/blastp_blastx_db/reference_proteomes.dmnd"
BLASTX="/lisc/data/scratch/course/2026w300485/databases/reference_proteomes.dmnd"
BLASTN="/lisc/opt/mirror/ncbi/2026-10-25/nt.nal"
BUSCO_LINEAGES="eukaryota_odb12,embryophyta_odb12,viridiplantae_odb12"

###############################################################################
# SAMPLESHEET
###############################################################################
cat > "${SAMPLESHEET}" <<EOF
sample,datatype,datafile,library_layout
${TAXON_ID},pacbio,${READS},SINGLE
EOF

###############################################################################
# PARAMS FILE
###############################################################################
cat > "${BASE}/params.yaml" <<EOF
align: true
use_work_dir_as_temp: true
EOF

###############################################################################
# CUSTOM PROCESS CONFIG
###############################################################################
cat > "${BASE}/custom.config" <<EOF
process {
    withName: 'DIAMOND_BLASTX' {
        cpus   = { 24 * task.attempt }
        time   = { 48.h * task.attempt }
    }
}
EOF

###############################################################################
# NEXTFLOW ENV
###############################################################################
export NXF_HOME="${BASE}/.nextflow_home"
export NXF_WORK="${WORKDIR}"
export NXF_TEMP="${TMPDIR}"
export TMPDIR="${TMPDIR}"

export NXF_SINGULARITY_CACHEDIR="${BASE}/singularity_cache"
mkdir -p "${NXF_SINGULARITY_CACHEDIR}"

###############################################################################
# RUN
###############################################################################
nextflow run sanger-tol/blobtoolkit \
    -r 0.10.2 \
    -profile singularity \
    -resume \
    -work-dir "${WORKDIR}" \
    -params-file "${BASE}/params.yaml" \
    -c "${BASE}/custom.config" \
    --input "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    --taxon "${TAXON}" \
    --fasta "${FASTA}" \
    --busco "${BUSCO_DB}" \
    --busco_lineages "${BUSCO_LINEAGES}" \
    --taxdump "${TAXDUMP}" \
    --blastn "${BLASTN}" \
    --blastp "${BLASTP}" \
    --blastx "${BLASTX}" \
    --image_format png \
    --publish_dir_mode copy
```

Before you save, change:

- ```--mail-user``` to your own email
- ```BASE``` to the folder you just created
- ```FASTA``` to the assembly you want to check

## 4. Run
The run takes many hours. If the job stops at the time limit, submit it again: ```-resume``` continues from where it stopped.

## 5. Results

```text
results/blobtoolkit/plots/      static plots (png)
results/blobtoolkit/<assembly>/ BlobDir for the interactive viewer
results/busco/                  BUSCO results
```

- **Snail plot:** assembly size, N50, GC, BUSCO completeness
- **Blob plot:** scaffolds by GC vs. coverage, coloured by taxon
- **Cumulative plot:** assembly span added by each taxon

**Keep a record of everything you have done and save the plots. This will all be included in your protocol.**
