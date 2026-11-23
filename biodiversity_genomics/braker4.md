**Biodiversity Genomics**

# Exercise: Predicting Genes with BRAKER4

Daria Shipilina and Norman Wickett

## 1. Introduction

BRAKER4 predicts protein-coding genes in a genome. It uses two kinds of evidence: RNA-seq reads (where genes are expressed) and proteins from other plants (what genes should look like). GeneMark-ETP makes a first guess from the evidence, AUGUSTUS learns from those guesses and predicts genes again, and TSEBRA merges both into the final gene set. Using RNA-seq and proteins together is implemented in ETP mode.

BRAKER4 is a Snakemake pipeline: a recipe that remembers which steps are already done. Every step runs inside a container (Singularity, now called Apptainer), a sealed kit with the exact tools and versions it needs.

Today each of you sets up your own copy of BRAKER4 and runs it on Slurm with a small test: a 5 Mb piece of the *Dianthus serotinus* genome and RNA-seq from two libraries.

Documentation: https://github.com/Gaius-Augustus/BRAKER4


## 2. Overview

This is where everything will be once you are done with the setup:

```text
/lisc/data/scratch/course/2026w300485/
|-- BRAKER4_SHARED/            prepared for you: test data, containers, 3 run files
`-- username/
    `-- BRAKER/
        |-- BRAKER4/           your copy of the BRAKER4 code
        `-- Dsero_test5Mb/     your run folder: your 3 run files and all results
```

## 3. Make Your BRAKER Folder and Clone BRAKER4

```bash
cd /lisc/data/scratch/course/2026w300485/<your_username>
mkdir BRAKER
cd BRAKER
git clone https://github.com/Gaius-Augustus/BRAKER4.git
ls BRAKER4
```

You should see a list of files, including `Snakefile`. This is the recipe for the whole pipeline.


## 4. The Shared Folder

```bash
ls -lh /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED
```

It holds:

- **Test data:** a 5 Mb piece of `scaffold_1` of the *D. serotinus* genome, and RNA-seq reads from libraries 19L and 19R (the first 4 million read pairs of each).
- **`samples.csv`, `config.ini`, `run_slurm.sh`:** the three files that describe and start a run.
- **`singularity_cache/`:** all BRAKER4 containers, already downloaded.

You copy only the three run files. The test data and containers stay in the shared folder (`samples.csv` already points there), so we don't make 20 copies of the same data. Don't change anything inside `BRAKER4_SHARED`.

## 5. Make Your Run Folder and Copy the Run Files

```bash
cd /lisc/data/scratch/course/2026w300485/<your_username>/BRAKER
mkdir Dsero_test5Mb
cd Dsero_test5Mb
cp /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED/samples.csv .
cp /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED/config.ini .
cp /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED/run_slurm.sh .
ls -l
```

You should see the three files.

## 6. Look Inside the Run Files

### samples.csv: what to annotate

| Column | What it holds |
| --- | --- |
| `sample_name` | `Dsero_test`; also the name of the output folder |
| `genome` | the 5 Mb test genome (not repeat-masked) |
| `protein_fasta` | plant proteins from OrthoDB (`Viridiplantae.fa.gz`) |
| `fastq_r1`, `fastq_r2` | read 1 and read 2 files; libraries are separated by `:` |
| `busco_lineage` | `eukaryota_odb12`, used to measure completeness |

Empty columns are options we don't use today.

### config.ini: how to annotate

```bash
cat config.ini
```

The most important settings:

| Setting | Meaning |
| --- | --- |
| `masking_tool = red` | BRAKER4 masks repeats itself with Red, because our genome is not masked |
| `skip_optimize_augustus = 0` | AUGUSTUS training is optimised: slower, but more accurate |
| `auto_mode_fallback = 0` | if the evidence is too thin for ETP mode, stop with an error instead of quietly switching to a simpler mode |
| `busco_download_path` | where the BUSCO datasets are; shared by everyone |

**Slurm settings.** The `[SLURM_ARGS]` section sets the resources every step asks Slurm for:

```ini
[SLURM_ARGS]
cpus_per_task = 16
mem_of_node = 64000
max_runtime = 240
```

- `cpus_per_task`: CPUs per job
- `mem_of_node`: memory in MB (64000 = 64 GB)
- `max_runtime`: time limit per job in minutes (240 = 4 hours)

Each step of the pipeline becomes its own Slurm job with these limits.

### run_slurm.sh: how to start

```bash
cat run_slurm.sh
```

```bash
#!/bin/bash
# Start from your run folder with: bash run_slurm.sh
COURSE=/lisc/data/scratch/course/2026w300485

# stop if we are not in a run folder
if [ ! -f config.ini ]; then echo "No config.ini here: cd into your run folder first"; exit 1; fi

module load snakemake/8.18.2-foss-2025b
export BRAKER4_CONFIG="$(pwd)/config.ini"
export APPTAINERENV_PYTHONPATH=""

snakemake \
    --snakefile ../BRAKER4/Snakefile \
    --executor slurm \
    --jobs 20 \
    --use-singularity \
    --singularity-prefix "$COURSE/BRAKER4_SHARED/singularity_cache" \
    --singularity-args "-B $COURSE" \
    --envvars APPTAINERENV_PYTHONPATH \
    --latency-wait 120 \
    --printshellcmds
```

| Line | What it does |
| --- | --- |
| `module load snakemake/...` | borrows Snakemake from the cluster's software library |
| `export BRAKER4_CONFIG=...` | tells BRAKER4 to use the `config.ini` in this folder |
| `export APPTAINERENV_PYTHONPATH=""` and `--envvars ...` | keep the cluster's Python settings out of the containers, also inside Slurm jobs |
| `--snakefile ../BRAKER4/Snakefile` | the recipe, from your own BRAKER4 copy one folder up |
| `--executor slurm` | sends every step to Slurm as a separate job |
| `--jobs 20` | at most 20 jobs at the same time |
| `--use-singularity` | runs each step inside its container |
| `--singularity-prefix` | where the containers are; shared, so nobody downloads them again |
| `--singularity-args "-B ..."` | the door that lets the containers see the course folder. Every folder with input files must be listed here, or steps fail at once |
| `--latency-wait 120` | waits up to 2 minutes for files written on compute nodes to appear |
| `--printshellcmds` | prints the command of each step |

## 7. Open a screen Session

The run takes about an hour. Snakemake must keep running the whole time, even after you close your laptop. screen keeps your terminal session alive on login01 when you disconnect, like leaving the lab with the experiment still running.

```bash
screen -S braker4
cd /lisc/data/scratch/course/2026w300485/username/BRAKER/Dsero_test5Mb        # <- CHANGE username
```

screen keys (on a Mac, use Control, not Command):

Leave the session running: press Control+A, release both, then press D.
Go back in: `screen -r braker4`. If screen says the session is still attached (for example after your connection dropped), use `screen -d -r braker4`.
List your sessions: `screen -ls`
Scroll up to read earlier output: Control+A, release, then Esc; move with the arrow keys; press Esc again to stop scrolling.
End the session once the run is finished: type exit inside it.

## 8. Load Snakemake and Set the Environment

```bash
module load snakemake/8.18.2-foss-2025b
export BRAKER4_CONFIG="$(pwd)/config.ini"
export APPTAINERENV_PYTHONPATH=""
echo $BRAKER4_CONFIG
```

The last line should print a path ending in `Dsero_test5Mb/config.ini`.

These settings only exist in this terminal window. In a new window, or after logging in again, run them again. (`run_slurm.sh` sets them for itself.)

Type the quotes in `""` as plain quotes with nothing in between. Curly quotes copied from a web page (`“”`) break this silently.

## 9. Dry Run

```bash
snakemake --snakefile ../BRAKER4/Snakefile -n
```

`-n` means dry run: Snakemake reads the recipe, checks that every input file exists and lists all steps, but runs nothing. The output should end with `This was a dry-run`. Above it, the `Job stats` table lists every step (rule) and how many times it will run. Which steps do you recognise from the lecture?

If you see `MissingInputException`, a file is not where Snakemake expects it. Check that you are in your run folder, and look at the paths in `samples.csv`.

## 10. Start the Run on Slurm

```bash
bash run_slurm.sh
```

Snakemake stays in your tmux window on the login node and sends each step to Slurm as its own job, which runs on a compute node. You can now leave tmux (Control+B, release, then D); the run keeps going.

## 11. Watch the Run

Go back in with `tmux attach -t braker4`. After each finished step Snakemake prints a line like `12 of 40 steps (30%) done`; the run is finished when both numbers match and it says `(100%) done`. To see your jobs in the Slurm queue, run `squeue -u $USER` (`$USER` is filled in with your username; `R` = running, `PD` = waiting). Expect about an hour, plus waiting time in the queue.

If a step fails, read the lines just above `Error in rule`, or open the log file it names (in `logs/Dsero_test/`). After fixing the problem, run `bash run_slurm.sh` again: finished steps are not repeated. If Snakemake says the folder is locked, first run `snakemake --snakefile ../BRAKER4/Snakefile --unlock`.

## 12. Results

The final results are in `output/Dsero_test/results/`:

```text
braker.gff3.gz, braker.gtf.gz   predicted genes (the annotation)
braker.aa.gz                    predicted protein sequences
braker.codingseq.gz             predicted coding sequences
braker_report.html              summary report with plots
quality_control/                BUSCO, gene statistics, training summary, runtime plot
```

Open `braker_report.html` first: right-click it in the VS Code Explorer, choose Download, and open it in your browser. To count the predicted genes:

```bash
zcat output/Dsero_test/results/braker.gff3.gz | awk '$3 == "gene"' | wc -l
```

The test genome is 5 Mb of an 848 Mb genome (less than 1%), so a small gene set and very low BUSCO completeness are expected. How long each step took is shown in `quality_control/runtime_plot.png`.

# Part 2: Annotating Your Own Genome

## 13. What Is Different from the Test

You run BRAKER4 exactly as in Part 1, with the same `run_slurm.sh`. Four things change:

- **Genome:** your own genome instead of the 5 Mb test piece.
- **Masking:** you soft-mask your genome yourself with your EDTA results, instead of letting BRAKER4 mask it with Red.
- **RNA-seq:** all RNA-seq libraries of your species instead of the two test libraries.
- **Time:** the whole run takes 2–4 days, so each job is allowed up to 3 days instead of 4 hours.

Everyone runs their own annotation. Start as early as you can: if something goes wrong on day 2, you still need time to rerun.

You will make a new run folder, copy `samples.csv` and `config.ini` (and `run_slurm.sh`) from the shared folder again, soft-mask your genome, and edit `samples.csv` and `config.ini` in VS Code. Sections 14–19 go through every step.

## 14. Make Your Run Folder and Copy the Run Files

Name the folder after your species, for example `Dserotinus_full`, `Dsuperbus_full` or `Pnigrum_full`. It **must** sit inside your `BRAKER` folder, next to `BRAKER4`, because `run_slurm.sh` looks for the code one folder up (`../BRAKER4/Snakefile`).

```bash
cd /lisc/data/scratch/course/2026w300485/username/BRAKER        # <- CHANGE username
mkdir Dserotinus_full                                           # <- CHANGE to your species
cd Dserotinus_full                                              # <- CHANGE to your species
cp /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED/samples.csv .
cp /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED/config.ini .
cp /lisc/data/scratch/course/2026w300485/BRAKER4_SHARED/run_slurm.sh .
ls -l
```

You should see the three files. `run_slurm.sh` needs no changes.

## 15. Your Genome

> **[? NORM]** Students choose their own genome. Is this right? Which file should they use (for example the scaffolded assembly after Hi-C, or the hifiasm primary assembly)? ?

> **[? NORM]** Should students remove anything from the genome before annotation, for example scaffolds shorter than some length, or organelle/contaminant scaffolds flagged by BlobToolKit? If yes, which cutoff, and we add the step here. ?

## 16. Soft-mask Your Genome with Your EDTA Results

Repeats should be masked before gene prediction, so that transposable elements are not mistaken for genes. BRAKER4 needs a **soft-masked** genome: repeats in lowercase letters (`acgt`), everything else in uppercase (`ACGT`). This works better than hard-masking (repeats replaced by `N`), because the gene finders can still read the sequence. You already annotated the transposable elements of your genome with EDTA; here you turn that annotation into a soft-masked genome. It takes a few minutes.


**1. Find your EDTA annotation.** EDTA writes the whole-genome TE annotation to a file ending in `.EDTA.TEanno.gff3`. Search your folder for it:

```bash
find /lisc/data/scratch/course/2026w300485/username -name "*.EDTA.TEanno.gff3"        # <- CHANGE username
```

If several files are listed, take the one from the EDTA run on your whole genome (the same genome file you annotate now).

**3. Save both paths in variables**, so you don't have to type them again:

```bash
GENOME=MY_GENOME_PATH        # <- CHANGE: full path of your genome (section 15)
EDTA=MY_EDTA_GFF3            # <- CHANGE: full path of the .EDTA.TEanno.gff3 file from step 2
ls -lh $GENOME $EDTA
```

`ls` must list both files. The variables only exist in this terminal; if you start again later, set them again.

**4. Check that the sequence names match.** The masking only works if EDTA used exactly the same sequence names as your genome:

```bash
grep ">" $GENOME | cut -d' ' -f1 | tr -d '>' | sort -u > genome_names.txt
awk -F'\t' '!/^#/ && NF>=5 {print $1}' $EDTA | sort -u > edta_names.txt
wc -l genome_names.txt edta_names.txt
comm -13 genome_names.txt edta_names.txt | wc -l
```

- `wc -l` shows the number of sequences in your genome and the number of sequences with at least one TE. The second number can be smaller (short scaffolds may have no TEs), but not larger.
- The last command counts names that are in the EDTA file but not in your genome. It must print `0`. If not, EDTA was run on another genome file or changed the names: stop here, go to "If you can't use EDTA" at the end of this section, and tell us.

**5. Convert the annotation to BED format.**

```bash
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=5 {print $1, $4-1, $5}' $EDTA > edta_repeats.bed
head -3 edta_repeats.bed
wc -l edta_repeats.bed
```

BED has three columns: sequence name, start, end. GFF3 counts positions from 1 and BED from 0, which is why the start is `$4-1`. `wc -l` shows how many TE pieces will be masked.

**6. Soft-mask the genome.**

```bash
module load BEDTools
bedtools maskfasta -fi $GENOME -bed edta_repeats.bed -fo genome_softmasked.fa -soft
```

`-soft` turns the repeat regions into lowercase instead of `N`. If `module load BEDTools` fails, find the exact module name with `module avail bedtools` (upper and lower case matter).

**7. Check the result.**

```bash
grep -c ">" genome_softmasked.fa
grep -v ">" genome_softmasked.fa | awk '{t += length($0); m += gsub(/[acgtn]/, "")} END {printf "%.1f%% of the genome is soft-masked\n", 100*m/t}'
```

- The number of sequences must be the same as in your genome (section 15).
- The masked percentage should be close to the total TE percentage in EDTA's summary file (ending in `.EDTA.TEanno.sum`, in the same folder as the GFF3). If it is `0.0%`, the masking did not work.

**8. Get the path**

```bash
realpath genome_softmasked.fa
exit
```

Keep the path

## 17. Your RNA-seq Reads

All RNA-seq reads are in `/lisc/data/scratch/course/2026w300485/RNAseq/`, one folder per library. You read them straight from there; **don't copy them** (they are tens of GB). See what is there:

```bash
ls /lisc/data/scratch/course/2026w300485/RNAseq/
```

Use only the RNA-seq of **your** species:

| Species | Folders in `RNAseq/` | Libraries |
| --- | --- | --- |
| *Dianthus serotinus* | `Dsero_19L`, `Dsero_19R`, `Dsero_2AB`, `Dsero_2BB` | 4 |
| *Dianthus superbus* | the three `Dsuper_...` folders | 3 |
| *Piper nigrum* | `Pnigrum_PRJEB38192` (all samples in one folder) | one per pair of files |

Every library has two files: read 1 (ends in `_1.fq.gz` or `_1.fastq.gz`) and read 2 (ends in `_2.fq.gz` or `_2.fastq.gz`). BRAKER4 needs all read-1 files in one list and all read-2 files in a second list, each separated by `:`. These commands make the two lists for you and save them as `reads_1.txt` and `reads_2.txt` in your run folder. Run **only the block for your species**, from inside your run folder.

*Dianthus serotinus:*

```bash
ls /lisc/data/scratch/course/2026w300485/RNAseq/Dsero_*/*_1.f*q.gz | paste -sd: | tr -d '\n' > reads_1.txt
ls /lisc/data/scratch/course/2026w300485/RNAseq/Dsero_*/*_2.f*q.gz | paste -sd: | tr -d '\n' > reads_2.txt
```

*Dianthus superbus:*

```bash
ls /lisc/data/scratch/course/2026w300485/RNAseq/Dsuper_*/*_1.f*q.gz | paste -sd: | tr -d '\n' > reads_1.txt
ls /lisc/data/scratch/course/2026w300485/RNAseq/Dsuper_*/*_2.f*q.gz | paste -sd: | tr -d '\n' > reads_2.txt
```

*Piper nigrum:*

```bash
ls /lisc/data/scratch/course/2026w300485/RNAseq/Pnigrum_PRJEB38192/*_1.f*q.gz | paste -sd: | tr -d '\n' > reads_1.txt
ls /lisc/data/scratch/course/2026w300485/RNAseq/Pnigrum_PRJEB38192/*_2.f*q.gz | paste -sd: | tr -d '\n' > reads_2.txt
```

What the commands do: `ls` lists every read-1 (or read-2) file with its full path (`*` matches any text, and `f*q.gz` matches both `fq.gz` and `fastq.gz`); `paste -sd:` joins the list into one line separated by `:`; `tr -d '\n'` removes the line break at the end, so the text pastes cleanly.

**Check the lists:**

```bash
tr ':' '\n' < reads_1.txt; echo
tr ':' '\n' < reads_2.txt; echo
```

Each command prints the files one per line. Both lists must have the same number of files, in the same order (the same library on the same line, only `_1` vs `_2` differs). `ls` sorts by name, so this is the case automatically.

## 18. Edit samples.csv in VS Code

The copied `samples.csv` still describes the test. You replace its content with a template and fill in five values.

**1. Open the file.** In the VS Code Explorer, go to `BRAKER` > your run folder and click `samples.csv`. Turn on word wrap so you can see long lines: **View > Word Wrap** (or Alt+Z; on a Mac, Option+Z).

**2. Replace the content with the template.** Select everything (Ctrl+A; Mac: Cmd+A), delete it, and paste these two lines:

```text
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
MY_SAMPLE_NAME,MY_GENOME_PATH,MY_MASKED_GENOME_PATH,/lisc/data/scratch/course/2026w300485/databases/Viridiplantae.fa.gz,,MY_READS_1,MY_READS_2,,,,,,eukaryota_odb12,
```

Line 1 is the header: never change it. Line 2 is your genome. Proteins (`Viridiplantae.fa.gz`) and the BUSCO lineage (`eukaryota_odb12`) stay as in the test for every species.

**3. Fill in the five placeholders.** For each one: double-click the placeholder in line 2 (VS Code selects the whole word, including the underscores), then paste your value (Ctrl+V; Mac: Cmd+V).

| Placeholder | Replace with | Where to get it |
| --- | --- | --- |
| `MY_SAMPLE_NAME` | a short name for your genome, e.g. `Dserotinus`; letters, numbers and underscores only. It becomes the output folder `output/MY_SAMPLE_NAME/` | type it |
| `MY_GENOME_PATH` | the full path of your genome (the plain, unmasked one) | section 15 (Copy Path) |
| `MY_MASKED_GENOME_PATH` | the full path of `genome_softmasked.fa`. If you skipped section 16, delete the placeholder so that two commas stay next to each other (`,,`); BRAKER4 then masks with Red | section 16, step 8 |
| `MY_READS_1` | the read-1 list | open `reads_1.txt` in VS Code, Ctrl+A, Ctrl+C (Mac: Cmd+A, Cmd+C) |
| `MY_READS_2` | the read-2 list | open `reads_2.txt` in VS Code, select all and copy the same way |

Rules for line 2:

- Don't add spaces or quotes anywhere.
- Don't delete or add commas. The empty places between two commas (`,,`) are options we don't use.
- The whole of line 2 must stay **one** line. With word wrap on it may look like several lines; the line number on the left shows only `1` and `2`.

**4. Save** with Ctrl+S (Mac: Cmd+S). An unsaved file shows a dot next to its name in the tab.

This is what a finished line 2 looks like (Daria's *D. serotinus* genome with all four libraries). It has no masked genome, so its third field is empty (`,,` right after the genome path); in your line, the path of `genome_softmasked.fa` goes there:

```text
Dserotinus,/lisc/data/scratch/course/2026w300485/dasha/Dserotinus_scaffolds_final.fa,,/lisc/data/scratch/course/2026w300485/databases/Viridiplantae.fa.gz,,/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_19L/Dsero_19L_MKRN260064944-1A_23M3G5LT4_L1_1.fq.gz:/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_19R/Dsero_19R_MKRN260064947-1A_23M3G5LT4_L1_1.fq.gz:/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_2AB/Dsero_2AB_MKRN260064940-1A_23M3G5LT4_L1_1.fq.gz:/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_2BB/Dsero_2BB_MKRN260064941-1A_23M3G5LT4_L1_1.fq.gz,/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_19L/Dsero_19L_MKRN260064944-1A_23M3G5LT4_L1_2.fq.gz:/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_19R/Dsero_19R_MKRN260064947-1A_23M3G5LT4_L1_2.fq.gz:/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_2AB/Dsero_2AB_MKRN260064940-1A_23M3G5LT4_L1_2.fq.gz:/lisc/data/scratch/course/2026w300485/RNAseq/Dsero_2BB/Dsero_2BB_MKRN260064941-1A_23M3G5LT4_L1_2.fq.gz,,,,,,eukaryota_odb12,
```

**5. Check the file** in the terminal, from your run folder:

```bash
awk -F',' '{print NF}' samples.csv
```

This counts the columns in each line. It must print `14` twice. Any other number means a comma was deleted or added.

```bash
paste -d' ' <(head -1 samples.csv | tr ',' '\n') <(tail -1 samples.csv | tr ',' '\n')
```

Each column next to its value, as in section 6. Check that `sample_name`, `genome`, `genome_masked`, `protein_fasta`, `fastq_r1`, `fastq_r2` and `busco_lineage` are filled (`genome_masked` is empty only if you skipped section 16).

```bash
awk -F',' 'NR==2 {print split($6,a,":"), split($7,b,":")}' samples.csv
```

The number of read-1 and read-2 files. Both numbers must be the same (4 for *D. serotinus*).

```bash
tail -1 samples.csv | tr ',:' '\n\n' | grep '^/' | xargs ls -lh
```

Lists every file named in line 2 with its size. If you see `No such file or directory`, that path has a typo: fix it in VS Code, save, and check again.

## 19. Edit config.ini in VS Code

The test allowed each job 4 hours. On a full genome, single steps (GeneMark-ETP, AUGUSTUS) take much longer, so you raise the limit to 3 days.

1. Click `config.ini` in your run folder in the VS Code Explorer.
2. Find the line `max_runtime = 240` in the `[SLURM_ARGS]` section (Ctrl+F; Mac: Cmd+F, and type `max_runtime`).
3. Change `240` to `4320` (minutes = 3 days), so the line reads:

```ini
max_runtime = 4320
```

4. Save (Ctrl+S; Mac: Cmd+S).

Leave everything else as it is: 16 CPUs and 64 GB per job, `masking_tool = red` (only used if `genome_masked` is empty), AUGUSTUS optimisation on, and `auto_mode_fallback = 0`, so the run stops with an error instead of quietly leaving ETP mode if your evidence is too thin.

Check in the terminal:

```bash
grep max_runtime config.ini
```

It must print `max_runtime = 4320`.

> **[? OPEN – Daria/Norm]** Switch off AUGUSTUS optimisation for the student runs (`skip_optimize_augustus = 1`)? It saves several hours per run and costs about 2–3% accuracy.

## 20. Dry Run and Start

Open a **new** screen session with a different name, so it doesn't mix with the test:

```bash
screen -S braker_full
cd /lisc/data/scratch/course/2026w300485/username/BRAKER/Dserotinus_full        # <- CHANGE username and folder
module load snakemake/8.18.2-foss-2025b
export BRAKER4_CONFIG="$(pwd)/config.ini"
export APPTAINERENV_PYTHONPATH=""
echo $BRAKER4_CONFIG
```

`echo` must print the path to the `config.ini` in your new run folder (not the test folder).

Dry run:

```bash
snakemake --snakefile ../BRAKER4/Snakefile -n
```

It must end with `This was a dry-run`. In the `Job stats` table, `hisat2_align` should now run once per RNA-seq library (4 times for *D. serotinus*), instead of twice as in the test. If you see `MissingInputException`, go back to the checks in section 18.

Start the run:

```bash
bash run_slurm.sh
```

Leave screen with Control+A, release, then D. The run keeps going.

## 21. Watch the Full Run

Check on your run once or twice a day:

```bash
screen -r braker_full
```

Read the last `... of ... steps (...%) done` line, then leave again with Control+A, D. Outside screen, `squeue -u $USER` shows your running (`R`) and waiting (`PD`) jobs.

- **The run takes 2–4 days.** GeneMark-ETP, the AUGUSTUS optimisation and the AUGUSTUS prediction take most of it. The progress line can stay on the same step for a day; that is normal as long as the job is listed in `squeue`.
- **If a step fails:** scroll up (Control+A, then Esc) and read the lines just above `Error in rule`, or open the log file it names in `logs/MY_SAMPLE_NAME/`. If you can't solve it, send us those lines.
- **If your screen session is gone** (`screen -ls` shows nothing, for example after the login node was restarted), the run is not lost. Start a new session and continue; finished steps are not repeated:

```bash
screen -S braker_full
cd /lisc/data/scratch/course/2026w300485/username/BRAKER/Dserotinus_full        # <- CHANGE username and folder
bash run_slurm.sh
```

- **If Snakemake says the folder is locked:** unlock it first, then start again:

```bash
module load snakemake/8.18.2-foss-2025b
export BRAKER4_CONFIG="$(pwd)/config.ini"
export APPTAINERENV_PYTHONPATH=""
snakemake --snakefile ../BRAKER4/Snakefile --unlock
bash run_slurm.sh
```

- Don't delete or move anything in your run folder while the run is going.

## 22. Results of the Full Run

The results are in `output/MY_SAMPLE_NAME/results/`, with the same files as in the test (section 12). Set your sample name once, then run the checks:

```bash
SAMPLE=Dserotinus                                   # <- CHANGE to your sample_name
zcat output/$SAMPLE/results/braker.gff3.gz | awk '$3 == "gene"' | wc -l
cat output/$SAMPLE/results/quality_control/busco_summary.txt
cat output/$SAMPLE/results/quality_control/compleasm_summary.txt
cat output/$SAMPLE/results/quality_control/gene_set_statistics.txt
```

- The first command counts the predicted genes.
- `busco_summary.txt` and `compleasm_summary.txt` show how complete the gene set is (how many of the expected single-copy genes were found).
- `gene_set_statistics.txt` gives numbers of genes, transcripts, exons and their lengths.

Download `braker_report.html` (right-click in the VS Code Explorer > Download) and open it in your browser. Compare it with your test report: number of genes, BUSCO completeness, the training summary and the runtime plot.

For your protocol, record:

- BRAKER4 version (commit ID and date from section 3)
- your genome file, its number of sequences and its size
- the RNA-seq libraries you used, the protein set (`Viridiplantae.fa.gz`) and the BUSCO lineage (`eukaryota_odb12`)
- the settings: ETP mode, masking (EDTA soft-masking and the percentage masked, or Red), 16 CPUs, 64 GB, any settings you changed
- how long the run took
- number of genes and BUSCO/compleasm completeness

**Keep a record of everything you have done and save the reports. This will all be included in your protocol.**
