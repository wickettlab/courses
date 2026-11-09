# Exercise: Hi-C Scaffolding with YaHS and Pretext

Biodiversity Genomics  
Daria Shipilina and Norman Wickett

10 November, 2026

## 1. Introduction

Today we use Hi-C data to turn our contigs into chromosome-scale scaffolds. There will be **two jobs**. Job 1 maps the Hi-C reads to your assembly and will run during the class. When it has finished, you run job 2: scaffolding with YaHS, preparing the pairs file, and making a contact map with Pretext. While the mapping is running, we will learn manual curation with Pretext maps we created for you before.

## 2. Mapping Hi-C Reads

Mapping finds the place in a reference where each read comes from; here the reference is your own assembly. Hi-C read pairs are special: the two reads of a pair come from parts of the genome that were close to each other in the nucleus, so they can land far apart, even on different contigs. This is the information YaHS will use to order and orient the contigs.

### 2.1 Select your best assembly

Look at the statistics from the previous exercises and select your best assembly. **One assembly per group, and one job per group.** Write down which assembly you chose and why in your `run_log.md`.

### 2.2 Hi-C data

All Hi-C reads are in:

```text
/lisc/data/scratch/course/2026w300485/HiC_data/
```

*Dianthus serotinus*:

```text
Dserotinus_DTG-OmniC-1367_MKDL260019622-1A_23NVKWLT4_L8_1.fq.gz
Dserotinus_DTG-OmniC-1367_MKDL260019622-1A_23NVKWLT4_L8_2.fq.gz
```

*Dianthus superbus*:

```text
Dsuperbus_HiC_R1.fq.gz
Dsuperbus_HiC_R2.fq.gz
```

*Piper nigrum*:

```text
PNIG_SRR8820208.lite.1_1.fastq
PNIG_SRR8820208.lite.1_2.fastq
```

### 2.3 Commands

```bash
module load SAMtools BWA

samtools faidx contigs.fa
bwa index contigs.fa

bwa mem -5SP -t 16 contigs.fa hic_R1.fastq.gz hic_R2.fastq.gz \
  | samtools view -b -F 0x904 - \
  | samtools sort -n -@ 16 -o hic.bam -
```

Replace `contigs.fa` with your selected assembly and `hic_R1`/`hic_R2` with the Hi-C files of your species.

**Note:** `-5SP` are the Hi-C settings: they allow the two reads of a pair to land far apart and handle reads that span a ligation junction. `samtools view -F 0x904` drops reads that didn't map and extra (secondary and supplementary) alignments of the same read. `samtools sort -n` puts the two reads of each pair next to each other.

### 2.4 Job 1

Write a submission script with the commands above. Request 16 CPUs and the following memory and time for your species:

| Species | Memory | Time |
|---|---|---|
| *Dianthus serotinus* | 50G | 24 hours |
| *Dianthus superbus* | 50G | 24 hours |
| *Piper nigrum* | 25G | 16 hours |
All 16 CPU

## 3. Scaffolding with YaHS

YaHS (Yet another Hi-C Scaffolding tool) takes your contigs and the Hi-C alignments. It first breaks contigs where Hi-C coverage drops, as these are possible assembly errors. It then cuts the contigs into chunks, counts the Hi-C links between chunks, and uses these counts to decide which contig ends belong together and in which orientation. It does this in several rounds, from small chunks to large ones, and joins neighbouring contigs with a stretch of `N`s.

**Run this only after job 1 has finished!** YaHS needs the complete BAM file.

```bash
module load SAMtools YaHS

yahs -o asm contigs.fa hic.bam

samtools faidx asm_scaffolds_final.fa
```

YaHS writes the scaffolds (`asm_scaffolds_final.fa`), an AGP file describing how contigs were joined into scaffolds (`asm_scaffolds_final.agp`), and the Hi-C links in binary format (`asm.bin`). `samtools faidx` indexes the new scaffolds; we need the scaffold names and lengths for the next step.

## 4. Preparing the Pairs File

**Note:** header: the pairs-format line, one line per scaffold with its length, and the column names.

```bash
gawk 'BEGIN { print "## pairs format v1.0" }
      { print "#chromsize:", $1, $2 }
      END { print "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2" }' \
      asm_scaffolds_final.fa.fai > asm.pairs
```

**Note:** the Hi-C links, moved from contig to scaffold positions with `juicer pre`.

```bash
juicer pre asm.bin asm_scaffolds_final.agp contigs.fa.fai \
  | gawk 'BEGIN { OFS = "\t" }
          $3 > 0 && $7 > 0 { print ".", $2, $3, $6, $7, ".", "." }' \
  | LC_ALL=C sort -k2,2d -k4,4d -S 8G -T . >> asm.pairs
```

Note the double `>>`: it adds the links to the file after the header. A single `>` would overwrite the header.

## 5. PretextMap and PretextSnapshot

```bash
module load Conda
conda activate /lisc/data/scratch/botany/wickett_lab/Dasha/D_superbus_assembly/envs/pretext

PretextMap -o asm.pretext < asm.pairs
PretextSnapshot --sequences "=full" --map asm.pretext --prefix asm. --folder .
```

PretextMap builds the contact map (`asm.pretext`) that we open in PretextView. PretextSnapshot saves a picture of the full map.

### 5.1 Job 2

Sections 3, 4 and 5 go together into job 2. Request 4 CPUs, 50G memory and 12 hours.

## 6. Manual Curation

While everything is mapping, we will use the Pretext files we created for you before to learn how manual curation works. From there, we will move on to your personal laptops and do the manual curation.

Install PretextView on your laptop:

https://github.com/sanger-tol/PretextView/releases

The Pretext files are here:
/lisc/data/scratch/course/2026w300485/HiC_data/

