# Biodiversity Genomes

## Day 4: Assembling Genomes with hifiasm

### Norman Wickett and Daria Shipilina

### 27 October, 2026

Last week you got to know the PacBio HiFi data for *Piper nigrum*,
*Dianthus superbus*, and *Dianthus serotinus*. You used **SeqKit** to
count how many reads and bases we have, **FastK** to count k-mers, and
**GenomeScope** and **Smudgeplot** to estimate genome size,
heterozygosity, repeat content, and ploidy. Today all of those numbers
become useful. We will assemble the reads into contigs with **hifiasm**,
and you will see that almost every decision we make when running an
assembler depends on what we learned about the genome last week.

**Keep your SeqKit, GenomeScope, and Smudgeplot results open.** You will
need them in Section 5.

## 1. How hifiasm works (the short version)

hifiasm (<https://github.com/chhylp123/hifiasm>) is currently the most
widely used assembler for PacBio HiFi reads. It is used by the Darwin
Tree of Life project, the Vertebrate Genomes Project, and ERGA, which is
why we use it in this course. You don't need to know every detail of the
algorithm to use it well, but you do need to know what it is doing in
each step, because each step has parameters that you may have to change.

**Step 1: Count k-mers and find the coverage peaks.** hifiasm starts by
counting k-mers (k = 51 by default) in all of the reads and drawing a
k-mer histogram. This is exactly what FastK and GenomeScope did last
week! From the histogram, hifiasm decides which peak contains the
**heterozygous** k-mers (found on only one of the two copies of a
chromosome) and which contains the **homozygous** k-mers (found on both
copies). The coverage of the homozygous peak is used later as a ruler
for many decisions, so if hifiasm gets this number wrong, many later
steps go wrong too.

**Step 2: Correct sequencing errors without erasing real differences.**
HiFi reads are very accurate, but not perfect. For each read, hifiasm
lines up all of the other reads that overlap it. If only one or two
reads disagree with the rest at a position, that is probably a
sequencing error, and it gets corrected. If several reads support each
of two different bases, that is probably a heterozygous site, and both
versions are kept. This is called **haplotype-aware error correction**,
and hifiasm runs three rounds of it by default.

*Question: Why would simply taking the majority base at each position be
a bad idea for a diploid genome?*

**Step 3: Build a string graph from the overlaps.** After error
correction, hifiasm finds all overlaps between the corrected reads and
builds a string graph, just like the overlap graphs we drew in the
lecture. Where the two copies of a chromosome are identical, the graph
is a single path. Where the two copies differ, the path splits into two
and joins again, forming a **bubble**:

Stretches of the graph without any branches are called **unitigs**.
hifiasm writes out two versions of this graph: the raw unitig graph
(`r_utg`), which keeps everything, and the processed unitig graph
(`p_utg`), where tiny bubbles caused by leftover errors have been
removed.

**Step 4: Walk through the graph to build contigs.** A contig is a path
through the graph. Every time a path reaches a bubble, hifiasm has to
choose one side. Here is the problem: if two bubbles are separated by a
homozygous stretch that is longer than the reads, nothing in the HiFi
data tells hifiasm which side of the first bubble belongs with which
side of the second. With HiFi reads alone, the two haplotypes are
therefore only **partially phased**: over long distances, a "haplotype
1" contig can switch back and forth between the two parental copies.
Hi-C data (or reads from the parents) can solve this, which is why
hifiasm can use them.

**Step 5: Purge duplicated haplotypes.** In regions where the two copies
of a chromosome are very different, they don't form a bubble at all.
Instead they are assembled as two separate contigs, so the same region
ends up in the assembly twice. These extra copies are called
**haplotigs**. hifiasm looks for
contigs that are very similar to another contig *and* have about half of
the homozygous coverage, and removes them from the primary assembly.
This is called **purging**, and it is one of the steps you will be
tuning today.

Steps 1 to 3 take most of the time (many hours for our genomes). Steps 4
and 5 take minutes to an hour. That is why hifiasm saves the corrected
reads and overlaps in three files ending in `.bin`. If these files are
in the directory where you run hifiasm, and you use the same output
prefix, hifiasm skips error correction and goes straight to building the
graph. We will use this trick today so you can try several parameter
sets in one day.

## 2. Finding hifiasm on LiSC

Find the hifiasm module, load it, and open
the help menu. Remember that LiSC is case-sensitive, so check how the
module name is written first:

``` text
module avail hifiasm
```

Once the module is loaded, try opening the help menu with `--help`. Did
it work? hifiasm uses `-h` instead:

``` text
hifiasm -h
```

The first line of the help menu gives you the **usage**:

``` text
Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>
```

This tells you that you first give hifiasm its options, then one or more
read files. Below the usage you'll see a very long list of options,
grouped into categories like "Overlap/Error correction", "Assembly",
"Purge-dups", and "Hi-C-partition". At the very bottom, hifiasm suggests
reading the full documentation at
<https://hifiasm.readthedocs.io/en/latest/>. Keep that page open in your
browser for the rest of the day.

## 3. What hifiasm needs as input

hifiasm needs very little to get started:

- **HiFi reads** in FASTA or FASTQ format. The files can be gzipped
  (`.gz`), and you can list several files one after the other. hifiasm
  ignores the quality scores in FASTQ files, so FASTA works just as
  well.

Our HiFi data are in the same place as last week:

``` text
ls -lh /lisc/data/scratch/course/2026w300485/course_data/HiFi_data
```


The simplest possible hifiasm command looks like this:

``` text
hifiasm -o <prefix> -t <number_of_cpus> <reads.fastq.gz>
```

`-o` sets the **prefix** that all output files will start with, and `-t`
sets the number of CPUs (this should match `--cpus-per-task` in your
submission script). Everything else in this document is about one
question: are hifiasm's defaults right for *our* genomes?

## 4. The parameters that matter

hifiasm has more than 50 options, and most of the time you will not need
to change any of them. But the defaults were designed with one kind of
genome in mind: diploid, moderately heterozygous, with a clean k-mer
histogram (think of a human). The further your genome is from that, the
more likely it is that some defaults need changing. Here are the
parameters that matter most, grouped by the question they answer.

### 4.1 Where is the homozygous peak? (`--hom-cov`)

In Step 1, hifiasm has to decide which peak in its k-mer histogram is
the homozygous peak. Its strategy is simple: it finds the tallest peak,
and then looks for a smaller neighbouring peak. If there is a smaller
peak to the *right*, the tallest peak is the heterozygous one and the
right peak is the homozygous one. If there is a smaller peak to the
*left*, the tallest peak is the homozygous one and the left peak is the
heterozygous one. You can follow this decision in the log file, which
contains lines like these:

``` text
[M::ha_analyze_count] highest: count[30] = ...
[M::ha_analyze_count] left: count[15] = ...
[M::ha_analyze_count] right: none
[M::ha_ft_gen] peak_hom: 30; peak_het: 15
```

You can predict where the peaks *should* be before you even run hifiasm.
GenomeScope reports `kcov`, the k-mer coverage of a single copy of the
genome. For long HiFi reads, k-mer coverage and read coverage are almost
the same, so for a **diploid** genome:

``` text
heterozygous peak (peak_het) ≈ kcov
homozygous peak   (peak_hom) ≈ 2 × kcov
```

You can check this using your SeqKit results. The total number of bases
divided by the haploid genome size gives the coverage of the homozygous
peak. For example, imagine a diploid plant with a 500 Mb genome and 30
Gb of HiFi reads: 30,000 Mb ÷ 500 Mb = 60, so we would expect `peak_hom`
around 60 and `peak_het` around 30.

Later in the log, hifiasm reports the homozygous coverage it actually
used:

``` text
[M::purge_dups] homozygous read coverage threshold: <X>
```

If this number is far from what you expected, the assembly will usually
be too large or too small. You can fix this in two ways:

- `--hom-cov <INT>` tells hifiasm the homozygous coverage directly.
- `--hg-size <size>` gives hifiasm your estimated haploid genome size
  (for example `--hg-size 760m`), so it can work out roughly where the
  homozygous peak should be and pick the closest real peak.

### 4.2 How hard should we purge? (`-l`, `-s`, `--purge-max`)

These three options control Step 5.

`-l` sets the **purge level**: 0 means no purging, 1 means only purge
haplotigs that are completely contained within another contig, 2 means
purge all types of haplotigs, and 3 (the default for HiFi-only
assemblies) means purge all types of haplotigs in the most aggressive
way. The hifiasm documentation recommends `-l0` for inbred or highly
homozygous genomes. Why do you think purging could do harm in a genome
that has almost no heterozygosity?

`-s` sets the **similarity threshold**: how similar two contigs must be
before one is treated as a duplicate of the other. The default is 0.55
with `-l3` and 0.75 with `-l1` or `-l2`. A *lower* value means that less
similar contigs also count as duplicates, so more gets purged. The
documentation suggests lowering `-s` for highly heterozygous genomes,
and when one of the two partially phased haplotypes is much larger than
the other.

`--purge-max` sets an **upper coverage limit** for purging: contigs with
higher coverage than this are never purged. hifiasm chooses this value
automatically and prints it in the log
(`purge duplication coverage threshold`). A larger value makes the
primary assembly more contiguous, but it can also collapse real repeats
or segmental duplications into one copy. A smaller value is more
cautious.

### 4.3 How many copies of each chromosome? (`--n-hap`)

`--n-hap` tells hifiasm the number of haplotypes, and the default is 2.
This is important to understand: **most of hifiasm is designed for
diploids.** The purging step assumes two haplotypes, and the partially
phased and fully phased (Hi-C) outputs only make sense for diploid
genomes. For polyploids, the hifiasm documentation recommends setting
`--n-hap` to the number of haplotypes, using only the **primary
assembly**, and then running extra rounds of purging with a separate
tool such as **purge_dups**. 

Polyploidy also raises a coverage question. hifiasm usually needs **at
least 13x coverage per haplotype**. For a diploid, coverage per
haplotype is the total coverage divided by 2. What is it for a
tetraploid? A hexaploid?

Note that `--n-hap` is also used during error correction, so it should
be set in the very first run.

### 4.4 Which assemblies do we want? (`--primary`, `--h1`, `--h2`)

hifiasm can produce three kinds of output, depending on the options and
the data you give it:

- **Default (HiFi only):** a primary assembly (`bp.p_ctg`) plus two
  partially phased haplotype assemblies (`bp.hap1.p_ctg` and
  `bp.hap2.p_ctg`).
- **With `--primary`:** a primary assembly (`p_ctg`) plus an alternate
  assembly (`a_ctg`) containing the bits that were purged or not chosen.
  The primary assembly is usually more contiguous, but the alternate
  assembly is very fragmented. Using `-l0` also produces this pair.
- **With Hi-C reads (`--h1 <R1> --h2 <R2>`):** two fully phased
  haplotype assemblies (`hic.hap1.p_ctg` and `hic.hap2.p_ctg`). Each
  contig now comes from only one parental copy. Note that hifiasm does
  *not* scaffold contigs into chromosomes, even with Hi-C data. We'll do
  that later in the course with YaHS.

Heng Li has written a very clear explanation of these different kinds of
assemblies:
<https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies>. Have a
look at the figure and decide which panel matches each of the outputs
above.

### 4.5 Cheat sheet

| Option         | Default           | What it controls                              | Works with saved `.bin` files?  |
|----------------|-------------------|-----------------------------------------------|---------------------------------|
| `-t`           | 1                 | number of CPUs                                | yes                             |
| `--hom-cov`    | auto              | homozygous coverage used as the "ruler"       | yes                             |
| `--hg-size`    | auto              | haploid genome size, to help find the peaks   | yes                             |
| `-l`           | 3                 | purge level (0 = no purging)                  | yes                             |
| `-s`           | 0.55 (with `-l3`) | how similar a duplicate must be to be purged  | yes                             |
| `--purge-max`  | auto              | coverage limit for purging                    | yes                             |
| `--n-hap`      | 2                 | number of haplotypes (ploidy)                 | partly; set it in the first run |
| `--primary`    | off               | primary/alternate output instead of hap1/hap2 | yes                             |
| `--h1`, `--h2` | none              | Hi-C read files for full phasing              | yes                             |
| `-f`           | 37                | memory filter for k-mer counting              | yes                             |
| `-D`, `-N`     | 5.0, 100          | repeat handling during error correction       | no; needs a fresh run (`-i`)    |

## 5. Exercise: Design a parameter set for your species

You will work in the same species groups as last week. Your job is to
decide how *you* would run hifiasm on your species, and to justify every
choice. There is no single correct answer, but every choice must be backed by a number or a biological fact. Please consider:

1.  **Ploidy.** Does your Smudgeplot agree with the chromosome counts?
    If two ploidy levels have been reported for your species, which one
    is our plant? What does this mean for `--n-hap`?
2.  **Coverage.** Do you have at least 13x per haplotype? If not, what
    problems would you expect?
3.  **Peaks.** Where do you expect `peak_het` and `peak_hom` in the
    hifiasm log? If your species is polyploid, sketch the k-mer
    histogram you expect and mark every peak. Which peak should hifiasm
    call homozygous, and which one might it pick by mistake?
4.  **Heterozygosity and purging.** Is your genome heterozygous enough
    that purging is needed? Should you keep the default `-l3`? Is the
    heterozygosity so high that `-s` might need to be lower? *Dianthus
    superbus* is gynodioecious (populations contain female and
    hermaphrodite plants). What does that suggest about how much it
    outcrosses, and about its heterozygosity?
5.  **Outputs.** Which assemblies do you want: partially phased
    hap1/hap2, primary/alternate? Why? Which ones would
    actually make sense for a polyploid?
6.  **Expected size.** How large should your primary assembly be? How
    large should each haplotype be? 

Present your plan to the class

See example below:

``` text
hifiasm -o maize_B73 -t 32 -l0 maize-hifi.fasta
```

| Option      | Value        | Why                                                                                                                                                         |
|-------------|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-l`        | 0            | B73 is an inbred line, so the two copies of each chromosome are almost identical. Purging would only remove real sequence.                                  |
| `--hom-cov` | not set      | 70 Gb of reads ÷ 2.4 Gb genome ≈ 29x, so `peak_hom` should be about 29. Check the log, and only set `--hom-cov` if hifiasm picked something very different. |
| `--n-hap`   | 2 (default)  | Maize is diploid.                                                                                                                                           |



Each group can run up to three assemblies:

1.  **Default:** only `-o` and `-t`.
2.  **Your plan** from Part C.
3.  **One test:** the same as your plan, except for **one** changed
    parameter, so you can see what that parameter does. For example,
    compare `-l0` with `-l3`, or set `--hom-cov` to half of your
    expected `peak_hom`.


## 6. Running hifiasm on LiSC

hifiasm writes a lot of files into the directory it is run from, so
first make a directory for your assembly in your personal folder and move into it:

``` text
mkdir hifiasm_<species>
cd hifiasm_<species>
```

Now create your submission script there. Piper nigrum full assembly needs a lot of
resources, so use 32 CPUs, 128G of memory, and an estimated run time of
one day. For both Dianthus species drop number of cpus to 16 and use 16 hours.
Your submission script should look something like this:

``` text
#!/bin/bash
#
#SBATCH --job-name=hifiasm_<species>
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=basic
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00

module load <hifiasm_module>

hifiasm -o <prefix> -t 32 <your_options> \
/lisc/data/scratch/course/2026w300485/course_data/HiFi_data/<your_species_reads> \
2> <prefix>.hifiasm.log
```

**Check with one of the instructors before submitting.** Then submit it
as usual:

``` text
sbatch <your_script>.sh
```

The assembly will run overnight, and you will get an email when it
starts and when it finishes. If you're curious, you can watch it work by
looking at the end of the log file:

``` text
tail <prefix>.hifiasm.log
```

## Further reading

- hifiasm documentation: <https://hifiasm.readthedocs.io/en/latest/>
- Cheng H. et al. (2021) Haplotype-resolved de novo assembly using
  phased assembly graphs with hifiasm. *Nature Methods* 18: 170–175.
- Heng Li's blog post on concepts in phased assemblies:
  <https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies>
- RCAC hifiasm tutorial (maize B73), including a comparison of purge
  levels:
  <https://docs.rcac.purdue.edu/lifesciences/tutorials/hifi-assembly/>
- Hands-on hifiasm exercise from the Eukaryotic Genome Assembly course:
  <https://eukaryotic-genome-assembly.github.io/handsOn_2/>
- Hu L. et al. (2019) The chromosome-scale reference genome of black
  pepper provides insight into piperine biosynthesis. *Nature
  Communications* 10: 4702.
