# Biodiversity Genomes

## Day 3: Exploring PacBio Data

### Norman Wickett and Daria Shipilina

### 20 October, 2026

Last week you were introduced to the concepts of different DNA sequencing
methods and ended the class submitting your very first job to
LiSC. You created a submission script to run FastQC on Illumina RNA-Seq
data (RNA-Seq data is produced by creating a sequencing library from the
expressed RNA from a certain tissue, otherwise known as a
**transcriptom**. In creating that submission script, you learned how to
load modules on LiSC, look at help menus to understand program options,
and specify the paths to the inputs and ouputs for a certain program. We
also talked about two ways to download the results from LiSC onto your
own “local” computer; this is particularly important for results that
can’t be viewed or interpreted from the command line or VS Code, such as
html files or image files. Those two ways were either by right-clicking
or control-clicking the file name in VS Code and choosing to download
it, or using the **secure copy** command, `scp`. The secure copy command
works just like the regular copy command, where you use the command,
then tell it what you want to copy, then tell it where you want it to be
copied. Because we copying something from the remote server, you must
enter your login credentials, so the command works like this:

``` text
scp <user_id>@login01.lisc.univie.ac.at:<path_to_what_you_want_to_download> \
<local_download_location>
```

Don’t forget that you have to run this command *from* your local
machine, not from a terminal that is connected to LiSC!

## 1. PacBio HiFi reads

At this point in the course, it’s time to look at some long-read
sequence data. Today, we’ll be starting our genome assembly adventures
by looking at the PacBio HiFi data for *Piper nigrum*, *Dianthus
superbus*, and *Dianthus serotinus*. The data can be found here:

``` text
/lisc/data/scratch/course/2026w300485/course_data/HiFi_data
```

One of the first things we might want to do is get some basic stats on
our sequencing data. We can do that using **SeqKit**. Go ahead and look
for the SeqKit module, load it, and view the help menu. SeqKit works by
calling the program `seqkit` and then giving it the specific command for
what you want it to do. There are many categories of commands, including
**Basic Operation**, **Format Conversion**, **Searching**, and so on.
Today we will use the `stats` command from the Basic Operation category.

For the exercises today you will be assigned one of the three species
and you will work in pairs (or a group of three if we are an odd
number). Working together, write a submission script for running SeqKit
stats on your data and view the output files. To do this, we can write a
simple submission script. Notice the command for running at as the last
line. **Remember that you must change this line so that it points to the
species you are working with.** Additionally, be careful using the `-o`
flag. In the submission script below, it assumes you are running this
from within the directory where you want the output to be written. This
will probably be your personal directory without our course directory.
You could of course put the full path to the location of the output
files. Your submission script should look something like this:

``` text
#!/bin/bash

#SBATCH --job-name=seqkit_stats
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL

module load SeqKit

seqkit stats -o Pnigrum.hifi-reads.stats \
/lisc/data/scratch/course/2026w300485/course_data/HiFi_data/PNIG01.hifi.fastq.gz
```

Look at the output files with your group. What do they tell us about the
data that we have received from the sequencing center? We will discuss
the differences between the three species in class.

## 2. FastK: A k-mer counter for HQ assembly data sets

Earlier today you will have already learned about k-mers and why they
are so important to the process of genome assembly, as well as the
process of estimating genome size, repeat content, heterozygosity, and
so on. If you want to learn even more about FastK, you can check out its
GitHub site here:
<a href="https://github.com/thegenemyers/FASTK">https://github.com/thegenemyers/FASTK</a>.
But we can also load the module and check the help menu. **Go ahead and
see if you can find the correct module name, load it, and open the help
menu.**

Unfortunetly for us, FastK has not adopted the `-h` option for the help
menu. When that fails, sometimes just typing the command works, so try
typing `FastK` *without* any options and see what happens.

The options we will be using with FastK are:

``` text
-v: Verbose mode, provides lots of output statistics
-T: Number of threads, or CPUs per task (should match SBATCH option)
-N: Path to output directory and output prefix
-P: Path to directory for working/temp files
-M: Memory (should match SBATCH option)
-k: k-mer size (will discuss in class)
-t: tells FastK to save the k-mer table itself, not just the histogram of k-mer counts. You will need this table for downstream analyses, please add -t4 flag to your script
```

For this job, we will use 8 CPUs per task, 16G of memory, and an
estimated run time of thirty minutes. **Try creating a submission script
for running FastK on one of the HiFi data files. Check with one of the
instructors before running it.**

### FastK Output Files

Let’s look at the output files. Were they created in the place that you
specified in your submission script? How many files were created? Which
ones can be viewed from within VS Code or using the `less` command? How
can you view these output files? Do you need to view them? What does
FastK do, and why is it important for our pre-assembly steps?

## 3. GenomeScope

As you will have seen in class, GenomeScope provides us with some
crucial information about our data before we start assembling the
genome. For the full documentation, you can visit
<a href="https://github.com/schatzlab/genomescope">https://github.com/schatzlab/genomescope</a>.
Running GenomeScope requires you to load the FastK and the GenomeScope
modules. It also requires that you first convert the output of FastK
into the appropriate input file for Genomescope. This means that we will
have *two* commands in our submission script:

``` text
Histex -G <path_to_.hist_file> > <path_to_output_directory>/<species_prefix>.histo
```

\*AND\*\*

``` text
genomescope.R -i <path_to_.histo_file> -o <output_file_name> -k <k-mer_size> -p <known_ploidy_of_species>
```

For our data, use 16 CPUs and 8G memory, and an estimated run time of
one hour. Use a k-mer size of 31. We will talk about the ploidy of these
species in class. Go ahead and create a submission script for running
GenomeScope. Run it by one of your instructors if you’re feeling unsure.

## 4. SmudgePlot

Smudgeplots are another useful way of exploring our data before trying
to assemble it into contigs. It is particularly useful for verifying if
the data is consistent with the known (or estimated) ploidy of your
species. Unlike everything else we have used so far, SmudgePlot is not
available as a module on LiSC. Therefore, we need to use a **conda
environment**. Conda environments are super useful when working on a
cluster where you don’t have admin rights. In this situation (which is
the situation we are in working on LiSC), you can’t install your own
software, update software required by newer versions of
programs/pipelines, or run different pieces of software that have
different system requirements. Say you want to run a program called
`ksrates` (we’ll run this later in the course) on Tuesday morning and
later in the day you want to run a program called `EDTA` (we’ll also
talk about this one later). If `ksrates` requires python version 3.8 to
run properly and `EDTA` requires python 3.11, you may have a problem if
LiSC has set python 3.8 as the default python (and doesn’t have python
3.11 available). And this is just the beginning… It’s possible that
those programs have *many* other system-level dependencies to work and
they may not be the same! The administrators of LiSC are *very* helpful
and by loading modules we can be sure that most of the software we want
to run will work. However, sometimes you may be interested in running
less frequently used software or you want to use different (or older)
versions, in which case conda environments are super useful.

Wouldn’t it be helpful if you could have a contained “tool kit” that had
everything you need to run a specific piece of software regardless of
the system-level (LiSC-level) settings? Conda environments allow you to
do exactly that. You can create an isolated workspace that contains
everything you need for a specific task. When it comes time to do that
task, you can just load the conda environment that you created for that
task.

Thankfully, you don’t have to create your own conda environments for
this course, as either LiSC or your instructors have done it for you.
But, you do need to *load* Conda in your submission script, so that you
can activate the conda environment that you need. For us, that will look
like:

``` text
module load Conda
conda activate smudgeplot-0.5.3
module load FastK #please activate FastK again as smudgeplot relies on FastK scripts
```

With that in your script, you will be able to run SmudgePlot just like
you would if there was an available module.

Running Smudgeplot requires a lot of memory. You may remember from the
lecture that storing information about these k-mers and their
relationships to one another can be memory intensive. So to run
SmudgePlot we will use 16 CPUs, 64G of memory, and an estimated run time
of ten hours. This won’t run before the end of class today, so we will
look at the results next week. Of course, you can look at them before
then on your own! The command for SmudgePlot will use the following
options:

``` text
smudgeplot hetmers \
-L 7 \ # Low coverage cutoff (k-mers below 7x excluded)
-t <num_cpus> \
-o <output_directory> \ # You will create this first
--verbose \
<path_to_.ktab_file>

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot all -o <output_directory> <path_to_.smu_file>

```

**Create a submission script for running SmudgePlot. If you are feeling
unsure about it, check with your instructors before submitting it to
LiSC.**
