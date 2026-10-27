# Biodiversity Genomes

## Day 2: Submitting Jobs to LiSC and Basic Sequence File Stats

### Norman Wickett and Daria Shipilina

### 13 October, 2026

Today we will be learning how to write small scripts that submit jobs to
LiSC. From last week you should remember that we submit jobs to LiSC
without running them directly. Instead, there is a system on LiSC that
runs the jobs in a way that uses resources most effectively. This is
super important because many users are running jobs on LiSC
simultaneously, so there is no way one individual can anticipate how
their job might impact other users. The system that LiSC uses to
schedule jobs is called **Slurm**.

## 1. Submission Scripts

LiSC is a high-performance computing cluster, which means it has A LOT
of power to run computationally intensive jobs, and run many of them at
the same time. If everyone tried to run their jobs directly from the
command line, LiSC wouldn’t be able control when they start and how many
can run at once, which is why it uses Slurm. Without this control, the
resources of LiSC can slow down and become less useful for everyone who
uses the HPCC. Instead of just typing out the command for the program
you want to run directly in the terminal, you must put it into something
called a submission script. Inside that submission script are estimates
about how much time and resources will be needed for your job,
**modules** (we’ll talk about these a bit more later) that are needed to
run specific programs, and the specific command(s) to run your job. When
you submit your submission script to Slurm, it evaluates the available
resources and how your job will impact those resources, and schedules
the best time for it to run.

A very simple submission script looks like this:

``` text

#!/bin/bash
#
#SBATCH --job-name=test_trim
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --partition=basic
#SBATCH --mail-type=ALL
#SBATCH --time=0-02:00:00

module load <your_module_here>

command used to run your program
```

So what do all of those lines mean?

``` text
#!/bin/bash -> this means that the script can be interpreted by bash
#
#SBATCH --job-name=test_trim -> a name so that you can identify the job
#SBATCH --cpus-per-task=16 -> how many processors are needed
#SBATCH --mem=32GB -> how much memory (RAM) is needed
#SBATCH --partition=basic -> just run on the "normal" LiSC
#SBATCH --mail-type=ALL -> get emails when anything happens with your job
#SBATCH --time=0-02:00:00 -> how much time is required (days,hours,minutes,seconds)
```

It can be hard to know exactly what to put for some of these values, so
I tend to overestimateto make sure the job will complete without any
issues. We will discuss what these should be - you don’t need to figure
it out by yourself!

What about those other lines?

``` text
module load <your_module_here>
```

To run something on LiSC you must first load a module so that your
program and all its dependent programs are ready to go. For Trimmomatic,
you would load the `Trimmomatic` module.

``` text
command used to run your program
```

The exact way you run each program will be different for every program,
so there’s no universal line to write here. You’ll see this as we run
the programs below. Often, you’ll want to run the same thing on multiple
files. In that case we can use a `for loop` to go through each file and
run your command. You’ll see this in the script for Trimmomatic…

But first you’re probably wondering how you actually submit the
submission script! It s pretty easy and just takes the `sbatch` command.
Let’s say you have a submission script called `trimmomatic.sh`. To
submit that, you just type:

``` text
sbatch path/to/your/script/trimmomatic.sh
```

Hit enter, and you’ll get a message that it has been submitted. Now sit
back and wait for the email that says your job has started!

## 2. Modules

We briefly mentioned modules above. For many commonly-used programs,
LiSC makes them available via modules. Before loading a module in a
submission script, it can be good to check how the module is written
(LiSC is case-sensitive) and verify the commands once the module is
loaded. Let’s say you want to check if there is a module available for
FastQC. Go ahead and try this:

``` text
module avail fastqc
```

Hopefully you see something like this:

``` text
-------------- /lisc/opt/sw/modules/all ----------------
   FastQC/0.11.9-GCCcore-12.3.0-Java-11    FastQC/0.12.1-GCCcore-14.3.0-Java-21
   FastQC/0.12.1-GCCcore-12.3.0-Java-11    FastQC/0.13.0-GCCcore-14.3.0-Java-21 (L,D)

  Where:
   L:  Module is loaded
   D:  Default Module

If the avail list is too long consider trying:

"module --default avail" or "ml -d av" to just list the default modules.
"module overview" or "ml ov" to display the number of modules for each name.

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules matching 
any of the "keys".
```

Those first two lines (e.g., `FastQC/0.11.9-GCCcore-12.3.0-Java-11`)
tell you where the program actually lives on LiSC. You don’t need that
information but you do need to notice how FastQC is written: `FastQC` -
this means that when you load the module you must type it exactly like
that. What happens if you type:

``` text
module load fastqc
```

Hopefully you’ll see a line that says, “Lmod has detected the following
error:”. This means that a module was not loaded. Now try:

``` text
module load FastQC
```

Now there is no error message and the module has been succesffuly
loaded. Once a module is loaded it is possible to actually use the
program. One super useful thing you can do *without* a submission script
is view the help menu for a program. Help menus are often accessed using
the `-h` option/flag. I realize this is a little confusing, but even
though LiSC is case sensitive and requires you to type `FastQC` exactly
as written, the actual command for running FastQC is simply `fastqc`
(there is a difference between the *module name* and the *actual
command*). To view the help menu of FastQC, type:

``` text
fastqc -h
```

What do you see?

Help menus do two important things: they give you the **usage** of the
command (the actual syntax of how you write the command, tell it which
files to use etc.) and all (or most) **options/flags** that you can use
with the command. Notice the first line of the FastQC help menu is:

``` text
usage: fastqc [-h] [-v] [--outdir OUTDIR] [-c CONTAMINANTS] [-a ADAPTERS] [-l LIMITS]
          [-d DIR] [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
              [--dup_length DUP_LENGTH] [--memory MEMORY] [-t THREADS] [-q]
          [--fast5 FAST5] [--phred64] [-j JAVA] [--nogroup] [--png] [--extract]
          [--delete] [-f FORMAT] [--debug] [files ...]
```

That tells you that after you type the command `fastqc` you give it the
options you want to use, followed by the files that you want to run
using FastQC. Following the usage statement is an explanation of all of
the options. Sometimes a program will tell you which flags are
*mandatory* and which flags are *optional*. But remember that we don’t
want to do this directly from the command line, but using a submission
script instead.

## 3. Submission Script Practice: FASTA and FASTQ stats

You will have already gone over what FASTA and FASTQ files look like in
class so you know that these files are simply text files that contain
some identifiers (lines that start with either `>` or `@`) followed by
DNA sequences (these can either be nucleotide sequences comprising A, C,
T, and G or amino acid sequences). For many reasons, however, you might
want to know some basic information about your sequence files before
doing anything with them (like trimming poor quality bases, or even
assembling the reads into contigs). For Illumina data, it is standard
practice to run something called **FastQC**
(<https://github.com/s-andrews/fastqc>). FastQC will give you an overall
assessment of the quality of your data over all positions in the
sequence reads and highlight anything suspicious with your libraries.
Luckily, you already know how to load the module for FastQC and some
basics about how to use it!

Let’s run FastQC on some Illumina data using a submission script that
you write yourself. The data we will be using are RNA-Seq data (sequence
data from the expressed genes) from *Dianthus*. Later in this course we
will be using these data to help annotate the *Dianthus* genome but for
now we can use these data to understand both how submission scripts
work, and how FastQC works.

To write our script we need to know where the data are located (the path
to the data), where the output should go, and which options we are
using. We’ll talk about the options in class, and why we are using them.
Because the default location of the output files is the same as the
input, you’ll also want to create a new directory in your personal
course directory for the output, otherwise everyone will be trying to
write the same output files at the same time in the same place. Go ahead
and make a directory called `fastqc_output`. Because of how FastQC write
files as it is working, you’ll also need to create a directory called
`fastqc_temp`.

Now let’s actually write the script. To creat a file, navigate to your
personal course directory in VS code using the explorer. Click on the
`New File...` icon and name the new file `fastqc_submission.sh`. You
should now see your new file in the file window, where you can edit it
as follows:

``` text
#!/bin/bash
#
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --partition=basic
#SBATCH --mail-type=ALL
#SBATCH --time=0-01:00:00

module load FastQC

fastqc --outdir /lisc/data/scratch/course/2026w300485/<your_directory>/fastqc_output \
--dir /lisc/data/scratch/course/2026w300485/<your_directory>/fastqc_temp -t 16 \
/lisc/data/scratch/course/2026w300485/course_data/RNA_Seq/Dsuper_6B_MKRN260064942-1A_23M3G5LT4_L1_1.fq.gz \
/lisc/data/scratch/course/2026w300485/course_data/RNA_Seq/Dsuper_6B_MKRN260064942-1A_23M3G5LT4_L1_2.fq.gz
```

Before we submit this file, **don’t forget to save it!** Now, to submit
the file, simply type:

``` text
sbatch fastqc_submission.sh
```

Once Slurm starts running your job, you will receive and email. You will
also receive an email once it has finished or if your job fails. In
class, we will talk about the results, as well as how to download and
view the images produced by FastQC.
