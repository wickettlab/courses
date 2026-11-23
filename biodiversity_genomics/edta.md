# Biodiversity Genomes

## Day 7: Annotating and Masking Transposable Elements

### Norman Wickett and Daria Shipilina

### 17 November, 2026

Before we can annotate your newly assembled and scaffolded genome sequences, we need to *mask* repetitive regions.
Today, we will focus on annotating and masking **transposable elements** using a program called EDTA (the **E**xtensive *de novo* **T**E **A**nnotator).
The full documentation for EDTA can be found here: <a href="https://github.com/oushujun/EDTA">https://github.com/oushujun/EDTA</a>.

## 1. Preparing your data for EDTA

In general, there is very little you need to do to prepare your data to run EDTA, as it simply takes a fasta file of your scaffolds as input.
However, EDTA can take a long time to run on your full assembly, so to test it in class we will run it only on one scaffold.
The first scaffold in your fasta file should be the longest, so let's print out the longest scaffold to its own file.
We can do this with the following bash one-liner:

``` text
awk '/^>/ { if (++n > 1) exit } 1' <your_scaffolds.fa> > <longest_scaffold.fa>
```

Now you should have a file named ``longest_scaffold.fa`` (or whatever you decided to call it) that you can use as the input for EDTA, which will run within our class period.

## 2. Running EDTA

Now let's try running EDTA.
To run this program, we must once again load a conda environment.
So the first thing you should do is load the Conda module.
Thankfully, LiSC has created an EDTA conda environment accessible to all LiSC users, so you can activate the environment by typing:

```text
conda activate edta
```

I like to do this directly from the command line just to test that everything is working before putting it into a submission script.
You will know it's working if you type:

```text
EDTA.pl -h
```

and you get a help menu. You'll notice that you are actually running a perl script (the .pl extension lets us know) that will link a bunch of other dependant programs, which is why having this installed as a conda environment is so useful.

The options we will use when running EDTA are:

```text
--genome <path_to_your_scaffold(s)>
--step all
--overwrite 1
--sensitive 1
--anno 1
--force 1
--threads 32
```

Take a minute to read through the help menu or documentation to familiarize yourself with what these options mean. **Go ahead and write a submission script to run EDTA on your longest scaffold.** I think if you use 32 cpus (should match the number of threads in your command), 64GB of memory, and four hours it should run before class ends today. For *Piper nigrum* you might want to increase the memory and time a bit, just to be sure.

Once your job has finished, you should have several files that include information about the identity and position of TEs in your scaffold, as well as a file that *masks* the TEs. In EDTA these will be indicated by Ns. You will also have a couple of figures that give you an overview of the TE landscape in your scaffold. Go ahead and download these files (they will be the .pdf files in the output) and have a look at them. We will discuss what they mean in class.

## 3. For Next Week

After you have successfully run EDTA on the longest scaffold, you can set it up to run on your full set of scaffolds. You'll need to adjust the expected resources, so increase the memory to 12GB and the time to two days. That should be enough!
