As you learn more about genome analysis in the course, you will
understand that processing and analysing enormous amounts of DNA
sequence data requires enormous amounts of computer power. To access
this computer power, a high-performance computing cluster (like LiSC) is
almost always required. In the previous exercises you should have
configured your access to LiSC through VS Code. Thankfully, VS Code
allows you to access and navigate LiSC in a way that is somewhat
familiar. However, to use LiSC to its full potential it is important to
understand how to navigate such Linux-based systems from the command
line. The command line is simply a way of interacting with a computer
operating system through text-based commands, rather than the “point and
click” system that most of us are used to. Command line computing is a
fundamental skill in bioinformatics and is used daily by scientists
studying biodiversity genomics. The goal of this exercise is to provide
you with some basic command line skills, which you will build upon over
the next few months of this course. Later, you will learn more about how
to use your command line skills to run analyses (like genome assembly,
genome annotation etc.) on LiSC.

Using the command line to interact with a computer or compute cluster
requires something that is called either a **command line interface**
(pretty self explanatory) or a **shell**. Shells consist of two “parts”
– the program displays your commands and then the “language” of the
shell that actually does what you type. The programs that display what
you type come in different “flavours,” like **terminal** on a Mac or the
**windows terminal** on a Windows PC. These programs then use a
“language” to interpret what you type. Probably the most common shell
“language” is called **bash** or some variation of bash. Thankfully, in
this class we are all using VS Code so we don’t all have to learn an
interface and shell language that is specific to each of our computers
and may be different from other students in the class.

When you open a terminal in VS code, it will usually open as a panel at
the bottom of your VS Code window. In this document, commands that enter
in the VS Code terminal will be represented by text like this:

``` text
Specific commands will be written in text like this.
```

If you ever need to enter your own information, instead of typing
exactly what is in the command, the information you need to provide will
be inside angle brackets, \< and \>. For example, if I want you do run a
command on your file called “my_data.txt” and I ask you to type:

``` text
command -i <your_file_here> -p
```

You would NOT type

``` text
command -i <your_file_here> -p
```

And instead you would type

``` text
command -i my_data.txt -p
```

## 1. Navigating a File System

Open VS Code and connect to LiSC (click on the offset \>\< symbol in the
lower left corner of the window). A new window may pop up as a new
connected window, and you will be asked to enter your password. If there
is no TERMINAL window open (you will know that there is one open because
there will be a window with the word TERMINAL underlined in blue) open a
new one by going to Terminal \> New Terminal. All of the instructions
below will be typed in that terminal window.

One of the first things we learn in this class is how to navigate a file
system using the command line. File systems are simply how files are
organized on your operating system. Two common operating systems are
Windows or MacOS. In these operating systems, files are organized into
folders, and you can see what is in those folders by double clicking on
them, or you can move them around by clicking on them and “physically”
move them from one folder to another. We can do the same thing when we
interact with a Linux file system, but instead of clicking and dragging,
we use text commands. In Linux, folders are called directories, so we
will use “directory” from now on. Using this system is particularly
useful for HPPCs because people can connect from many different
operating systems and everyone can use this “universal language” to
navigate on the remote server.

Now that you know why we use the command line, let’s try out some
commands.

The first command to try is:

``` text
pwd    
```

This command represents “print working directory” and essentially tells
you “where you are” currently in the file system. The output should look
something like this:

``` text
/lisc/home/user/wickett
```

By default, when you log into LiSC, you will be “in” your home
directory. In the example above, I’ve logged in as myself, so I am in
the `wickett` directory.

Go ahead and try it. What does it show you? Hopefull it will look like:

``` text
/lisc/home/user/[your user name]
```

You should notice that your location in the file system is given by a
series of words separated by forward slashes. Each word is a directory
that is nested within the directory to its left. If you were using an
operating system like Windows or MacOS, you would have a folder called
**lisc**, in which there would be a folder called **home**, which itself
would contain a folder called **user**, and in that folder there would
be a folder named after you. Just like in Windows or MacOS, there can be
many folders within the folder that is “one level up.” So in our example
there could be many folders named after many different people, all
existing within the **user** folder.

The people who run LiSC are super helpful and, in addition to created an
account for each student in the course, they created a course directory.
In a minute, you will learn how to move to different directories but
`pwd` will always show you your current directory, that is, your current
location in the file system. For example, if I was in our course
directory it would look something like this:

``` text
/lisc/data/scratch/course/2026w300485/sequences
```

This means that I am currently in the directory `sequences` inside our
course directory. We will be doing all of our work in the directory that
was set up for this course - `/lisc/data/scratch/course/2026w300485/` -
so it will be good to get used to that location.

Just to reinforce the idea about how directories are nested within
directories that are “one level up,” the `sequences` directory is nested
within the `2026w300485` directory, which is nested within the `course`
directory, which is nested within the `scratch` directory, and so on.
When using a graphical interface like Windows or MacOS, you would click
on `course` and see all of the directories contained within, then you
would click on `2026w300485` and see what is in that directory. And so
on, and so on…

This long string of directories is called the **PATH**. This is a super
useful thing to know, because many analyses require that you provide the
PATH to input data or where output data should be created. Whenever we
talk about the “full PATH,” this is what is meant. It is also possible
to use a “relative PATH,” but we’ll deal with that in a minute.

Remember that by default you will be in your home directory when you
login to LiSC. For example, when I log into LiSC I would be in
`/lisc/home/user/wickett/`. But, you will also remember that we will do
all of our work in the course directory. That means that I have to move
from my home directory into the course directory. Since we are not using
a graphical interface for navigating our file system, if I want to move
to another directory I can use:

``` text
cd
```

This represents “change directory.” For example, if I want to move to
our course directory I would type:

``` text
cd /lisc/data/scratch/course/2026w300485/
```

If I typed `pwd` afterwards, it should tell me that I am in our course
directory by returning the full PATH to my current location in the file
system. When the output shows up on the screen after you type a command,
we say that the output of the command is, “printed to the screen.”

Go ahead and move into our course directory using `cd`. Once you have
done this, it would be useful to see what exists inside that directory.
Are there more directories? Are there any files? The way we can view the
contents of the directory is with:

``` text
ls
```

Try typing `ls`. What do you see printed to the screen? At the very
least you should see something that looks like:

``` text
course_data envs README.txt
```

That simply means that there are two directories, `course_data` and
`envs`, and a text file called `README.txt`. Depending on how different
systems are set up, you may see that the directories are printed in a
different color than the text file. This can be helpful for
understanding the contents of a directory.

## 2. Working with Files (and directories and executables)

This is a bit of computer science “philosophy,” and it may be a bit
weird to think about at first, but I think it helps understand
Unix/Linux file systems (like the one used on LiSC). If this feels a
little strange and confusing, don’t get discouraged! We will discuss
this and go over everything in class and you can ask questions any time.
Unix/Linux essentially treats *everything* as a file. Directories,
scripts, “programs,” DNA sequences etc. etc. are all files. What makes
them different is something called **permissions**. Permissions are a
code that is given to a file that provides and identity or function.

To understand how permissions work, it helps to see what they look like
on the screen. To do that we can use the `ls` command that prints the
contents of a directory to the screen, which we learned about earlier.
If we add something after `ls` to make the command more specific, we are
passing **arguments** to the `ls` command. Arguments can be
**options/flags** or positional arguments. Options or flags are passed
to the command using either one or two dashes. Usually, if the argument
is a full word, you use two dashes, but if it is abbreviated, it uses a
single dash.

For example, if you wanted a command to print the long version of an
output, you could pass the flag “long” to the command. This could either
be passed as `--long` or `-l`. The command `ls` is a great example of
this. Try typing the following:

``` text
ls -l
```

Hopefully what is printed to the screen looks a little different than
before!

``` text
drwxr-xr-x 2 wickett login 4096 Oct  1 13:19 course_data
drwxr-xr-x 2 wickett login 4096 Oct  1 13:19 envs
-rw-r--r-- 1 wickett login  113 Oct  1 13:20 README.txt
```

When you use `ls -l` to view the contents of your current
location/directory, you get information about the owner of the files,
the size of the files, when the files were created, and then that
10-character string of letters and dashes at the very beginning. The
first of those characters tells you what kind of file it is. If there is
a `d` that means it’s a directory, so you could use `cd` to move into it
and then see the files contained therein. If the first character is a
dash, `-`, that means it’s a regular file. Regular files are generally
text files that will be read by some program, and the program interprets
the text. We will look at files later, and hopefully it will make more
sense. The last nine characters are read in groups of three and specify
who is allowed to view (read) or edit (write) the file (or run the file
if it is a script). We don’t really need to worry about that for now.

You will often be asked to do something with a file, so let’s practice a
few standard commands. Since we are working in a shared space, it is
important that you don’t do anything to files that need to be accessed
by other people (this is usually not possible due to those permissions I
mentioned above, but it’s still good practice to be careful). To do
this, the first thing you should do is create a directory in which you
will do all of your work. The directory should be named as your last
name. The command to make a directory is `mkdir` - makes sense, right!

Make sure you are in the course directory by typing `pwd`. You should
see:

``` text
/lisc/data/scratch/course/2026w300485
```

To make your own directory, type:

``` text
mkidr <your_last_name>
```

Don’t forget that you need to actually type your last name instead of
`<your_last_name>`!! Now that you’ve created a directory, go ahead and
use `cd` to move into that directory. Verify that you are there by
typing `pwd` and checking that the path is correct. If you want, you
could type `ls` and you should notice there is nothing in your
directory. Yet.

The next several sections of this exercise will be a very brief overview
of some standard commands that you will often use. We will practice them
in class, so please ask questions if something doesn’t make sense!

## 3. Relative Paths

Okay, okay. I know I said we would practice some commands, but first
it’s important to learn about **relative paths**. There are two
important ways we can use a kind of shorthand to refer to the
location/path of files. They are a single dot or two dots:

``` text

.

or

..

```

The single dot is a way of specifying your current location. So if you
wanted to do something in your current directory and need to specify the
path, you can replace the full path with just that single dot. The
double dot is super useful and it points to one directory *above* your
current directory. So if you are in the directory you just created,
using the double dot refers to the course directory in which your
directory is contained. This is what I mean by “one level above.” *Why
is this important?* Well, if you wanted to, say, copy a file from the
course directory into your own directory, instead of typing the full
path to the course directory, you could just use `..` and for putting it
in your directory (assuming that’s where you currently are) you can just
use `.` - you’ll get used to this, I promise!

## 4. Copying Files

Copying a file from one location to another is super useful, especially
for this course. For example, there may be a file with DNA sequence data
in the `course_data` directory. Instead of everyone working on the file
in that location, you can copy the file into your own directory. Let’s
give it a try! There is a file in the `course_data` directory called
`fake_data.fa`. To copy that file into your own directory using relative
paths (you should be in your own directory as you do this), type:

``` text
cp ../course_data/fake_data.fa .
```

That is saying “copy from the directory called `course_data` that is one
level above my current location, a file called `fake_data.fa` and put it
in my current location. There may also be times when you want to copy an
entire directory, along with its contents. To do this you use
the”recursive” flag, `-r`. So to try that you could type:

``` text
cp -r ../course_data/fake_directory/ .
```

## 5. Moving Files

Moving files is a little different that copying files becuase when you
move a file you are not leaving a copy in the original location. With
the previous example you could use the same command but replace `cp`
with `mv`. The `mv`command means “move” (makes sense!). But if you moved
the `fake_data.fa` file to your directory, it would no longer be in the
`course_data` directory, so it’s important to be careful with this
command.

The move command is also useful for another reason: *it will rename
files without moving them* (or move files *and* rename them
simultaneously). One way to think about this is that you can “move” a
file from one name to another, even if it’s location changes. Or you
could move a file from one location to another and specify that its name
changes as well.

**Create a file in your own directory using VS Code. Name it with your
name and give it the `txt` file extension.** You can write whatever you
want in this file, even something like, “Hello, my name is Norm.” It’s
not super important what you write. If you’re not sure how to creat a
file, we will go over it in class. **Create a new directory within your
directory as well. Call this new directory `test_directory`.**

Now, let’s say I just created a file called `norm.txt` in my own
directory. If I want to move it into the `test_directory` I would type:

``` text
mv ./norm.txt ./test_directory/ ## Your file won't be called norm.txt of course! ##
```

Check to make sure you understand the use of relative paths in that
command. Now that the file is in the `test_directory` you could also
rename it without moving it:

``` text
mv ./test_directory/norm.txt ./test_directory/new_norm.txt
```

Use `cd` to move into the `test_directory` and view the files with the
`ls` command. You should see that the file with its original name
doesn’t exist anymore, but it exists with its new name. It’s important
to note that you can also rename files when you copy them in the same
way we have done here. If, for example, you were copying the
`fake_data.fa` file that we did previously, you could also give the copy
a new name by specifying the new name in the command.

## 6. Removing Files

We’re not going to spend too much time talking about removing files, but
you can remove files with the `rm` command. For example, you could
remove the `norm.txt` file by typing:

``` text
rm norm.txt
```

It is **VERY IMPORTANT** to be careful with the remove command though.
When you normally delete a file on your computer, you could retrieve it
from the recycling or trash bin if you made a mistake. But the remove
command will permanently delete the file. While LiSC has some backup
systems in place, it is best practice to use the remove command
carefully.

## 7. Viewing Files

One great thing about VS Code is that it’s easy to view files right
there in the file viewer panel. However, there may be times when you
want to view files the “old fashioned way” from the command line. One
common way to quickly view the contents of a file is to use the `less`
command. Try typing this in the terminal:

``` text
less <your_file_name>.txt
```

What happened? If it worked, it should display the contents of the file
on the screen, but do it “outside” of the terminal prompt. To get back
to the terminal, you can simply type `q`. This is convenient when you
just want to check whether a file contains what you expect. Two other
useful commands for viewing the contents of a file are `head` and
`tail`, which simply print the first or last 10 lines of a file to the
screen. You could try this with the `fake_data.fa` file that you copied
into your directory:

``` text
head ./fake_data.fa
```

or

``` text
tail ./fake_data.fa
```

What do you see? Does it make sense?
