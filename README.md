A Guide to Multivariate Bayesian Analyses in Nursing Research
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 Repository Superceded

**Note:** As of December 2023, this repository has been superceded by a
new and improved version. The updated tutorial is now available at
[github.com/lwheinsberg/mvNUR](https://github.com/lwheinsberg/mvNUR).
Please visit the new repository for the latest content and improvements.

Thank you for your interest and understanding!

------------------------------------------------------------------------

OLD CONTENT RETAINED BELOW

# 2 Copyright information

Copyright 2023, University of Pittsburgh. All Rights Reserved. License:
[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# 3 Overview

*I have a Bayesian inference joke but the first three people I told it
to didn’t laugh and now I’m not so sure it’s funny.* - @JSEllenberg

This repository was created to support a presentation at the
International Society of Nurses in Genetics:

Heinsberg LW. Multivariate Bayesian Approaches for Analyzing Correlated
Phenotypes in Nursing Research. (Expert Lecturer Abstract, Podium).
Presented at the International Society of Nurses in Genetics, November
2023, Providence, Rhode Island.

which is currently being adapted for (hopeful) publication as a
manuscript entitled:

Heinsberg LW, Davis TS, Maher D, Bender CM, Conley YP, Weeks DE. A Guide
to Multivariate Bayesian Analyses in Nursing Research. In preparation
for submission to Biological Research for Nursing.

The goal of this repository is to provide detailed code and a fully
synthetic example data set to guide nurse scientists and other
researchers in conducting multivariate Bayesian analyses to examine
associations between correlated phenotypes and single nucleotide
polymorphisms (SNPs, i.e., genetic variants). This guide will focus on
the application of two multivariate Bayesian software programs: (1)
`bnlearn` and (2) `mvBIMBAM`.

**Approach 1:** `bnlearn` - bnlearn is an R package for learning the
graphical structure of Bayesian networks,estimating their parameters,
and performing some useful inference as described in:

“Multiple Quantitative Trait Analysis Using Bayesian Networks” by
Scutari, Howell, Balding, Mackay (Genetics, 2014).

and at:

[Link](https://www.bnlearn.com/)

Additional details about the method and its interpretation can be found
in `01_mvNUR_bnlearn.Rmd`.

**Approach 2:** `mvBIMBAM` - mvBIMBAM implements a terminal-based
Bayesian approach for genetic association analysis of multiple related
phenotypes, as described in:

Shim H, Chasman DI, Smith JD, Mora S, Ridker PM, Nickerson DA, Krauss
RM, Stephens M. A multivariate genome-wide association analysis of 10
LDL subfractions, and their response to statin treatment, in 1868
Caucasians. PLoS One. 2015 Apr 21;10(4):e0120758. doi:
10.1371/journal.pone.0120758. PMID: 25898129; PMCID: PMC4405269.

Stephens M. A unified framework for association analysis with multiple
related phenotypes. PLoS One. 2013 Jul 5;8(7):e65245. doi:
10.1371/journal.pone.0065245. Erratum in: PLoS One. 2019 Mar
19;14(3):e0213951. PMID: 23861737; PMCID: PMC3702528.

Additional details about the method and its interpretation can be found
in `02_mvNUR_mvBIMBAM.Rmd`.

While `bnlearn` is ran directly in R, `mvBIMBAM` is a terminal-based
program with no point/click desktop app. As such, please note that this
example code, particularly the mvBIMBAM approach, requires at least an
introductory understanding of R and Unix.

# 4 Installation

## 4.1 bnlearn

You can install bnlearn from
[CRAN](https://cran.r-project.org/web/packages/bnlearn/index.html)
using:

    install.packages("bnlearn")

or

    install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")

Check your R version if your are having any issues installing. As of
October 2023, bnlearn requires R version 4.2.0 or higher.

## 4.2 mvBIMBAM

mvBIMBAM can be installed from
[GitHub](https://github.com/heejungshim/mvBIMBAM). We have expanded on
these instructions below. THe following section assumes that you have an
introductory understanding of navigating your computer using the
terminal/Unix commands. If you are new to Unix navigation via the
terminal, please visit this
[link](https://macpaw.com/how-to/use-terminal-on-mac) for a brief
introduction.

The binary executable file for Linux is in the `mvBIMBAM/bin/` directory
while a zip file of the binary executable file, source code, and example
input files for Mac is in `mvBIMBAM\forMAC` directory. For transparency
and to provide as much help as possible, we have expanded on the Mac
instructions below.

### 4.2.1 Mac

The instructions for installing mvBIMBAM on a Mac are included here.

1)  Check for installation requirements.

Note that you need to have GSL (GNU Scientific Library) installed on
your machine before you can install mvBIMBAM. GSL is an open-source
software library that provides a wide range of mathematical and
statistical functions for scientific and numerical computing. If GSL is
not installed, it is simplest to install it using homebrew.

Before moving to the instructions below, using the terminal, check to
see if you have GSL installed on your machine by asking for the version
number:

    gsl-config --version

For example, on my machine, this command returns “2.7.1”.

If GSL is not installed, the simplest way to install it is via homebrew.
Homebrew is a popular package manager for macOS and Linux. It provides a
convenient way to install, manage, and update software packages and
libraries on your computer. Homebrew simplifies the process of setting
up development tools and applications.

First, install homebrew using the following terminal command:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

Confirm that it has been installed, and check for the path:

    brew --version
    brew --prefix

For example, on my machine, the first command returns “Homebrew 4.0.23”
and the second returns “/opt/homebrew”. Note that this tutorial was
developed on an Apple silicon machine. If you have an Apple Intel, the
default homebrew install path will be “usr/local/Cellar” or “usr/local”.

Now let’s move on to installing GSL:

    brew install gsl

Confirm install and path:

    gsl-config --version 
    gsl--config --prefix

Voila! (Hopefully :))

2)  Now let’s move back to the task in hand - installing mvBIMBAM
    itself! First, download the mvBIMBAM .zip file from
    [GitHub/heejungshim/mvBIMBAM/forMac](https://github.com/heejungshim/mvBIMBAM/tree/master/forMac).

3)  Next, unzip the BIMBAM_multi_pheno_MAC file by double clicking the
    zipped folder.

4)  Using the terminal, navigate to the newly unzipped mvBIMBAM folder
    (e.g., cd Desktop/BIMBAM_multi_pheno_MAC).

<!-- -->

    cd Desktop/BIMBAM_multi_pheno_MAC
    pwd 

For example, on my machine, when I use the `pwd` command (print working
directory), my folder is located at:
/Users/username/Desktop/BIMBAM_multi_pheno_MAC

5)  Now let’s attempt to run `make all` terminal command which will
    activate the software build automation tool and actually intall it
    on your machine.

<!-- -->

    make all

6)  *In case of terminal error in step 5*: The motivation behind step 6,
    which is NOT included in the source instructions for installing
    mvBIMBAM, is that when running the `make all` Unix command, I (and
    subsequently Dr. Tara Davis, the tutorial test user) received the
    following error:

<!-- -->

    cd src && /Applications/Xcode.app/Contents/Developer/usr/bin/make
    g++ -static-libgcc -DIMPUTATION   -O3 control.o fpmath.o indiv.o diploid.o haploid.o model.o param.o  fp.o -lm libgsl.a libgslcblas.a  -o bimbam
    clang: error: unsupported option '-static-libgcc'
    make[1]: *** [fp] Error 1
    make: *** [all] Error 2

which suggests there is an issue with the build process. If you
encounter this error, you can “hack” a solution to this install problem
by modifying the Makefile (located within the `src` folder) in a text
editor as detailed below. Note that there are separate instructions for
modifying the Makefile for an Apple Silicon (e.g., 2021 M1) vs. an Apple
Intel.

If you are not sure which Mac you have, you can click on the apple
symbol in the upper left corner of your machine. If you have an Apple
silicon, the Chip will be listed as Apple M1 or M2. If you have the
Apple Intel, you will see the processor listed as Intel Core i5, i7, or
similar.

\*\*\* APPLE INTEL INSTRUCTIONS \*\*\*

**CHANGE A** Add GSL Compiler and Linker Flags: In makefile located at
`/Users/username/Desktop/BIMBAM_multi_pheno_MAC/src/Makefile`, the
compiler and linker flags can be updated to include the necessary GSL
headers and libraries.

    CFLAGS += -I/usr/local/Cellar/gsl/2.7.1/include
    LDFLAGS += -L/usr/local/Cellar/gsl/2.7.1/lib

These flags include the GSL headers and specify the location of the GSL
libraries on my machine. Note that, depending on your computer itself
and the locations/versions of your GSL, these lines will likely need to
be adapted.

To check the version and path of your GSL, you can use the following
commands:

    gsl-config --prefix
    gsl-config --version 

Replace `usr/local/Cellar` with the location returned by the first
command and replace `2.7.1` with the version returned by the second
command.

If you do not have GSL installed, return to step 1.

**CHANGE B** Update LIBS Variable. Specifically, modify:

    LIBS += -lm libgsl.a libgslcblas.a

to:

    LIBS += -lm -lgsl -lgslcblas

This directly links the GSL libraries for math, GSL core, and GSL CBLAS.

**CHANGE C** Update Compilation Command: Remove the static-libgcc flag.
Specifically, modify this line:

    fp: fp.o $(OBJS); $(CC) -static-libgcc $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam

to

    fp: fp.o $(OBJS); $(CC) $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam

The -static-libgcc flag is related to linking the GCC (GNU Compiler
Collection) runtime libraries statically. However, this flag is not
supported by the Clang compiler, which is commonly used on macOS
(including Apple Silicon M1 Macs). By removing it, the build rule
becomes compatible with Clang and more common macOS build setups.

A copy of the original Makefile, modified Makefile, and word document
comparing the two can be found in the `Install_Problems` folder.

Note that this problem has been opened as an issue on the mvBIMBAM
GitHub page <https://github.com/heejungshim/mvBIMBAM/issues/3>. Please
check there for updates from the program architects.

\*\*\* APPLE SILICON \*\*\*

**CHANGE A** Add GSL Compiler and Linker Flags: In makefile located at
`/Users/username/Desktop/BIMBAM_multi_pheno_MAC/src/Makefile`, the
compiler and linker flags can be updated to include the necessary GSL
headers and libraries.

    CFLAGS += -I/opt/homebrew/Cellar/gsl/2.7.1/include
    LDFLAGS += -L/opt/homebrew/Cellar/gsl/2.7.1/lib

These flags include the GSL headers and specify the location of the GSL
libraries on my machine. Note that, depending on your computer itself
and the locations/versions of your GSL, these lines will likely need to
be adapted.

To check the version and path of your GSL, you can use the following
commands:

    gsl-config --prefix
    gsl-config --version 

Replace `opt/homebrew/Cellar` with the location returned by the first
command and replace `2.7.1` with the version returned by the second
command.

If you do not have GSL installed, return to step 1.

**CHANGE B** Update LIBS Variable: Specifically, modify:

    LIBS += -lm libgsl.a libgslcblas.a

to:

    LIBS += -lm -lgsl -lgslcblas

This directly links the GSL libraries for math, GSL core, and GSL CBLAS.

**CHANGE C** Update Compilation Command: The compilation command for the
“fp” target should be updated to remove the static-libgcc flag and add a
linker. Specifically, modify:

    fp: fp.o $(OBJS); $(CC) -static-libgcc $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam

to

    fp: fp.o $(OBJS); $(CC) $(CFLAGS) $(OBJS) fp.o $(LDFLAGS) $(LIBS) -o bimbam

The -static-libgcc flag is related to linking the GCC (GNU Compiler
Collection) runtime libraries statically. However, this flag is not
supported by the Clang compiler, which is commonly used on macOS
(including Apple Silicon M1 Macs). By removing it, the build rule
becomes compatible with Clang and more common macOS build setups.

Note also that Apple Silicon also requires the extra `$(LDFLAGS)` which
is used to specify linker-related flags and options that may be specific
to the architecture of M1 (not needed by Intel).

A copy of the original Makefile, modified Makefile, and word document
comparing the two can be found in the `Install_Problems` folder.

Note that this problem has been opened as an issue on the mvBIMBAM
GitHub page <https://github.com/heejungshim/mvBIMBAM/issues/3>. Please
check there for updates from the program architects.

7)  Run Unix command `make all` to create the binary executable file in
    the mvBIMBAM directory. Hint: If you already attempted to use
    `make all` and it failed, use `make clean` first which will remove
    artifacts from the earlier attempts.

<!-- -->

    make clean 
    make all

I feel like there should be some confetti that shoots out of your
computer to help you celebrate successful install of mvBIMBAM but,
unfortunately, that can’t be programmed just yet and it is a little
difficult to tell if the program installed correctly. To confirm
successful install, follow instructions in step 8.

8)  Move the executable file to system’s PATH environment variable, such
    as /usr/local/bin, so the program can be called from any location
    (vs. only the folder where the executable file resides). To do this,
    navigate to the location of the Unix executable file via the
    terminal and then use the move (`mv`) command.

<!-- -->

    sudo mv bimbam /usr/local/bin/

Note that Confirm install/move worked:

    which bimbam

For example, on my machine, bimbam is located at `/usr/local/bin/bimbam`
now.

# 5 Synthetic data set overview

To facilitate hands-on learning, a synthetic data set was created for
use with this example code. The synthetic data set was created from a
study of women with breast cancer (BrCa). This data set contains no real
data, but mimics the correlation structure of the true data set. The
original data set contained data for 110 participants and included data
for 6 symptoms, 3 demographic variables, and 6 SNPs. We created a
synthetic version of this data set with fake data for 770 participants
for use in this tutorial.

Data set name: `BrCa_synthetic.csv`

Variable descriptions:

Symptoms/Outcomes  
-`ID`: Participant ID  
-`EMO_tscore`: Anxiety (PROMIS Anxiety T Score)  
-`bdito`: Depression (Beck’s Depression Inventory)  
-`FAT_tscore`: Fatigue (PROMIS Fatigue Score)  
-`paohcif`: Cognitive Function (Patient’s Assessment of Own Functioning
Inventory - Cognitive Function)  
-`EPSscore`: Sleepiness (Epworth Daytime Sleepiness)  
-`worst_pain`: Pain (Self-reported numeric rating of “average pain”)

Demographics/Covariates  
-`age`: Age in years  
-`education`: Education in years  
-`race`: Race (Self-identified, 0=White, 1=Black)

SNPs/Predictors  
-rs \#s (SNPs)

``` r
# Read in a mapping of variable name to informative label
dict <- read.csv("data/BrCa_synthetic_SimpleDict.csv")
pander::pander(dict, "Data dictionary") 
```

|  Variable  |       Label        |                                                                                                                                                                                                                                        Description                                                                                                                                                                                                                                         |
|:----------:|:------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| Emo_tscore |      Anxiety       | PROMIS Emotional Distress Anxiety (Short Form 8a) - 8 items measuring emotional distress and anxiety (panic, fearfulness, worry, dread, tension, nervousness, restlessness, and somatic symptoms including racing heart and dizziness) in the last 7 days. This instrument generates T-scores, which are standard scores with a mean of 50 and standard deviation of 10 in the U.S. general population. Higher T-scores indicate worse distress and anxiety. T-scores range from 10 to 90. |
|   bdito    |     Depression     | PROMIS Emotional Distress Anxiety (Short Form 8a) - 8 items measuring emotional distress and anxiety (panic, fearfulness, worry, dread, tension,nervousness, restlessness, and somatic symptoms including racing heart and dizziness) in the last 7 days. This instrument generates T-scores, which are standard scores with a mean of 50 and standard deviation of 10 in the U.S. general population. Higher T-scores indicate worse distress and anxiety. T-scores range from 10 to 90.  |
| FAT_tscore |      Fatigue       |                                                         PROMIS Fatigue (Short Form 8a) - 8 items measuring the experience of fatigue and the interference of fatigue on daily activities over the past 7 days. This instrument generates T-scores, which are standard scores with a mean of 50 and standard deviation of 10 in the U.S. general population. Higher T-scores indicate worse fatigue. T-scores range from 10 to 90.                                                          |
|  paohcif   | Cognitive function |                    Patient Assessment of Own Functioning Inventory (PAOFI) - 31 items measuring self-reported ‘recent’ function consisting of 5 dimensions: memory, language and communication, use of hands, sensory-perceptual, higher level cognitive and intellectual functioning. Higher scores indicate worse perceived functioning. Scores range from 0 to 155. For this study, we pulled the cognitive function dimension for which scores range from 0 to 45.                     |
|  EPSscore  |     Sleepiness     |                                                                                                                                                         Epworth Sleepiness Scale - 8 item scale measuring self-reported ‘recent’ daytime sleepiness. Higher score indicate greater daytime sleepiness. Scores range from 0 to 24.                                                                                                                                                          |
|    pain    |        Pain        |                                                 Brief Pain Inventory (Short Form) - a widely used clinical tool for assessing pain (worst pain in the last 24 hours $$0-10$$, least pain in the last 24 hours $$0-10$$, pain on average $$0-10$$, paing at the time of the interview $$0-10$$, pain relief from treatment(s) in the last 24 hours $$0%-100%$$). Higher scores indicate higher levels of pain, each item ranges from 0-10.                                                  |
|   rs4880   |     Variant 1      |                                                                                                                                                                                                                          SOD2 gene, hg19 position chr6: 160113872                                                                                                                                                                                                                          |
| rs5746136  |     Variant 2      |                                                                                                                                                                                                                          SOD2 gene, hg19 position chr6: 160103084                                                                                                                                                                                                                          |
| rs1041740  |     Variant 3      |                                                                                                                                                                                                                          SOD1 gene, hg19 position chr21: 33040162                                                                                                                                                                                                                          |
| rs10432782 |     Variant 4      |                                                                                                                                                                                                                          SOD1 gene, hg19 position chr21: 33036391                                                                                                                                                                                                                          |
| rs4135225  |     Variant 5      |                                                                                                                                                                                                                          TXN gene, hg19 position chr9: 113006691                                                                                                                                                                                                                           |
| rs7522705  |     Variant 6      |                                                                                                                                                                                                                          PRDX1 gene, hg19 position chr1: 45992300                                                                                                                                                                                                                          |
|    race    |        Race        |                                                                                                                                                                                                                           Self-identified race, 0=White, 1=Black                                                                                                                                                                                                                           |
|    age     |        Age         |                                                                                                                                                                                                                                        Age in years                                                                                                                                                                                                                                        |
| education  |     Education      |                                                                                                                                                                                                                              Self-reported years of education                                                                                                                                                                                                                              |

Data dictionary (continued below)

|          Type          |                                                                                                     Citation                                                                                                      |
|:----------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|        Numeric         | Pilkonis, P.A. et al. (2011) Item banks for measureing emotional distress from the patient-reported outcomes measurement information system (PROMIS): Depression, Anxiety, and Anger, Assesssment, 18(3), 263-283 |
|        Numeric         |                                                         Beck, A.T. et al. (1996) Beck Depression Inventory-II. San Antonio: The Psychological Corporation                                                         |
|        Numeric         |                              Cell, D. et al. (2016) PROMIS fatigue item bank had clinical validity across diverse chronic conditions. Journal of Clinical Epidemiology, 73, 128-134                               |
|        Numeric         |                                                 Chelune, G. J. et al. (1986) Neuropsychological and personality correlates of patients’ complaints of disability.                                                 |
|        Numeric         |                                   Cleeland, C.S. et al. (1994). Pain assessment: Global use of the Brief Pain Inventory. Annals, Academy of Medicine, Sigapore, 23(2), 129-138                                    |
|        Numeric         |                                                                         Cleeland, C.S. et al. (2009). The Brief Pain Inventory User Guide                                                                         |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
| Numeric, encoded value |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |

# 6 Examples

To prepare the data before running the examples, see
`00_mvNUR_DataPrep.Rmd`. After running that markdown, examples of each
method can be found in `01_mvNUR_bnlearn.Rmd` and
`02_mvNUR_mvBIMBAM.Rmd`.

# 7 Contact information

If you have any questions or comments, please feel free to contact me!

Lacey W. Heinsberg, PhD, RN: <law145@pitt.edu>

# 8 Acknowledgments

I’d like to express my gratitude to the following for their support and
contributions to this repository:

- Support from the National Institutes of Health under award number
  K99HD107030 made this project possible, and for that, I’m truly
  grateful.
- Special thanks to Dr. Tara Davis and Mr. Dylan Maher for being “test
  users” and providing their invaluable feedback on this guide.
- My deepest gratitude to Dr. Daniel Weeks for … everything. If it
  weren’t for his guidance and inspiration, I would likely never have
  ventured into the mystifying world of Bayesian statistics. —-\>
  Disclaimer: As of 10-31-2023, the public debut of this repository for
  ISONG, Dan has not had a chance to work his magic yet. Knowing the
  immeasurable impacts Dan makes on anything he touches, I’m certain
  that this guide will ascend to new heights once he provides his
  feedback. Stay tuned for a new and improved version of this guide as
  we move our accompanying tutorial manuscript forward soon.
