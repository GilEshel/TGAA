# TGAA
A Python program to analyze multiple sequence alignment and compare amino acid properties between two groups of species

Key assumptions:
1. Alignment file is in a fasta format.
2. Each sequence header contains the species name (I use a short species label*, e.g. Aratha, instead of Arabidopsis_thaliana), followed by a '#' symbol and a sequence id (e.g. '>Aratha#AT5G28770')
3. You define two groups of species using a species groups file (e.g. Example_files/species_groups_file.txt)

*If you give short species labels instead of full species names, you can provide (using the '-s' option) a file like Example_files/species_labels_names.txt with the full species name per short species label - for the outputs... 

The program will analyze the alignment file, given the two defined groups and will:
1. Calculate the sum of pairs alignment scores and hydrophobicity per group, between-groups and for both groups - per site
2. Calculate the sum of pairs alignment scores and hydrophobicity per group, between-groups and for both groups - per sliding window (average or median of 10 sites, by default) - use the '-sw' option to set sliding window size.
3. Will compare the amino acid composition and properties (e.g. physico-chemical) within and between groups, to find sites that diverged between the two groups (and are conserved within groups).


## Running TGAA
\# You can invoke the help option to see all required and optional arguments:

**python tgaa.py -h**

\# Running TGAA analysis on the provided alignment file:

**pythonw tgaa.py -f Example_files/Alignment.faa -g Example_files/species_groups_file.txt -1 Atacama -2 Not_Atacama -o Alignment_tgaa_out.txt -a -p -m BLOSUM62 -H Kyte-Doolittle -t protein -gp=-11 -ag Murphy -sw 10 -s Example_files/species_labels_names.txt

Use the 'pythonw' when you want to plot ('-p') figures. Otherwise you can just use 'python'

### Dependencies:
- Python (https://www.python.org/downloads/)
- NumPy (https://numpy.org/). 'pip install numpy'.
- Matplotlib (https://matplotlib.org/). 'python -m pip install -U pip'; 'python -m pip install -U matplotlib'. Only if you want to plot figures (using the '-p' option)
- I only tested this program on Mac and Linux (Not sure if works on Windows - let me know...)

### Some notes before you start
1. Examples for input files can be found in the Example_files folder.
2. Will add a manual in the near future.
