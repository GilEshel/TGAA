# -*- coding: utf-8 -*-
#!/usr/bin/env python

#########################################################################################
### The MIT License (MIT)
### 
### Copyright (c) 2019 Gil Eshel
### 
### Permission is hereby granted, free of charge, to any person obtaining a copy
### of this software and associated documentation files (the "Software"), to deal
### in the Software without restriction, including without limitation the rights
### to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
### copies of the Software, and to permit persons to whom the Software is
### furnished to do so, subject to the following conditions:
### 
### The above copyright notice and this permission notice shall be included in all
### copies or substantial portions of the Software.
###
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
### AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
### OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
### SOFTWARE.
#########################################################################################
### DATE CREATED: Jan 23, 2019
### AUTHOR: Gil Eshel
### CONTACT1: ge30@nyu.edu
### CONTACT2: giltu1@gmail.com
###
### CITATION: To be added
### 
#########################################################################################
###	usage: python tgaa.py [-h] [-f FASTA_MSA_FILE]
###	                                             [-g SPECIES_GROUPS_FILE]
###	                                             [-1 GROUP1] [-2 GROUP2]
###	                                             [-t SEQUENCE_TYPE]
###	                                             [-m AA_SUB_MATRIX]
###	                                             [-H HYDROPHOBICITY_SCALE]
###	                                             [-gp GAP_PENALTY]
###	                                             [-ag AA_GROUPING]
###	                                             [-sw SLIDING_WINDOW_SIZE]
###	                                             [-o OUTPUT]
###	                                             [-s SPECIES_LABELS_NAMES_FILE]
###	                                             [-sd SPECIES_LABELS_NAMES_DELIMITER]
###	                                             [-p] [-a]
###	
###	TGAA is a program to parse amino acid multiple sequence
###	alignment (in a fasta format) , and identify aa positions that are conserved
###	in one group of species, but diverged in a different group of species (i.e.
###	positions with the same aa in one group, but different aa in the second group
###	- does not have to be the same aa in all species from second group)
###	
###	optional arguments:
###	  -h, --help            show this help message and exit
###	  -f FASTA_MSA_FILE, --fasta-msa-file FASTA_MSA_FILE
###	                        [required] Specify a multiple sequence alignment file,
###	                        in a fasta format
###	  -g SPECIES_GROUPS_FILE, --species-groups-file SPECIES_GROUPS_FILE
###	                        [required] A two-column (tab-delimited) file
###	                        containing a list of species labels (column 1), and
###	                        their group membership (column 2), first defined group
###	                        is considered the background (based on the first taxon
###	                        group membership). Expecting no column headers. e.g.
###	                        sp1 1 sp2 1 sp3  \sp4 2 . These species labels should
###	                        be added to the sequence identifiers, followed by '#',
###	                        e.g. Aratha#AT5G49450 for an Arabidopsis gene
###	  -1 GROUP1, --group1 GROUP1
###	                        [optional] Indicate a name for GROUP1. Default is
###	                        'GROUP1')
###	  -2 GROUP2, --group2 GROUP2
###	                        [optional] Indicate a name for GROUP2. Default is
###	                        'GROUP2')
###	  -t SEQUENCE_TYPE, --sequence-type SEQUENCE_TYPE
###	                        Indicate the sequence type, either 'protein' or
###	                        'codon'. Default is 'protein'
###	  -m AA_SUB_MATRIX, --aa-sub-matrix AA_SUB_MATRIX
###	                        Indicate the substitution matrix to use when scoring
###	                        (sum of pairs) the aa alignment. Choose among these
###	                        options: PAM1, PAM10, PAM30, PAM70, BLOSUM45,
###	                        BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90. Default:
###	                        BLOSUM62
###	  -H HYDROPHOBICITY_SCALE, --hydrophobicity-scale HYDROPHOBICITY_SCALE
###	                        Indicate the hydrophobicity scale that should be used
###	                        for plotting and mean calculations. Choose among these
###	                        options: Kyte-Doolittle, Hopp-Woods, Cornette,
###	                        Eisenberg , Rose, Janin, Engelman_GES. Default: Kyte-
###	                        Doolittle
###	  -gp GAP_PENALTY, --gap-penalty GAP_PENALTY
###	                        Indicate the gap penalty to use when calculating the
###	                        sum of pairs. Use '=', i.e.: -gp=-11. Default: -11
###	  -ag AA_GROUPING, --aa-grouping AA_GROUPING
###	                        Indicate the amino acid physiochemical classification
###	                        method desired. Select from the following: Taylor
###	                        (using Taylor 1986 classification), Murphy (using
###	                        Murphy 2000, 10 letter reduced alphabet), Katzir
###	                        (using Katzir 2006, 7 groups). Default: Taylor
###	  -sw SLIDING_WINDOW_SIZE, --sliding-window-size SLIDING_WINDOW_SIZE
###	                        For sliding window calculations, indicate the size of
###	                        the window (number of amino acid adjacent sites). Will
###	                        calculate the window mean and median hydrophobicity
###	                        and sum of pairs values, and plot them. Default: 10
###	  -o OUTPUT, --output OUTPUT
###	                        [optional] Specify the name of the output file.
###	                        Default: 'conserved_aa_diverging_groups.txt'
###	  -s SPECIES_LABELS_NAMES_FILE, --species-labels-names-file SPECIES_LABELS_NAMES_FILE
###	                        [optional] A two-column file containing the list of
###	                        species labels (column 1), and their species full name
###	                        (column 2). If provided, the full names will be used
###	                        to describe the species, instead of the species
###	                        labels. Expecting no column headers. e.g. sp1\\t
###	                        Arabidopsis thaliana\\nsp2\\tAnastatica hierochuntica\\n.
###	                        These species labels should be added to the sequence
###	                        identifiers, followed by '#', e.g. Aratha#AT5G49450
###	                        for an Arabidopsis gene
###	  -sd SPECIES_LABELS_NAMES_DELIMITER, --species-labels-names-delimiter SPECIES_LABELS_NAMES_DELIMITER
###	                        Indicate the delimiter of the species labels to full
###	                        names file [e.g ' ', ',', '\t'. Also 't', 'w', 'c' will
###	                        work]. Default is '\t')
###	  -p, --plot            Flag to plot sum of pairs and hydrophobicity plots. If
###	                        the --plot flag is indicated, various graphs will be
###	                        plotted. It requires that the python matplotlib
###	                        library will be perinstalled, and that this script
###	                        will be exacuted using 'pythonw
###	                        tgaa.py ' instead of 'python
###	                        tgaa.py '
###	  -a, --output-all-aa-grouping-file
###	                        Flag to output a file with all possible aa grouping
###	                        classifications, tp indicate if substitutions within
###	                        and between groups have a potential effect of
###	                        structure/function. It will output the following:
###	                        Taylor (using Taylor 1986 classification), Murphy
###	                        (using Murphy 2000, 10 letter reduced alphabet),
###	                        Katzir (using Katzir 2006, 7 groups), Grantham (using
###	                        Grantham 1974, physico-chemical distance matrix,
###	                        d>=100 indicate radical substitution), Zhang (using
###	                        Zhang 2000, charge classification, polarity
###	                        classification and polarity and volumn
###	                        classification), and Betts and Russell (2003), to
###	                        catch substitution of amino acids that doesn't
###	                        substitute particularly well with any other amino
###	                        acids.
###	
#########################################################################################

import sys, argparse, datetime
from itertools import islice	# for the sliding window function
def main():
	# Defining the arguments:
	parser = argparse.ArgumentParser(description="TGAA is a program to parse amino acid multiple sequence alignment (in a fasta format)\n, and identify aa positions that are conserved in one group of species, but diverged in a different \ngroup of species (i.e. positions with the same aa in one group, but different aa in the second group\n - does not have to be the same aa in all species from second group)", prog="python tgaa.py")
	parser.add_argument("-f","--fasta-msa-file", help="[required] Specify a multiple sequence alignment file, in a fasta format")
	parser.add_argument("-g","--species-groups-file", help="[required] A two-column (tab-delimited) file containing a list of species labels (column 1), and their group membership (column 2), first defined group is considered the background (based on the first taxon group membership). Expecting no column headers. e.g. sp1\t1\nsp2\t1\nsp3\t\2\n\sp4\t2\n. These species labels should be added to the sequence identifiers, followed by '#', e.g. Aratha label: Aratha#AT5G49450 for an Arabidopsis gene")
#	parser.add_argument("-g","--species-groups-file", help="[required] A two-line file containing a list of species labels for group1 and group2, respectively. e.g. GROUP1=sp1,sp2,sp3 and GROUP2=sp4,sp5,sp6. Species labels should be added to the sequence identifiers, followed by '#', e.g. Aratha#AT5G49450 for an Arabidopsis gene")
	parser.add_argument("-1","--group1",help="[optional] Indicate a name for GROUP1. Default is 'GROUP1')", default="GROUP1")
	parser.add_argument("-2","--group2",help="[optional] Indicate a name for GROUP2. Default is 'GROUP2')", default="GROUP2")
	parser.add_argument("-t","--sequence-type",help="Indicate the sequence type, either 'protein' or 'codon'. Default is 'protein'", default="protein")
	parser.add_argument("-m","--aa-sub-matrix",help="Indicate the substitution matrix to use when scoring (sum of pairs) the aa alignment. Choose among these options: PAM1, PAM10, PAM30, PAM70, BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90. Default: BLOSUM62", default="BLOSUM62")
	parser.add_argument("-H","--hydrophobicity-scale",help="Indicate the hydrophobicity scale that should be used for plotting and mean calculations. Choose among these options: Kyte-Doolittle, Hopp-Woods, Cornette, Eisenberg , Rose, Janin, Engelman_GES. Default: Kyte-Doolittle", default="Kyte-Doolittle")
	parser.add_argument("-gp","--gap-penalty", default=-11, help="Indicate the gap penalty to use when calculating the sum of pairs. Use '=', i.e.: -gp=-11. Default: -11") #, action='store_true', required=False, type=int,action="store_true", nargs='?'
	parser.add_argument("-ag","--aa-grouping", default="Taylor", help="Indicate the amino acid physiochemical classification method desired. Select from the following: Taylor (using Taylor 1986 classification), Murphy (using Murphy 2000, 10 letter reduced alphabet), Katzir (using Katzir 2006, 7 groups). Default: Taylor")
	parser.add_argument("-sw","--sliding-window-size", default=10, help="For sliding window calculations, indicate the size of the window (number of amino acid adjacent sites). Will calculate the window mean and median hydrophobicity and sum of pairs values, and plot them. Default: 10")
	parser.add_argument("-o","--output", default="conserved_aa_diverging_groups.txt", help="[optional] Specify the name of the output file. Default: 'conserved_aa_diverging_groups.txt'")
	parser.add_argument("-s","--species-labels-names-file", help="[optional] A two-column file containing the list of species labels (column 1), and their species full name (column 2). If provided, the full names will be used to describe the species, instead of the species labels. Expecting no column headers. e.g. sp1\\tArabidopsis thaliana\\nsp2\\tAnastatica hierochuntica\\n. These species labels should be added to the sequence identifiers, followed by '#', e.g. Aratha#AT5G49450 for an Arabidopsis gene")
	parser.add_argument("-sd","--species-labels-names-delimiter",help="Indicate the delimiter of the species labels to full names file [e.g ' ', ',', '\\t'. Also 't', 'w', 'c' will work]. Default is '\\t')", default="\t")
	parser.add_argument("-p","--plot",action='store_true',help="Flag to plot sum of pairs and hydrophobicity plots. If the --plot flag is indicated, various graphs will be plotted. It requires that the python matplotlib library will be perinstalled, and that this script will be exacuted using 'pythonw tgaa.py ' instead of 'python tgaa.py'")
	parser.add_argument("-a","--output-all-aa-grouping-file", action='store_true', help="Flag to output a file with all possible aa grouping classifications, tp indicate if substitutions within and between groups have a potential effect of structure/function. It will output the following: Taylor (using Taylor 1986 classification), Murphy (using Murphy 2000, 10 letter reduced alphabet), Katzir (using Katzir 2006, 7 groups), Grantham (using Grantham 1974, physico-chemical distance matrix, d>=100 indicate radical substitution), Zhang (using Zhang 2000, charge classification, polarity classification and polarity and volumn classification), and Betts and Russell (2003), to catch substitution of amino acids that doesn't substitute particularly well with any other amino acids.")
	if len(sys.argv[1:])==0:
		#parser.print_help()
		parser.print_usage() # for just the usage line
		parser.exit()
	args = parser.parse_args()
	#
	## Useful dictionaries:
	#
	codon = {
	'ATT':'I','ATC':'I','ATA':'I','ATH':'I',
	'CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L','YTR':'L','CTN':'L',
	'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GTN':'V',
	'TTT':'F','TTC':'F','TTY':'F',
	'ATG':'M',
	'TGT':'C','TGC':'C','UGY':'C',
	'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
	'GGT':'G','GGC':'G','GGA':'G','GGG':'G','GGN':'G',
	'CCT':'P','CCC':'P','CCA':'P','CCG':'P','CCN':'P',
	'ACT':'T','ACC':'T','ACA':'T','ACG':'T','ACN':'T',
	'TCT':'S','TCC':'S','TCA':'S','TCG':'S','TCN':'S','AGT':'S','AGC':'S','AGY':'S',
	'TAT':'Y','TAC':'Y','TAY':'Y',
	'TGG':'W',
	'CAA':'Q','CAG':'Q','CAR':'Q',
	'AAT':'N','AAC':'N','AAY':'N',
	'CAT':'H','CAC':'H','CAY':'H',
	'GAA':'E','GAG':'E','GAR':'E',
	'GAT':'D','GAC':'D','GAY':'D',
	'AAA':'K','AAG':'K','AAR':'K',
	'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R','CGN':'R','AGR':'R',
	'TAA':'*','TAG':'*','TGA':'*','TRA':'*','TAR':'*',
	'---':'-'
	}
	#'TAA':'Stop','TAG':'Stop','TGA':'Stop',
	#
	aa_names = {
	'I':'Isoleucine',
	'L':'Leucine',
	'V':'Valine',
	'F':'Phenylalanine',
	'M':'Methionine',
	'C':'Cysteine',
	'A':'Alanine',
	'G':'Glycine',
	'P':'Proline',
	'T':'Threonine',
	'S':'Serine',
	'Y':'Tyrosine',
	'W':'Tryptophan',
	'Q':'Glutamine',
	'N':'Asparagine',
	'H':'Histidine',
	'E':'Glutamate',
	'D':'Aspartate',
	'K':'Lysine',
	'R':'Arginine',
	'U':'Selenocysteine',	# unusual aa
	'O':'Pyrrolysine',	# unusual aa
	'B':'Asparagine or Aspartate',
	'Z':'Glutamine or Glutamate'
	}
	#
	#
	# Amino acid classification into 7 groups, based on the book "Encyclopedic Reference of Genomics and Proteomics in Molecular Medicine" from 2006. Chapter: "Amino Acids: Physicochemical Properties", page 59
	# Their description:
	# "A more specific classification of amino acids takes into
	# consideration the chemical nature of their side chains,
	# as outlined in Fig. 2. 
	# group 1 includes amino acids with aliphatic side-chains (Ile, Val, Leu, Ala, Gly), 
	# group 2, aromatic side-chains (Phe, Trp, Tyr), 
	# group 3, basic (positively-charged) side-chains (Lys, Arg, His), 
	# group 4, acidic (negatively-charged) side-chains or their corresponding amides (Asp, Glu, Asn, Gln) and 
	# group 5, aliphatic hydroxyl side chains (Ser, Thr). 
	# The sixth group consists of proline alone, since it has a secondary amine group, which affects the protein backbone conformation in a unique way and 
	# the seventh group consists of Cys and Met, which have sulfur-containing side chains."
	#
	aa_physchem_prop_grouping_ERGPMM_book = {
	'I':'Aliphatic side-chain',
	'L':'Aliphatic side-chain',
	'V':'Aliphatic side-chain',
	'F':'Aromatic side-chain',
	'M':'Sulfur-containing side chain',
	'C':'Sulfur-containing side chain',
	'A':'Aliphatic side-chain',
	'G':'Aliphatic side-chain',
	'P':'Proline_as_it_has_unique_effect',
	'T':'Aliphatic hydroxyl side chain',
	'S':'Aliphatic hydroxyl side chain',
	'Y':'Aromatic side-chain',
	'W':'Aromatic side-chain',
	'Q':'Acidic (negatively-charged) side-chain',
	'N':'Acidic (negatively-charged) side-chain',
	'H':'Basic (positively-charged) side-chain',
	'E':'Acidic (negatively-charged) side-chain',
	'D':'Acidic (negatively-charged) side-chain',
	'K':'Basic (positively-charged) side-chain',
	'R':'Basic (positively-charged) side-chain'
	}
	#
	prop_7_aa_grouping_ERGPMM_book = {
	'Aliphatic side-chain': ['I','V','L','A','G'], # Isoleucine (Ile), Valine (Val), Leucine (Leu), Alanine (Ala), Glycine (Gly)
	'Aromatic side-chain': ['F','W','Y'],	# Phenylalanine (Phe), Tryptophan (Trp), Tyrosine (Tyr)
	'Basic (positively-charged) side-chain': ['K','R','H'],	# Lysine (Lys), Arginine (Arg), Histidine (His)
	'Acidic (negatively-charged) side-chain': ['D','E','N','Q'],	# Aspartic acid (Asp), Glutamic acid (Glu), Asparagine (Asn), Glutamine (Gln)
	'Aliphatic hydroxyl side chain': ['S','T'],	# Serine (Ser), Threonine (Thr)
	'Proline_as_it_has_unique_effect':['P'],	# Proline (Pro)
	'Sulfur-containing side chain':['C','M']	# Cysteine (Cys), Methionine (Met)
	}
	#
	aa_physchem_prop_grouping_Murphy_2000 = { # Based on the Murphy simplified aa alphabet paper: https://doi.org/10.1093/protein/13.3.149
	'I':'LVIM',
	'L':'LVIM',
	'V':'LVIM',
	'F':'FWY',
	'M':'LVIM',
	'C':'C',
	'A':'A',
	'G':'G',
	'P':'P',
	'T':'ST',
	'S':'ST',
	'Y':'FWY',
	'W':'FWY',
	'Q':'EDNQ',
	'N':'EDNQ',
	'H':'H',
	'E':'EDNQ',
	'D':'EDNQ',
	'K':'KR',
	'R':'KR'
	}
	#
	prop_10_aa_grouping_Murphy_2000 = {	# Based on the Murphy simplified aa alphabet paper: https://doi.org/10.1093/protein/13.3.149
	'LVIM': ['L','V','I','M'],
	'C': ['C'],
	'A':['A'],
	'G':['G'],
	'ST':['S','T'],
	'P':['P'],
	'FWY':['F','W','Y'],
	'EDNO':['E','D','N','O'],
	'KR':['K','R'],
	'H':['H']
	}
	#
	# Taylor (1986) classification (mixture of classifications - complex...) - for each group define all the properties based on the amino acids (union), and compare with the other group to find how they differ
	aa_physchem_prop_Taylor = {
	'I':['Aliphatic','Hydrophobic'],
	'L':['Aliphatic','Hydrophobic'],
	'V':['Small','Aliphatic','Hydrophobic'],
	'F':['Aromatic','Hydrophobic'],
	'M':['Hydrophobic'],
	'C':['Tiny','Small','Hydrophobic','Polar'],	# for C, it depend on the context (C with S-H or C with S-S)
	'A':['Tiny','Small','Hydrophobic'],
	'G':['Tiny','Small','Hydrophobic'],
	'P':['Small'],
	'T':['Small','Hydrophobic','Polar'],
	'S':['Tiny','Small','Polar'],
	'Y':['Aromatic','Hydrophobic','Polar'],
	'W':['Aromatic','Hydrophobic','Polar'],
	'Q':['Polar'],
	'N':['Small','Polar'],
	'H':['Aromatic','Hydrophobic','Positive','Charged','Polar'],
	'E':['Negative','Charged','Polar'],
	'D':['Small','Negative','Charged','Polar'],
	'K':['Hydrophobic','Positive','Charged','Polar'],
	'R':['Positive','Charged','Polar']
	}
	#
	prop_aa_Taylor = {
	'Tiny':['G','A','C','S'],
	'Small': ['P','A','G','S','C','N','D','T','V'],
	'Aliphatic': ['I','V','L'],
	'Aromatic': ['F','Y','H','W'],
	'Hydrophobic': ['A','G','C','T','V','I','L','M','F','Y','W','H','K'],
	'Positive': ['H','K','R'],
	'Negative': ['D','E'],
	'Charged': ['D','E','H','K','R'],
	'Polar': ['Y','W','T','H','K','R','C','S','N','D','E','Q']
	}
	#
	# Grantham (1974) physico-chemical distance matrix (based on Dagan et al., 2002 - consider d<100 as a conservative amino acid replacement, and d>=100 as radical amino acid replacement)
	physico_chemical_distance_Grantham = {
	'I':{'I':0,'L':5,'V':29,'F':21,'M':10,'C':198,'A':94,'G':135,'P':95,'T':89,'S':142,'Y':33,'W':61,'Q':109,'N':149,'H':94,'E':134,'D':168,'K':102,'R':97},
	'L':{'I':5,'L':0,'V':32,'F':22,'M':15,'C':198,'A':96,'G':138,'P':98,'T':92,'S':145,'Y':36,'W':61,'Q':113,'N':153,'H':99,'E':138,'D':172,'K':107,'R':102},
	'V':{'I':29,'L':32,'V':0,'F':50,'M':21,'C':192,'A':64,'G':109,'P':68,'T':69,'S':124,'Y':55,'W':88,'Q':96,'N':133,'H':84,'E':121,'D':152,'K':97,'R':96},
	'F':{'I':21,'L':22,'V':50,'F':0,'M':28,'C':205,'A':113,'G':153,'P':114,'T':103,'S':155,'Y':22,'W':40,'Q':116,'N':158,'H':100,'E':140,'D':177,'K':102,'R':97},
	'M':{'I':10,'L':15,'V':21,'F':28,'M':0,'C':196,'A':84,'G':127,'P':87,'T':81,'S':135,'Y':36,'W':67,'Q':101,'N':142,'H':87,'E':126,'D':160,'K':95,'R':91},
	'C':{'I':198,'L':198,'V':192,'F':205,'M':196,'C':0,'A':195,'G':159,'P':169,'T':149,'S':112,'Y':194,'W':215,'Q':154,'N':139,'H':174,'E':170,'D':154,'K':202,'R':180},
	'A':{'I':94,'L':96,'V':64,'F':113,'M':84,'C':195,'A':0,'G':60,'P':27,'T':58,'S':99,'Y':112,'W':148,'Q':91,'N':111,'H':86,'E':107,'D':126,'K':106,'R':112},
	'G':{'I':135,'L':138,'V':109,'F':153,'M':127,'C':159,'A':60,'G':0,'P':42,'T':59,'S':56,'Y':147,'W':184,'Q':87,'N':80,'H':98,'E':98,'D':94,'K':127,'R':125},
	'P':{'I':95,'L':98,'V':68,'F':114,'M':87,'C':169,'A':27,'G':42,'P':0,'T':38,'S':74,'Y':110,'W':147,'Q':76,'N':91,'H':77,'E':93,'D':108,'K':103,'R':103},
	'T':{'I':89,'L':92,'V':69,'F':103,'M':81,'C':149,'A':58,'G':59,'P':38,'T':0,'S':58,'Y':92,'W':128,'Q':42,'N':65,'H':47,'E':65,'D':85,'K':78,'R':71},
	'S':{'I':142,'L':145,'V':124,'F':155,'M':135,'C':112,'A':99,'G':56,'P':74,'T':58,'S':0,'Y':144,'W':177,'Q':68,'N':46,'H':89,'E':80,'D':65,'K':121,'R':110},
	'Y':{'I':33,'L':36,'V':55,'F':22,'M':36,'C':194,'A':112,'G':147,'P':110,'T':92,'S':144,'Y':0,'W':37,'Q':99,'N':143,'H':83,'E':122,'D':160,'K':85,'R':77},
	'W':{'I':61,'L':61,'V':88,'F':40,'M':67,'C':215,'A':148,'G':184,'P':147,'T':128,'S':177,'Y':37,'W':0,'Q':130,'N':174,'H':115,'E':152,'D':181,'K':110,'R':101},
	'Q':{'I':109,'L':113,'V':96,'F':116,'M':101,'C':154,'A':91,'G':87,'P':76,'T':42,'S':68,'Y':99,'W':130,'Q':0,'N':46,'H':24,'E':29,'D':61,'K':53,'R':43},
	'N':{'I':149,'L':153,'V':133,'F':158,'M':142,'C':139,'A':111,'G':80,'P':91,'T':65,'S':46,'Y':143,'W':174,'Q':46,'N':0,'H':68,'E':42,'D':23,'K':94,'R':86},
	'H':{'I':94,'L':99,'V':84,'F':100,'M':87,'C':174,'A':86,'G':98,'P':77,'T':47,'S':89,'Y':83,'W':115,'Q':24,'N':68,'H':0,'E':40,'D':81,'K':32,'R':29},
	'E':{'I':134,'L':138,'V':121,'F':140,'M':126,'C':170,'A':107,'G':98,'P':93,'T':65,'S':80,'Y':122,'W':152,'Q':29,'N':42,'H':40,'E':0,'D':45,'K':56,'R':54},
	'D':{'I':168,'L':172,'V':152,'F':177,'M':160,'C':154,'A':126,'G':94,'P':108,'T':85,'S':65,'Y':160,'W':181,'Q':61,'N':23,'H':81,'E':45,'D':0,'K':101,'R':96},
	'K':{'I':102,'L':107,'V':97,'F':102,'M':95,'C':202,'A':106,'G':127,'P':103,'T':78,'S':121,'Y':85,'W':110,'Q':53,'N':94,'H':32,'E':56,'D':101,'K':0,'R':26},
	'R':{'I':97,'L':102,'V':96,'F':97,'M':91,'C':180,'A':112,'G':125,'P':103,'T':71,'S':110,'Y':77,'W':101,'Q':43,'N':86,'H':29,'E':54,'D':96,'K':26,'R':0}	
	}
	#
	# Zhang (2000) charge classification:
	charge_classification_Zhang = {
	'Positive': ['H','K','R'],
	'Negative': ['D','E'],
	'Neutral': ['A','N','C','Q','G','I','L','M','F','P','S','T','W','Y','V']
	}
	aa_charge_classification_Zhang = {
	'I':'Neutral',
	'L':'Neutral',
	'V':'Neutral',
	'F':'Neutral',
	'M':'Neutral',
	'C':'Neutral',
	'A':'Neutral',
	'G':'Neutral',
	'P':'Neutral',
	'T':'Neutral',
	'S':'Neutral',
	'Y':'Neutral',
	'W':'Neutral',
	'Q':'Neutral',
	'N':'Neutral',
	'H':'Positive',
	'E':'Negative',
	'D':'Negative',
	'K':'Positive',
	'R':'Positive'
	}
	polarity_classification_Zhang = {
	'Polar': ['R','N','D','C','Q','E','G','H','K','S','T','Y'],
	'Nonpolar': ['A','I','L','M','F','P','W','V']
	}
	aa_polarity_classification_Zhang = {
	'I':'Nonpolar',
	'L':'Nonpolar',
	'V':'Nonpolar',
	'F':'Nonpolar',
	'M':'Nonpolar',
	'C':'Polar',
	'A':'Nonpolar',
	'G':'Polar',
	'P':'Nonpolar',
	'T':'Polar',
	'S':'Polar',
	'Y':'Polar',
	'W':'Nonpolar',
	'Q':'Polar',
	'N':'Polar',
	'H':'Polar',
	'E':'Polar',
	'D':'Polar',
	'K':'Polar',
	'R':'Polar'
	}
	polarity_and_volumn_classification_Zhang = {
	'Special': ['C'],
	'Neutral_and_small': ['A','G','P','S','T'],
	'Polar_and_relatively_small': ['N','D','Q','E'],
	'Polar_and_relatively_large': ['H','K','R'],
	'Nonpolar_and_relatively_small': ['I','L','M','V'],
	'Nonpolar_and_relatively_large': ['F','W','Y']
	}
	aa_polarity_and_volumn_classification_Zhang = {
	'I':'Nonpolar_and_relatively_small',
	'L':'Nonpolar_and_relatively_small',
	'V':'Nonpolar_and_relatively_small',
	'F':'Nonpolar_and_relatively_large',
	'M':'Nonpolar_and_relatively_small',
	'C':'Special',
	'A':'Neutral_and_small',
	'G':'Neutral_and_small',
	'P':'Neutral_and_small',
	'T':'Neutral_and_small',
	'S':'Neutral_and_small',
	'Y':'Nonpolar_and_relatively_large',
	'W':'Nonpolar_and_relatively_large',
	'Q':'Polar_and_relatively_small',
	'N':'Polar_and_relatively_small',
	'H':'Polar_and_relatively_large',
	'E':'Polar_and_relatively_small',
	'D':'Polar_and_relatively_small',
	'K':'Polar_and_relatively_large',
	'R':'Polar_and_relatively_large'
	}
	#
	# BETTS and RUSSELL (2003), amino acids that does not substitute particularly well with any other amino acids:
	Uniq_aa_Betts_and_Russell = ['H','C','G','P','W']
	# 'H' is the most common amino acid in protein active or binding sites. It is commonly replaced by 'C' (which usually act together with 'H').
	# 'C' has a role that is very dependent on cellular location, so it is hard to predict the substitution effect - in extracellular proteins (outside of the cell), cysteines are frequently involved in disulphide bonds, where pairs of cysteines are oxidized to form a covalent bond (that stabilize the protein structure - so if one 'C' of the pair is lost, that have a major effect). Inside the cells or in the membrane (except within extracellular domains), the disulphide bonds are rare... In the intracellular proteins, 'C' can be important for metal binding (like 'H'), so in active sites. Or it has no role...
	# 'G' can be substitue for other small amino acids, but because it only has an hydrogen as its side chain, it can be placed in places where other aa can't be placed, and it allows more conformational flexibility. If a substitution in a conserved 'G' site occurs, it could have a drastic impact on function.
	# 'P' can sometimes be substitue for other small amino acids, but it is unique as it is the only amino acid where the side chain is connected to the protein backbone twice, forming a five-membered ring (imino acid). proline is unable to occupy many of the main-chain conformations easily adopted by all other amino acids (the opposite of glycine). It Often found in very tight turns in protein structures (i.e. where the polypeptide chain must change direction). It can also function to introduce kinks into α-helices, since it is unable to adopt a normal helical conformation. Despite being aliphatic the preference for turn structure means that prolines are usually found on the protein surface.
	# 'W' can be replaced by other aromatic residues, but it is unique in terms of chemistry and size, meaning that often replacement by anything could be disastrous.
	# 'R' is a positively-charged, polar amino acid. It thus most prefers to substitute for the other positively-charged amino acid, lysine, although in some circumstances it will also tolerate a change to other polar amino acids. Note that a change from arginine to lysine is not always neutral. In certain structural or functional contexts, such a mutation can be devastating to function. - only check if 'R' is substitute with 'K'
	# Molecular weights (adapted from http://biotools.nubic.northwestern.edu/proteincalc.html)	
	molecWeight = {
	'I':131.18,
	'L':131.18,
	'V':117.15,
	'F':165.19,
	'M':149.21,
	'C':121.16,
	'A':89.09,
	'G':75.07,
	'P':115.13,
	'T':119.12,
	'S':105.09,
	'Y':181.19,
	'W':204.23,
	'Q':146.15,
	'N':132.12,
	'H':155.16,
	'E':147.13,
	'D':133.10,
	'K':146.19,
	'R':174.20,
	'B':132.61,
	'Z':146.64
	}
	#
	#
	# Hydrophobicity scales, from the CLC genomics program: http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/602/index.php?manual=Hydrophobicity_scales.html
	# "Several hydrophobicity scales have been published for various uses. Many of the commonly used hydrophobicity scales are described below.
	# Kyte-Doolittle scale: The Kyte-Doolittle scale is widely used for detecting hydrophobic regions in proteins. Regions with a positive value are hydrophobic. This scale can be used for identifying both surface-exposed regions as well as transmembrane regions, depending on the window size used. Short window sizes of 5-7 generally work well for predicting putative surface-exposed regions. Large window sizes of 19-21 are well suited for finding transmembrane domains if the values calculated are above 1.6 [Kyte and Doolittle, 1982]. These values should be used as a rule of thumb and deviations from the rule may occur.
	# Engelman scale: The Engelman hydrophobicity scale, also known as the GES-scale, is another scale which can be used for prediction of protein hydrophobicity [Engelman et al., 1986]. As the Kyte-Doolittle scale, this scale is useful for predicting transmembrane regions in proteins.
	# Eisenberg scale: The Eisenberg scale is a normalized consensus hydrophobicity scale which shares many features with the other hydrophobicity scales [Eisenberg et al., 1984].
	# Hopp-Woods scale: Hopp and Woods developed their hydrophobicity scale for identification of potentially antigenic sites in proteins. This scale is basically a hydrophilic index where apolar residues have been assigned negative values. Antigenic sites are likely to be predicted when using a window size of 7 [Hopp and Woods, 1983].
	# Cornette scale: Cornette et al. computed an optimal hydrophobicity scale based on 28 published scales [Cornette et al., 1987]. This optimized scale is also suitable for prediction of alpha-helices in proteins.
	# Rose scale: The hydrophobicity scale by Rose et al. is correlated to the average area of buried amino acids in globular proteins [Rose et al., 1985]. This results in a scale which is not showing the helices of a protein, but rather the surface accessibility.
	# Janin scale: This scale also provides information about the accessible and buried amino acid residues of globular proteins [Janin, 1979].
	# Welling scale: Welling et al. used information on the relative occurrence of amino acids in antigenic regions to make a scale which is useful for prediction of antigenic regions. This method is better than the Hopp-Woods scale of hydrophobicity which is also used to identify antigenic regions.
	# Kolaskar-Tongaonkar: A semi-empirical method for prediction of antigenic regions has been developed [Kolaskar and Tongaonkar, 1990]. This method also includes information of surface accessibility and flexibility and at the time of publication the method was able to predict antigenic determinants with an accuracy of 75%.
	# Surface Probability: Display of surface probability based on the algorithm by [Emini et al., 1985]. This algorithm has been used to identify antigenic determinants on the surface of proteins.
	# Chain Flexibility: Display of backbone chain flexibility based on the algorithm by [Karplus and Schulz, 1985]. It is known that chain flexibility is an indication of a putative antigenic determinant.
	hydrophobicity_Kyte_Doolittle = {
	'I':4.5,
	'L':3.8,
	'V':4.2,
	'F':2.8,
	'M':1.9,
	'C':2.5,
	'A':1.8,
	'G':-0.4,
	'P':-1.6,
	'T':-0.7,
	'S':-0.8,
	'Y':-1.3,
	'W':-0.9,
	'Q':-3.5,
	'N':-3.5,
	'H':-3.2,
	'E':-3.5,
	'D':-3.5,
	'K':-3.9,
	'R':-4.5
	}
	#
	hydrophobicity_Hopp_Woods = {
	'I':-1.8,
	'L':-1.8,
	'V':-1.5,
	'F':-2.5,
	'M':-1.3,
	'C':-1,
	'A':-0.5,
	'G':0,
	'P':0,
	'T':-0.4,
	'S':0.3,
	'Y':-2.3,
	'W':-3.4,
	'Q':0.2,
	'N':0.2,
	'H':-0.5,
	'E':3,
	'D':3,
	'K':3,
	'R':3
	}
	#
	hydrophobicity_Cornette = {
	'I':4.8,
	'L':5.7,
	'V':4.7,
	'F':4.4,
	'M':4.2,
	'C':4.1,
	'A':0.2,
	'G':0,
	'P':-2.2,
	'T':-1.9,
	'S':-0.5,
	'Y':3.2,
	'W':1,
	'Q':-2.8,
	'N':-0.5,
	'H':0.5,
	'E':-1.8,
	'D':-3.1,
	'K':-3.1,
	'R':1.4
	}
	#
	hydrophobicity_Eisenberg = {
	'I':1.38,
	'L':1.06,
	'V':1.08,
	'F':1.19,
	'M':0.64,
	'C':0.29,
	'A':0.62,
	'G':0.48,
	'P':0.12,
	'T':-0.05,
	'S':-0.18,
	'Y':0.26,
	'W':0.81,
	'Q':-0.85,
	'N':-0.78,
	'H':-0.4,
	'E':-0.74,
	'D':-0.9,
	'K':-1.5,
	'R':-2.53
	}
	#
	hydrophobicity_Rose = {
	'I':0.88,
	'L':0.85,
	'V':0.86,
	'F':0.88,
	'M':0.85,
	'C':0.91,
	'A':0.74,
	'G':0.72,
	'P':0.64,
	'T':0.7,
	'S':0.66,
	'Y':0.76,
	'W':0.85,
	'Q':0.62,
	'N':0.63,
	'H':0.78,
	'E':0.62,
	'D':0.62,
	'K':0.52,
	'R':0.64
	}
	#
	hydrophobicity_Janin = {
	'I':0.7,
	'L':0.5,
	'V':0.6,
	'F':0.5,
	'M':0.4,
	'C':0.9,
	'A':0.3,
	'G':0.3,
	'P':-0.3,
	'T':-0.2,
	'S':-0.1,
	'Y':-0.4,
	'W':0.3,
	'Q':-0.7,
	'N':-0.5,
	'H':-0.1,
	'E':-0.7,
	'D':-0.6,
	'K':-1.8,
	'R':-1.4
	}
	#
	hydrophobicity_Engelman_GES = {
	'I':3.1,
	'L':2.8,
	'V':2.6,
	'F':3.7,
	'M':3.4,
	'C':2,
	'A':1.6,
	'G':1,
	'P':-0.2,
	'T':1.2,
	'S':0.6,
	'Y':-0.7,
	'W':1.9,
	'Q':-4.1,
	'N':-4.8,
	'H':-3,
	'E':-8.2,
	'D':-9.2,
	'K':-8.8,
	'R':-12.3
	}
	#
	# AA substitution matrices
	PAM1 = {	# for alignments with 99% identity
	'I':{'I':9872,'L':9,'V':33,'F':7,'M':12,'C':2,'A':2,'G':0,'P':0,'T':7,'S':1,'Y':1,'W':0,'Q':1,'N':3,'H':0,'E':2,'D':1,'K':2,'R':2},
	'L':{'I':22,'L':9947,'V':15,'F':13,'M':45,'C':0,'A':3,'G':1,'P':3,'T':3,'S':1,'Y':2,'W':4,'Q':6,'N':3,'H':4,'E':1,'D':0,'K':2,'R':1},
	'V':{'I':57,'L':11,'V':9901,'F':1,'M':17,'C':3,'A':13,'G':3,'P':3,'T':10,'S':2,'Y':2,'W':0,'Q':2,'N':1,'H':3,'E':2,'D':1,'K':1,'R':2},
	'F':{'I':8,'L':6,'V':0,'F':9946,'M':4,'C':0,'A':1,'G':1,'P':0,'T':1,'S':2,'Y':28,'W':3,'Q':0,'N':1,'H':2,'E':0,'D':0,'K':0,'R':1},
	'M':{'I':5,'L':8,'V':4,'F':1,'M':9874,'C':0,'A':1,'G':0,'P':0,'T':2,'S':1,'Y':0,'W':0,'Q':2,'N':0,'H':0,'E':0,'D':0,'K':4,'R':1},
	'C':{'I':1,'L':0,'V':2,'F':0,'M':0,'C':9973,'A':1,'G':0,'P':1,'T':1,'S':5,'Y':3,'W':0,'Q':0,'N':0,'H':1,'E':0,'D':0,'K':0,'R':1},
	'A':{'I':6,'L':4,'V':18,'F':2,'M':6,'C':3,'A':9867,'G':21,'P':22,'T':32,'S':35,'Y':2,'W':0,'Q':8,'N':9,'H':2,'E':17,'D':10,'K':2,'R':2},
	'G':{'I':0,'L':1,'V':5,'F':1,'M':1,'C':1,'A':21,'G':9935,'P':3,'T':3,'S':21,'Y':0,'W':0,'Q':3,'N':12,'H':1,'E':7,'D':11,'K':2,'R':1},
	'P':{'I':1,'L':2,'V':2,'F':1,'M':1,'C':1,'A':13,'G':2,'P':9926,'T':4,'S':12,'Y':0,'W':0,'Q':8,'N':2,'H':5,'E':3,'D':1,'K':2,'R':5},
	'T':{'I':11,'L':2,'V':9,'F':1,'M':6,'C':1,'A':22,'G':2,'P':5,'T':9871,'S':32,'Y':2,'W':0,'Q':3,'N':13,'H':1,'E':2,'D':4,'K':8,'R':2},
	'S':{'I':2,'L':1,'V':2,'F':3,'M':4,'C':11,'A':28,'G':16,'P':17,'T':38,'S':9840,'Y':2,'W':5,'Q':4,'N':34,'H':2,'E':6,'D':7,'K':7,'R':11},
	'Y':{'I':1,'L':1,'V':1,'F':21,'M':0,'C':3,'A':1,'G':0,'P':0,'T':1,'S':1,'Y':9945,'W':2,'Q':0,'N':3,'H':4,'E':1,'D':0,'K':0,'R':0},
	'W':{'I':0,'L':0,'V':0,'F':1,'M':0,'C':0,'A':0,'G':0,'P':0,'T':0,'S':1,'Y':1,'W':9976,'Q':0,'N':0,'H':0,'E':0,'D':0,'K':0,'R':2},
	'Q':{'I':1,'L':3,'V':1,'F':0,'M':4,'C':0,'A':3,'G':1,'P':6,'T':2,'S':2,'Y':0,'W':0,'Q':9876,'N':4,'H':23,'E':27,'D':5,'K':6,'R':9},
	'N':{'I':3,'L':1,'V':1,'F':1,'M':0,'C':0,'A':4,'G':6,'P':2,'T':9,'S':20,'Y':4,'W':1,'Q':4,'N':9822,'H':21,'E':6,'D':36,'K':13,'R':1},
	'H':{'I':0,'L':1,'V':1,'F':2,'M':0,'C':1,'A':1,'G':0,'P':3,'T':1,'S':1,'Y':4,'W':1,'Q':20,'N':18,'H':9912,'E':1,'D':3,'K':1,'R':8},
	'E':{'I':3,'L':1,'V':2,'F':0,'M':1,'C':0,'A':10,'G':4,'P':3,'T':2,'S':4,'Y':1,'W':0,'Q':35,'N':7,'H':2,'E':9865,'D':56,'K':4,'R':0},
	'D':{'I':1,'L':0,'V':1,'F':0,'M':0,'C':0,'A':6,'G':6,'P':1,'T':3,'S':5,'Y':0,'W':0,'Q':6,'N':42,'H':4,'E':53,'D':9859,'K':3,'R':0},
	'K':{'I':4,'L':1,'V':1,'F':0,'M':20,'C':0,'A':2,'G':2,'P':3,'T':11,'S':8,'Y':1,'W':0,'Q':12,'N':25,'H':2,'E':7,'D':6,'K':9926,'R':37},
	'R':{'I':3,'L':1,'V':1,'F':1,'M':4,'C':1,'A':1,'G':0,'P':4,'T':1,'S':6,'Y':0,'W':8,'Q':10,'N':1,'H':10,'E':0,'D':0,'K':19,'R':9913}
	}
	#
	PAM10 = {	# for alignments with 90% identity
	'I':{'I':9,'L':-4,'V':-1,'F':-5,'M':-3,'C':-9,'A':-8,'G':-17,'P':-12,'T':-5,'S':-10,'Y':-9,'W':-20,'Q':-11,'N':-8,'H':-13,'E':-8,'D':-11,'K':-9,'R':-8,'X':-8,'*':-23},
	'L':{'I':-4,'L':7,'V':-5,'F':-5,'M':-2,'C':-21,'A':-9,'G':-14,'P':-10,'T':-10,'S':-12,'Y':-10,'W':-9,'Q':-8,'N':-10,'H':-9,'E':-13,'D':-19,'K':-11,'R':-12,'X':-9,'*':-23},
	'V':{'I':-1,'L':-5,'V':8,'F':-12,'M':-4,'C':-9,'A':-5,'G':-9,'P':-9,'T':-6,'S':-10,'Y':-10,'W':-22,'Q':-10,'N':-12,'H':-9,'E':-10,'D':-11,'K':-13,'R':-11,'X':-8,'*':-23},
	'F':{'I':-5,'L':-5,'V':-12,'F':9,'M':-7,'C':-19,'A':-12,'G':-12,'P':-13,'T':-12,'S':-9,'Y':-1,'W':-7,'Q':-19,'N':-12,'H':-9,'E':-20,'D':-21,'K':-20,'R':-12,'X':-12,'*':-23},
	'M':{'I':-3,'L':-2,'V':-4,'F':-7,'M':12,'C':-20,'A':-8,'G':-12,'P':-11,'T':-7,'S':-8,'Y':-17,'W':-19,'Q':-7,'N':-15,'H':-17,'E':-10,'D':-17,'K':-4,'R':-7,'X':-9,'*':-23},
	'C':{'I':-9,'L':-21,'V':-9,'F':-19,'M':-20,'C':10,'A':-10,'G':-13,'P':-11,'T':-11,'S':-6,'Y':-7,'W':-22,'Q':-20,'N':-17,'H':-10,'E':-20,'D':-21,'K':-20,'R':-11,'X':-13,'*':-23},
	'A':{'I':-8,'L':-9,'V':-5,'F':-12,'M':-8,'C':-10,'A':7,'G':-4,'P':-4,'T':-3,'S':-3,'Y':-11,'W':-20,'Q':-7,'N':-7,'H':-11,'E':-5,'D':-6,'K':-10,'R':-10,'X':-6,'*':-23},
	'G':{'I':-17,'L':-14,'V':-9,'F':-12,'M':-12,'C':-13,'A':-4,'G':7,'P':-10,'T':-10,'S':-4,'Y':-20,'W':-21,'Q':-10,'N':-6,'H':-13,'E':-7,'D':-6,'K':-10,'R':-13,'X':-8,'*':-23},
	'P':{'I':-12,'L':-10,'V':-9,'F':-13,'M':-11,'C':-11,'A':-4,'G':-10,'P':8,'T':-7,'S':-4,'Y':-20,'W':-20,'Q':-6,'N':-9,'H':-7,'E':-9,'D':-12,'K':-10,'R':-7,'X':-8,'*':-23},
	'T':{'I':-5,'L':-10,'V':-6,'F':-12,'M':-7,'C':-11,'A':-3,'G':-10,'P':-7,'T':8,'S':-2,'Y':-9,'W':-19,'Q':-9,'N':-5,'H':-11,'E':-9,'D':-8,'K':-6,'R':-10,'X':-7,'*':-23},
	'S':{'I':-10,'L':-12,'V':-10,'F':-9,'M':-8,'C':-6,'A':-3,'G':-4,'P':-4,'T':-2,'S':7,'Y':-10,'W':-8,'Q':-8,'N':-2,'H':-9,'E':-7,'D':-7,'K':-7,'R':-6,'X':-6,'*':-23},
	'Y':{'I':-9,'L':-10,'V':-10,'F':-1,'M':-17,'C':-7,'A':-11,'G':-20,'P':-20,'T':-9,'S':-10,'Y':10,'W':-8,'Q':-18,'N':-7,'H':-6,'E':-11,'D':-17,'K':-12,'R':-14,'X':-11,'*':-23},
	'W':{'I':-20,'L':-9,'V':-22,'F':-7,'M':-19,'C':-22,'A':-20,'G':-21,'P':-20,'T':-19,'S':-8,'Y':-8,'W':13,'Q':-19,'N':-11,'H':-10,'E':-23,'D':-21,'K':-18,'R':-5,'X':-16,'*':-23},
	'Q':{'I':-11,'L':-8,'V':-10,'F':-19,'M':-7,'C':-20,'A':-7,'G':-10,'P':-6,'T':-9,'S':-8,'Y':-18,'W':-19,'Q':9,'N':-7,'H':-2,'E':-1,'D':-6,'K':-6,'R':-4,'X':-8,'*':-23},
	'N':{'I':-8,'L':-10,'V':-12,'F':-12,'M':-15,'C':-17,'A':-7,'G':-6,'P':-9,'T':-5,'S':-2,'Y':-7,'W':-11,'Q':-7,'N':9,'H':-2,'E':-5,'D':-1,'K':-4,'R':-9,'X':-6,'*':-23},
	'H':{'I':-13,'L':-9,'V':-9,'F':-9,'M':-17,'C':-10,'A':-11,'G':-13,'P':-7,'T':-11,'S':-9,'Y':-6,'W':-10,'Q':-2,'N':-2,'H':10,'E':-9,'D':-7,'K':-10,'R':-4,'X':-8,'*':-23},
	'E':{'I':-8,'L':-13,'V':-10,'F':-20,'M':-10,'C':-20,'A':-5,'G':-7,'P':-9,'T':-9,'S':-7,'Y':-11,'W':-23,'Q':-1,'N':-5,'H':-9,'E':8,'D':0,'K':-7,'R':-15,'X':-8,'*':-23},
	'D':{'I':-11,'L':-19,'V':-11,'F':-21,'M':-17,'C':-21,'A':-6,'G':-6,'P':-12,'T':-8,'S':-7,'Y':-17,'W':-21,'Q':-6,'N':-1,'H':-7,'E':0,'D':8,'K':-8,'R':-17,'X':-9,'*':-23},
	'K':{'I':-9,'L':-11,'V':-13,'F':-20,'M':-4,'C':-20,'A':-10,'G':-10,'P':-10,'T':-6,'S':-7,'Y':-12,'W':-18,'Q':-6,'N':-4,'H':-10,'E':-7,'D':-8,'K':7,'R':-2,'X':-8,'*':-23},
	'R':{'I':-8,'L':-12,'V':-11,'F':-12,'M':-7,'C':-11,'A':-10,'G':-13,'P':-7,'T':-10,'S':-6,'Y':-14,'W':-5,'Q':-4,'N':-9,'H':-4,'E':-15,'D':-17,'K':-2,'R':9,'X':-9,'*':-23},
	'X':{'I':-8,'L':-9,'V':-8,'F':-12,'M':-9,'C':-13,'A':-6,'G':-8,'P':-8,'T':-7,'S':-6,'Y':-11,'W':-16,'Q':-8,'N':-6,'H':-8,'E':-8,'D':-9,'K':-8,'R':-9,'X':-8,'*':-23,'X':-8,'*':-23},
	'*':{'I':-23,'L':-23,'V':-23,'F':-23,'M':-23,'C':-23,'A':-23,'G':-23,'P':-23,'T':-23,'S':-23,'Y':-23,'W':-23,'Q':-23,'N':-23,'H':-23,'E':-23,'D':-23,'K':-23,'R':-23,'X':-23,'*':1}
	}
	#
	PAM30 = {	# for alignments with 75% identity
	'I':{'I':8,'L':-1,'V':2,'F':-2,'M':-1,'C':-6,'A':-5,'G':-11,'P':-8,'T':-2,'S':-7,'Y':-6,'W':-14,'Q':-8,'N':-5,'H':-9,'E':-5,'D':-7,'K':-6,'R':-5,'X':-1,'*':-17},
	'L':{'I':-1,'L':7,'V':-2,'F':-3,'M':1,'C':-15,'A':-6,'G':-10,'P':-7,'T':-7,'S':-8,'Y':-7,'W':-6,'Q':-5,'N':-7,'H':-6,'E':-9,'D':-12,'K':-8,'R':-8,'X':-1,'*':-17},
	'V':{'I':2,'L':-2,'V':7,'F':-8,'M':-1,'C':-6,'A':-2,'G':-5,'P':-6,'T':-3,'S':-6,'Y':-7,'W':-15,'Q':-7,'N':-8,'H':-6,'E':-6,'D':-8,'K':-9,'R':-8,'X':-1,'*':-17},
	'F':{'I':-2,'L':-3,'V':-8,'F':9,'M':-4,'C':-13,'A':-8,'G':-9,'P':-10,'T':-9,'S':-6,'Y':2,'W':-4,'Q':-13,'N':-9,'H':-6,'E':-14,'D':-15,'K':-14,'R':-9,'X':-1,'*':-17},
	'M':{'I':-1,'L':1,'V':-1,'F':-4,'M':11,'C':-13,'A':-5,'G':-8,'P':-8,'T':-4,'S':-5,'Y':-11,'W':-13,'Q':-4,'N':-9,'H':-10,'E':-7,'D':-11,'K':-2,'R':-4,'X':-1,'*':-17},
	'C':{'I':-6,'L':-15,'V':-6,'F':-13,'M':-13,'C':10,'A':-6,'G':-9,'P':-8,'T':-8,'S':-3,'Y':-4,'W':-15,'Q':-14,'N':-11,'H':-7,'E':-14,'D':-14,'K':-14,'R':-8,'X':-1,'*':-17},
	'A':{'I':-5,'L':-6,'V':-2,'F':-8,'M':-5,'C':-6,'A':6,'G':-2,'P':-2,'T':-1,'S':0,'Y':-8,'W':-13,'Q':-4,'N':-4,'H':-7,'E':-2,'D':-3,'K':-7,'R':-7,'X':-1,'*':-17},
	'G':{'I':-11,'L':-10,'V':-5,'F':-9,'M':-8,'C':-9,'A':-2,'G':6,'P':-6,'T':-6,'S':-2,'Y':-14,'W':-15,'Q':-7,'N':-3,'H':-9,'E':-4,'D':-3,'K':-7,'R':-9,'X':-1,'*':-17},
	'P':{'I':-8,'L':-7,'V':-6,'F':-10,'M':-8,'C':-8,'A':-2,'G':-6,'P':8,'T':-4,'S':-2,'Y':-13,'W':-14,'Q':-3,'N':-6,'H':-4,'E':-5,'D':-8,'K':-6,'R':-4,'X':-1,'*':-17},
	'T':{'I':-2,'L':-7,'V':-3,'F':-9,'M':-4,'C':-8,'A':-1,'G':-6,'P':-4,'T':7,'S':0,'Y':-6,'W':-13,'Q':-5,'N':-2,'H':-7,'E':-6,'D':-5,'K':-3,'R':-6,'X':-1,'*':-17},
	'S':{'I':-7,'L':-8,'V':-6,'F':-6,'M':-5,'C':-3,'A':0,'G':-2,'P':-2,'T':0,'S':6,'Y':-7,'W':-5,'Q':-5,'N':0,'H':-6,'E':-4,'D':-4,'K':-4,'R':-3,'X':-1,'*':-17},
	'Y':{'I':-6,'L':-7,'V':-7,'F':2,'M':-11,'C':-4,'A':-8,'G':-14,'P':-13,'T':-6,'S':-7,'Y':10,'W':-5,'Q':-12,'N':-4,'H':-3,'E':-8,'D':-11,'K':-9,'R':-10,'X':-1,'*':-17},
	'W':{'I':-14,'L':-6,'V':-15,'F':-4,'M':-13,'C':-15,'A':-13,'G':-15,'P':-14,'T':-13,'S':-5,'Y':-5,'W':13,'Q':-13,'N':-8,'H':-7,'E':-17,'D':-15,'K':-12,'R':-2,'X':-1,'*':-17},
	'Q':{'I':-8,'L':-5,'V':-7,'F':-13,'M':-4,'C':-14,'A':-4,'G':-7,'P':-3,'T':-5,'S':-5,'Y':-12,'W':-13,'Q':8,'N':-3,'H':1,'E':1,'D':-2,'K':-3,'R':-2,'X':-1,'*':-17},
	'N':{'I':-5,'L':-7,'V':-8,'F':-9,'M':-9,'C':-11,'A':-4,'G':-3,'P':-6,'T':-2,'S':0,'Y':-4,'W':-8,'Q':-3,'N':8,'H':0,'E':-2,'D':2,'K':-1,'R':-6,'X':-1,'*':-17},
	'H':{'I':-9,'L':-6,'V':-6,'F':-6,'M':-10,'C':-7,'A':-7,'G':-9,'P':-4,'T':-7,'S':-6,'Y':-3,'W':-7,'Q':1,'N':0,'H':9,'E':-5,'D':-4,'K':-6,'R':-2,'X':-1,'*':-17},
	'E':{'I':-5,'L':-9,'V':-6,'F':-14,'M':-7,'C':-14,'A':-2,'G':-4,'P':-5,'T':-6,'S':-4,'Y':-8,'W':-17,'Q':1,'N':-2,'H':-5,'E':8,'D':2,'K':-4,'R':-9,'X':-1,'*':-17},
	'D':{'I':-7,'L':-12,'V':-8,'F':-15,'M':-11,'C':-14,'A':-3,'G':-3,'P':-8,'T':-5,'S':-4,'Y':-11,'W':-15,'Q':-2,'N':2,'H':-4,'E':2,'D':8,'K':-4,'R':-10,'X':-1,'*':-17},
	'K':{'I':-6,'L':-8,'V':-9,'F':-14,'M':-2,'C':-14,'A':-7,'G':-7,'P':-6,'T':-3,'S':-4,'Y':-9,'W':-12,'Q':-3,'N':-1,'H':-6,'E':-4,'D':-4,'K':7,'R':0,'X':-1,'*':-17},
	'R':{'I':-5,'L':-8,'V':-8,'F':-9,'M':-4,'C':-8,'A':-7,'G':-9,'P':-4,'T':-6,'S':-3,'Y':-10,'W':-2,'Q':-2,'N':-6,'H':-2,'E':-9,'D':-10,'K':0,'R':8,'X':-1,'*':-17},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-1,'A':-1,'G':-1,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-17},
	'*':{'I':-17,'L':-17,'V':-17,'F':-17,'M':-17,'C':-17,'A':-17,'G':-17,'P':-17,'T':-17,'S':-17,'Y':-17,'W':-17,'Q':-17,'N':-17,'H':-17,'E':-17,'D':-17,'K':-17,'R':-17,'X':-17,'*':1}
	}
	#
	PAM70 = {	# for alignments with 55% identity
	'I':{'I':7,'L':1,'V':3,'F':0,'M':1,'C':-4,'A':-2,'G':-6,'P':-5,'T':-1,'S':-4,'Y':-4,'W':-9,'Q':-5,'N':-3,'H':-6,'E':-4,'D':-5,'K':-4,'R':-3,'X':-1,'*':-11},
	'L':{'I':1,'L':6,'V':0,'F':-1,'M':2,'C':-10,'A':-4,'G':-7,'P':-5,'T':-4,'S':-6,'Y':-4,'W':-4,'Q':-3,'N':-5,'H':-4,'E':-6,'D':-8,'K':-5,'R':-6,'X':-1,'*':-11},
	'V':{'I':3,'L':0,'V':6,'F':-5,'M':0,'C':-4,'A':-1,'G':-3,'P':-3,'T':-1,'S':-3,'Y':-5,'W':-10,'Q':-4,'N':-5,'H':-4,'E':-4,'D':-5,'K':-6,'R':-5,'X':-1,'*':-11},
	'F':{'I':0,'L':-1,'V':-5,'F':8,'M':-2,'C':-8,'A':-6,'G':-7,'P':-7,'T':-6,'S':-4,'Y':4,'W':-2,'Q':-9,'N':-6,'H':-4,'E':-9,'D':-10,'K':-9,'R':-7,'X':-1,'*':-11},
	'M':{'I':1,'L':2,'V':0,'F':-2,'M':10,'C':-9,'A':-3,'G':-6,'P':-5,'T':-2,'S':-3,'Y':-7,'W':-8,'Q':-2,'N':-5,'H':-6,'E':-4,'D':-7,'K':0,'R':-2,'X':-1,'*':-11},
	'C':{'I':-4,'L':-10,'V':-4,'F':-8,'M':-9,'C':9,'A':-4,'G':-6,'P':-5,'T':-5,'S':-1,'Y':-2,'W':-11,'Q':-9,'N':-7,'H':-5,'E':-9,'D':-9,'K':-9,'R':-5,'X':-1,'*':-11},
	'A':{'I':-2,'L':-4,'V':-1,'F':-6,'M':-3,'C':-4,'A':5,'G':0,'P':0,'T':1,'S':1,'Y':-5,'W':-9,'Q':-2,'N':-2,'H':-4,'E':-1,'D':-1,'K':-4,'R':-4,'X':-1,'*':-11},
	'G':{'I':-6,'L':-7,'V':-3,'F':-7,'M':-6,'C':-6,'A':0,'G':6,'P':-3,'T':-3,'S':0,'Y':-9,'W':-10,'Q':-4,'N':-1,'H':-6,'E':-2,'D':-1,'K':-5,'R':-6,'X':-1,'*':-11},
	'P':{'I':-5,'L':-5,'V':-3,'F':-7,'M':-5,'C':-5,'A':0,'G':-3,'P':7,'T':-2,'S':0,'Y':-9,'W':-9,'Q':-1,'N':-3,'H':-2,'E':-3,'D':-4,'K':-4,'R':-2,'X':-1,'*':-11},
	'T':{'I':-1,'L':-4,'V':-1,'F':-6,'M':-2,'C':-5,'A':1,'G':-3,'P':-2,'T':6,'S':2,'Y':-4,'W':-8,'Q':-3,'N':0,'H':-4,'E':-3,'D':-2,'K':-1,'R':-4,'X':-1,'*':-11},
	'S':{'I':-4,'L':-6,'V':-3,'F':-4,'M':-3,'C':-1,'A':1,'G':0,'P':0,'T':2,'S':5,'Y':-5,'W':-3,'Q':-3,'N':1,'H':-3,'E':-2,'D':-1,'K':-2,'R':-1,'X':-1,'*':-11},
	'Y':{'I':-4,'L':-4,'V':-5,'F':4,'M':-7,'C':-2,'A':-5,'G':-9,'P':-9,'T':-4,'S':-5,'Y':9,'W':-3,'Q':-8,'N':-3,'H':-1,'E':-6,'D':-7,'K':-7,'R':-7,'X':-1,'*':-11},
	'W':{'I':-9,'L':-4,'V':-10,'F':-2,'M':-8,'C':-11,'A':-9,'G':-10,'P':-9,'T':-8,'S':-3,'Y':-3,'W':13,'Q':-8,'N':-6,'H':-5,'E':-11,'D':-10,'K':-7,'R':0,'X':-1,'*':-11},
	'Q':{'I':-5,'L':-3,'V':-4,'F':-9,'M':-2,'C':-9,'A':-2,'G':-4,'P':-1,'T':-3,'S':-3,'Y':-8,'W':-8,'Q':7,'N':-1,'H':2,'E':2,'D':0,'K':-1,'R':0,'X':-1,'*':-11},
	'N':{'I':-3,'L':-5,'V':-5,'F':-6,'M':-5,'C':-7,'A':-2,'G':-1,'P':-3,'T':0,'S':1,'Y':-3,'W':-6,'Q':-1,'N':6,'H':1,'E':0,'D':3,'K':0,'R':-3,'X':-1,'*':-11},
	'H':{'I':-6,'L':-4,'V':-4,'F':-4,'M':-6,'C':-5,'A':-4,'G':-6,'P':-2,'T':-4,'S':-3,'Y':-1,'W':-5,'Q':2,'N':1,'H':8,'E':-2,'D':-1,'K':-3,'R':0,'X':-1,'*':-11},
	'E':{'I':-4,'L':-6,'V':-4,'F':-9,'M':-4,'C':-9,'A':-1,'G':-2,'P':-3,'T':-3,'S':-2,'Y':-6,'W':-11,'Q':2,'N':0,'H':-2,'E':6,'D':3,'K':-2,'R':-5,'X':-1,'*':-11},
	'D':{'I':-5,'L':-8,'V':-5,'F':-10,'M':-7,'C':-9,'A':-1,'G':-1,'P':-4,'T':-2,'S':-1,'Y':-7,'W':-10,'Q':0,'N':3,'H':-1,'E':3,'D':6,'K':-2,'R':-6,'X':-1,'*':-11},
	'K':{'I':-4,'L':-5,'V':-6,'F':-9,'M':0,'C':-9,'A':-4,'G':-5,'P':-4,'T':-1,'S':-2,'Y':-7,'W':-7,'Q':-1,'N':0,'H':-3,'E':-2,'D':-2,'K':6,'R':2,'X':-1,'*':-11},
	'R':{'I':-3,'L':-6,'V':-5,'F':-7,'M':-2,'C':-5,'A':-4,'G':-6,'P':-2,'T':-4,'S':-1,'Y':-7,'W':0,'Q':0,'N':-3,'H':0,'E':-5,'D':-6,'K':2,'R':8,'X':-1,'*':-11},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-1,'A':-1,'G':-1,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-11},
	'*':{'I':-11,'L':-11,'V':-11,'F':-11,'M':-11,'C':-11,'A':-11,'G':-11,'P':-11,'T':-11,'S':-11,'Y':-11,'W':-11,'Q':-11,'N':-11,'H':-11,'E':-11,'D':-11,'K':-11,'R':-11,'X':-11,'*':1}
	}
	#
	BLOSUM45 = {	# for alignments with 45% identity
	'I':{'I':5,'L':2,'V':3,'F':0,'M':2,'C':-3,'A':-1,'G':-4,'P':-2,'T':-1,'S':-2,'Y':0,'W':-2,'Q':-2,'N':-2,'H':-3,'E':-3,'D':-4,'K':-3,'R':-3,'X':-1,'*':-5},
	'L':{'I':2,'L':5,'V':1,'F':1,'M':2,'C':-2,'A':-1,'G':-3,'P':-3,'T':-1,'S':-3,'Y':0,'W':-2,'Q':-2,'N':-3,'H':-2,'E':-2,'D':-3,'K':-3,'R':-2,'X':-1,'*':-5},
	'V':{'I':3,'L':1,'V':5,'F':0,'M':1,'C':-1,'A':0,'G':-3,'P':-3,'T':0,'S':-1,'Y':-1,'W':-3,'Q':-3,'N':-3,'H':-3,'E':-3,'D':-3,'K':-2,'R':-2,'X':-1,'*':-5},
	'F':{'I':0,'L':1,'V':0,'F':8,'M':0,'C':-2,'A':-2,'G':-3,'P':-3,'T':-1,'S':-2,'Y':3,'W':1,'Q':-4,'N':-2,'H':-2,'E':-3,'D':-4,'K':-3,'R':-2,'X':-1,'*':-5},
	'M':{'I':2,'L':2,'V':1,'F':0,'M':6,'C':-2,'A':-1,'G':-2,'P':-2,'T':-1,'S':-2,'Y':0,'W':-2,'Q':0,'N':-2,'H':0,'E':-2,'D':-3,'K':-1,'R':-1,'X':-1,'*':-5},
	'C':{'I':-3,'L':-2,'V':-1,'F':-2,'M':-2,'C':12,'A':-1,'G':-3,'P':-4,'T':-1,'S':-1,'Y':-3,'W':-5,'Q':-3,'N':-2,'H':-3,'E':-3,'D':-3,'K':-3,'R':-3,'X':-1,'*':-5},
	'A':{'I':-1,'L':-1,'V':0,'F':-2,'M':-1,'C':-1,'A':5,'G':0,'P':-1,'T':0,'S':1,'Y':-2,'W':-2,'Q':-1,'N':-1,'H':-2,'E':-1,'D':-2,'K':-1,'R':-2,'X':-1,'*':-5},
	'G':{'I':-4,'L':-3,'V':-3,'F':-3,'M':-2,'C':-3,'A':0,'G':7,'P':-2,'T':-2,'S':0,'Y':-3,'W':-2,'Q':-2,'N':0,'H':-2,'E':-2,'D':-1,'K':-2,'R':-2,'X':-1,'*':-5},
	'P':{'I':-2,'L':-3,'V':-3,'F':-3,'M':-2,'C':-4,'A':-1,'G':-2,'P':9,'T':-1,'S':-1,'Y':-3,'W':-3,'Q':-1,'N':-2,'H':-2,'E':0,'D':-1,'K':-1,'R':-2,'X':-1,'*':-5},
	'T':{'I':-1,'L':-1,'V':0,'F':-1,'M':-1,'C':-1,'A':0,'G':-2,'P':-1,'T':5,'S':2,'Y':-1,'W':-3,'Q':-1,'N':0,'H':-2,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-5},
	'S':{'I':-2,'L':-3,'V':-1,'F':-2,'M':-2,'C':-1,'A':1,'G':0,'P':-1,'T':2,'S':4,'Y':-2,'W':-4,'Q':0,'N':1,'H':-1,'E':0,'D':0,'K':-1,'R':-1,'X':-1,'*':-5},
	'Y':{'I':0,'L':0,'V':-1,'F':3,'M':0,'C':-3,'A':-2,'G':-3,'P':-3,'T':-1,'S':-2,'Y':8,'W':3,'Q':-1,'N':-2,'H':2,'E':-2,'D':-2,'K':-1,'R':-1,'X':-1,'*':-5},
	'W':{'I':-2,'L':-2,'V':-3,'F':1,'M':-2,'C':-5,'A':-2,'G':-2,'P':-3,'T':-3,'S':-4,'Y':3,'W':15,'Q':-2,'N':-4,'H':-3,'E':-3,'D':-4,'K':-2,'R':-2,'X':-1,'*':-5},
	'Q':{'I':-2,'L':-2,'V':-3,'F':-4,'M':0,'C':-3,'A':-1,'G':-2,'P':-1,'T':-1,'S':0,'Y':-1,'W':-2,'Q':6,'N':0,'H':1,'E':2,'D':0,'K':1,'R':1,'X':-1,'*':-5},
	'N':{'I':-2,'L':-3,'V':-3,'F':-2,'M':-2,'C':-2,'A':-1,'G':0,'P':-2,'T':0,'S':1,'Y':-2,'W':-4,'Q':0,'N':6,'H':1,'E':0,'D':2,'K':0,'R':0,'X':-1,'*':-5},
	'H':{'I':-3,'L':-2,'V':-3,'F':-2,'M':0,'C':-3,'A':-2,'G':-2,'P':-2,'T':-2,'S':-1,'Y':2,'W':-3,'Q':1,'N':1,'H':10,'E':0,'D':0,'K':-1,'R':0,'X':-1,'*':-5},
	'E':{'I':-3,'L':-2,'V':-3,'F':-3,'M':-2,'C':-3,'A':-1,'G':-2,'P':0,'T':-1,'S':0,'Y':-2,'W':-3,'Q':2,'N':0,'H':0,'E':6,'D':2,'K':1,'R':0,'X':-1,'*':-5},
	'D':{'I':-4,'L':-3,'V':-3,'F':-4,'M':-3,'C':-3,'A':-2,'G':-1,'P':-1,'T':-1,'S':0,'Y':-2,'W':-4,'Q':0,'N':2,'H':0,'E':2,'D':7,'K':0,'R':-1,'X':-1,'*':-5},
	'K':{'I':-3,'L':-3,'V':-2,'F':-3,'M':-1,'C':-3,'A':-1,'G':-2,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-2,'Q':1,'N':0,'H':-1,'E':1,'D':0,'K':5,'R':3,'X':-1,'*':-5},
	'R':{'I':-3,'L':-2,'V':-2,'F':-2,'M':-1,'C':-3,'A':-2,'G':-2,'P':-2,'T':-1,'S':-1,'Y':-1,'W':-2,'Q':1,'N':0,'H':0,'E':0,'D':-1,'K':3,'R':7,'X':-1,'*':-5},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-1,'A':-1,'G':-1,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-5},
	'*':{'I':-5,'L':-5,'V':-5,'F':-5,'M':-5,'C':-5,'A':-5,'G':-5,'P':-5,'T':-5,'S':-5,'Y':-5,'W':-5,'Q':-5,'N':-5,'H':-5,'E':-5,'D':-5,'K':-5,'R':-5,'X':-5,'*':1}
	}
	#
	BLOSUM50 = {	# for alignments with 50% identity
	'I':{'I':5,'L':2,'V':4,'F':0,'M':2,'C':-2,'A':-1,'G':-4,'P':-3,'T':-1,'S':-3,'Y':-1,'W':-3,'Q':-3,'N':-3,'H':-4,'E':-4,'D':-4,'K':-3,'R':-4,'X':-1,'*':-5},
	'L':{'I':2,'L':5,'V':1,'F':1,'M':3,'C':-2,'A':-2,'G':-4,'P':-4,'T':-1,'S':-3,'Y':-1,'W':-2,'Q':-2,'N':-4,'H':-3,'E':-3,'D':-4,'K':-3,'R':-3,'X':-1,'*':-5},
	'V':{'I':4,'L':1,'V':5,'F':-1,'M':1,'C':-1,'A':0,'G':-4,'P':-3,'T':0,'S':-2,'Y':-1,'W':-3,'Q':-3,'N':-3,'H':-4,'E':-3,'D':-4,'K':-3,'R':-3,'X':-1,'*':-5},
	'F':{'I':0,'L':1,'V':-1,'F':8,'M':0,'C':-2,'A':-3,'G':-4,'P':-4,'T':-2,'S':-3,'Y':4,'W':1,'Q':-4,'N':-4,'H':-1,'E':-3,'D':-5,'K':-4,'R':-3,'X':-1,'*':-5},
	'M':{'I':2,'L':3,'V':1,'F':0,'M':7,'C':-2,'A':-1,'G':-3,'P':-3,'T':-1,'S':-2,'Y':0,'W':-1,'Q':0,'N':-2,'H':-1,'E':-2,'D':-4,'K':-2,'R':-2,'X':-1,'*':-5},
	'C':{'I':-2,'L':-2,'V':-1,'F':-2,'M':-2,'C':13,'A':-1,'G':-3,'P':-4,'T':-1,'S':-1,'Y':-3,'W':-5,'Q':-3,'N':-2,'H':-3,'E':-3,'D':-4,'K':-3,'R':-4,'X':-1,'*':-5},
	'A':{'I':-1,'L':-2,'V':0,'F':-3,'M':-1,'C':-1,'A':5,'G':0,'P':-1,'T':0,'S':1,'Y':-2,'W':-3,'Q':-1,'N':-1,'H':-2,'E':-1,'D':-2,'K':-1,'R':-2,'X':-1,'*':-5},
	'G':{'I':-4,'L':-4,'V':-4,'F':-4,'M':-3,'C':-3,'A':0,'G':8,'P':-2,'T':-2,'S':0,'Y':-3,'W':-3,'Q':-2,'N':0,'H':-2,'E':-3,'D':-1,'K':-2,'R':-3,'X':-1,'*':-5},
	'P':{'I':-3,'L':-4,'V':-3,'F':-4,'M':-3,'C':-4,'A':-1,'G':-2,'P':10,'T':-1,'S':-1,'Y':-3,'W':-4,'Q':-1,'N':-2,'H':-2,'E':-1,'D':-1,'K':-1,'R':-3,'X':-1,'*':-5},
	'T':{'I':-1,'L':-1,'V':0,'F':-2,'M':-1,'C':-1,'A':0,'G':-2,'P':-1,'T':5,'S':2,'Y':-2,'W':-3,'Q':-1,'N':0,'H':-2,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-5},
	'S':{'I':-3,'L':-3,'V':-2,'F':-3,'M':-2,'C':-1,'A':1,'G':0,'P':-1,'T':2,'S':5,'Y':-2,'W':-4,'Q':0,'N':1,'H':-1,'E':-1,'D':0,'K':0,'R':-1,'X':-1,'*':-5},
	'Y':{'I':-1,'L':-1,'V':-1,'F':4,'M':0,'C':-3,'A':-2,'G':-3,'P':-3,'T':-2,'S':-2,'Y':8,'W':2,'Q':-1,'N':-2,'H':2,'E':-2,'D':-3,'K':-2,'R':-1,'X':-1,'*':-5},
	'W':{'I':-3,'L':-2,'V':-3,'F':1,'M':-1,'C':-5,'A':-3,'G':-3,'P':-4,'T':-3,'S':-4,'Y':2,'W':15,'Q':-1,'N':-4,'H':-3,'E':-3,'D':-5,'K':-3,'R':-3,'X':-1,'*':-5},
	'Q':{'I':-3,'L':-2,'V':-3,'F':-4,'M':0,'C':-3,'A':-1,'G':-2,'P':-1,'T':-1,'S':0,'Y':-1,'W':-1,'Q':7,'N':0,'H':1,'E':2,'D':0,'K':2,'R':1,'X':-1,'*':-5},
	'N':{'I':-3,'L':-4,'V':-3,'F':-4,'M':-2,'C':-2,'A':-1,'G':0,'P':-2,'T':0,'S':1,'Y':-2,'W':-4,'Q':0,'N':7,'H':1,'E':0,'D':2,'K':0,'R':-1,'X':-1,'*':-5},
	'H':{'I':-4,'L':-3,'V':-4,'F':-1,'M':-1,'C':-3,'A':-2,'G':-2,'P':-2,'T':-2,'S':-1,'Y':2,'W':-3,'Q':1,'N':1,'H':10,'E':0,'D':-1,'K':0,'R':0,'X':-1,'*':-5},
	'E':{'I':-4,'L':-3,'V':-3,'F':-3,'M':-2,'C':-3,'A':-1,'G':-3,'P':-1,'T':-1,'S':-1,'Y':-2,'W':-3,'Q':2,'N':0,'H':0,'E':6,'D':2,'K':1,'R':0,'X':-1,'*':-5},
	'D':{'I':-4,'L':-4,'V':-4,'F':-5,'M':-4,'C':-4,'A':-2,'G':-1,'P':-1,'T':-1,'S':0,'Y':-3,'W':-5,'Q':0,'N':2,'H':-1,'E':2,'D':8,'K':-1,'R':-2,'X':-1,'*':-5},
	'K':{'I':-3,'L':-3,'V':-3,'F':-4,'M':-2,'C':-3,'A':-1,'G':-2,'P':-1,'T':-1,'S':0,'Y':-2,'W':-3,'Q':2,'N':0,'H':0,'E':1,'D':-1,'K':6,'R':3,'X':-1,'*':-5},
	'R':{'I':-4,'L':-3,'V':-3,'F':-3,'M':-2,'C':-4,'A':-2,'G':-3,'P':-3,'T':-1,'S':-1,'Y':-1,'W':-3,'Q':1,'N':-1,'H':0,'E':0,'D':-2,'K':3,'R':7,'X':-1,'*':-5},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-1,'A':-1,'G':-1,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-5},
	'*':{'I':-5,'L':-5,'V':-5,'F':-5,'M':-5,'C':-5,'A':-5,'G':-5,'P':-5,'T':-5,'S':-5,'Y':-5,'W':-5,'Q':-5,'N':-5,'H':-5,'E':-5,'D':-5,'K':-5,'R':-5,'X':-5,'*':1}
	}
	#
	BLOSUM62 = {	# for alignments with 62% identity
	'I':{'I':4,'L':2,'V':3,'F':0,'M':1,'C':-1,'A':-1,'G':-4,'P':-3,'T':-1,'S':-2,'Y':-1,'W':-3,'Q':-3,'N':-3,'H':-3,'E':-3,'D':-3,'K':-3,'R':-3,'X':-1,'*':-4},
	'L':{'I':2,'L':4,'V':1,'F':0,'M':2,'C':-1,'A':-1,'G':-4,'P':-3,'T':-1,'S':-2,'Y':-1,'W':-2,'Q':-2,'N':-3,'H':-3,'E':-3,'D':-4,'K':-2,'R':-2,'X':-1,'*':-4},
	'V':{'I':3,'L':1,'V':4,'F':-1,'M':1,'C':-1,'A':0,'G':-3,'P':-2,'T':0,'S':-2,'Y':-1,'W':-3,'Q':-2,'N':-3,'H':-3,'E':-2,'D':-3,'K':-2,'R':-3,'X':-1,'*':-4},
	'F':{'I':0,'L':0,'V':-1,'F':6,'M':0,'C':-2,'A':-2,'G':-3,'P':-4,'T':-2,'S':-2,'Y':3,'W':1,'Q':-3,'N':-3,'H':-1,'E':-3,'D':-3,'K':-3,'R':-3,'X':-1,'*':-4},
	'M':{'I':1,'L':2,'V':1,'F':0,'M':5,'C':-1,'A':-1,'G':-3,'P':-2,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':0,'N':-2,'H':-2,'E':-2,'D':-3,'K':-1,'R':-1,'X':-1,'*':-4},
	'C':{'I':-1,'L':-1,'V':-1,'F':-2,'M':-1,'C':9,'A':0,'G':-3,'P':-3,'T':-1,'S':-1,'Y':-2,'W':-2,'Q':-3,'N':-3,'H':-3,'E':-4,'D':-3,'K':-3,'R':-3,'X':-2,'*':-4},
	'A':{'I':-1,'L':-1,'V':0,'F':-2,'M':-1,'C':0,'A':4,'G':0,'P':-1,'T':0,'S':1,'Y':-2,'W':-3,'Q':-1,'N':-2,'H':-2,'E':-1,'D':-2,'K':-1,'R':-1,'X':0,'*':-4},
	'G':{'I':-4,'L':-4,'V':-3,'F':-3,'M':-3,'C':-3,'A':0,'G':6,'P':-2,'T':-2,'S':0,'Y':-3,'W':-2,'Q':-2,'N':0,'H':-2,'E':-2,'D':-1,'K':-2,'R':-2,'X':-1,'*':-4},
	'P':{'I':-3,'L':-3,'V':-2,'F':-4,'M':-2,'C':-3,'A':-1,'G':-2,'P':7,'T':-1,'S':-1,'Y':-3,'W':-4,'Q':-1,'N':-2,'H':-2,'E':-1,'D':-1,'K':-1,'R':-2,'X':-2,'*':-4},
	'T':{'I':-1,'L':-1,'V':0,'F':-2,'M':-1,'C':-1,'A':0,'G':-2,'P':-1,'T':5,'S':1,'Y':-2,'W':-2,'Q':-1,'N':0,'H':-2,'E':-1,'D':-1,'K':-1,'R':-1,'X':0,'*':-4},
	'S':{'I':-2,'L':-2,'V':-2,'F':-2,'M':-1,'C':-1,'A':1,'G':0,'P':-1,'T':1,'S':4,'Y':-2,'W':-3,'Q':0,'N':1,'H':-1,'E':0,'D':0,'K':0,'R':-1,'X':0,'*':-4},
	'Y':{'I':-1,'L':-1,'V':-1,'F':3,'M':-1,'C':-2,'A':-2,'G':-3,'P':-3,'T':-2,'S':-2,'Y':7,'W':2,'Q':-1,'N':-2,'H':2,'E':-2,'D':-3,'K':-2,'R':-2,'X':-1,'*':-4},
	'W':{'I':-3,'L':-2,'V':-3,'F':1,'M':-1,'C':-2,'A':-3,'G':-2,'P':-4,'T':-2,'S':-3,'Y':2,'W':11,'Q':-2,'N':-4,'H':-2,'E':-3,'D':-4,'K':-3,'R':-3,'X':-2,'*':-4},
	'Q':{'I':-3,'L':-2,'V':-2,'F':-3,'M':0,'C':-3,'A':-1,'G':-2,'P':-1,'T':-1,'S':0,'Y':-1,'W':-2,'Q':5,'N':0,'H':0,'E':2,'D':0,'K':1,'R':1,'X':-1,'*':-4},
	'N':{'I':-3,'L':-3,'V':-3,'F':-3,'M':-2,'C':-3,'A':-2,'G':0,'P':-2,'T':0,'S':1,'Y':-2,'W':-4,'Q':0,'N':6,'H':1,'E':0,'D':1,'K':0,'R':0,'X':-1,'*':-4},
	'H':{'I':-3,'L':-3,'V':-3,'F':-1,'M':-2,'C':-3,'A':-2,'G':-2,'P':-2,'T':-2,'S':-1,'Y':2,'W':-2,'Q':0,'N':1,'H':8,'E':0,'D':-1,'K':-1,'R':0,'X':-1,'*':-4},
	'E':{'I':-3,'L':-3,'V':-2,'F':-3,'M':-2,'C':-4,'A':-1,'G':-2,'P':-1,'T':-1,'S':0,'Y':-2,'W':-3,'Q':2,'N':0,'H':0,'E':5,'D':2,'K':1,'R':0,'X':-1,'*':-4},
	'D':{'I':-3,'L':-4,'V':-3,'F':-3,'M':-3,'C':-3,'A':-2,'G':-1,'P':-1,'T':-1,'S':0,'Y':-3,'W':-4,'Q':0,'N':1,'H':-1,'E':2,'D':6,'K':-1,'R':-2,'X':-1,'*':-4},
	'K':{'I':-3,'L':-2,'V':-2,'F':-3,'M':-1,'C':-3,'A':-1,'G':-2,'P':-1,'T':-1,'S':0,'Y':-2,'W':-3,'Q':1,'N':0,'H':-1,'E':1,'D':-1,'K':5,'R':2,'X':-1,'*':-4},
	'R':{'I':-3,'L':-2,'V':-3,'F':-3,'M':-1,'C':-3,'A':-1,'G':-2,'P':-2,'T':-1,'S':-1,'Y':-2,'W':-3,'Q':1,'N':0,'H':0,'E':0,'D':-2,'K':2,'R':5,'X':-1,'*':-4},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-2,'A':0,'G':-1,'P':-2,'T':0,'S':0,'Y':-1,'W':-2,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-4},
	'*':{'I':-4,'L':-4,'V':-4,'F':-4,'M':-4,'C':-4,'A':-4,'G':-4,'P':-4,'T':-4,'S':-4,'Y':-4,'W':-4,'Q':-4,'N':-4,'H':-4,'E':-4,'D':-4,'K':-4,'R':-4,'X':-4,'*':1}
	}
	#
	BLOSUM80 = {	# for alignments with 80% identity
	'I':{'I':5,'L':1,'V':3,'F':-1,'M':1,'C':-2,'A':-2,'G':-5,'P':-4,'T':-1,'S':-3,'Y':-2,'W':-3,'Q':-3,'N':-4,'H':-4,'E':-4,'D':-4,'K':-3,'R':-3,'X':-1,'*':-6},
	'L':{'I':1,'L':4,'V':1,'F':0,'M':2,'C':-2,'A':-2,'G':-4,'P':-3,'T':-2,'S':-3,'Y':-2,'W':-2,'Q':-3,'N':-4,'H':-3,'E':-4,'D':-5,'K':-3,'R':-3,'X':-1,'*':-6},
	'V':{'I':3,'L':1,'V':4,'F':-1,'M':1,'C':-1,'A':0,'G':-4,'P':-3,'T':0,'S':-2,'Y':-2,'W':-3,'Q':-3,'N':-4,'H':-4,'E':-3,'D':-4,'K':-3,'R':-3,'X':-1,'*':-6},
	'F':{'I':-1,'L':0,'V':-1,'F':6,'M':0,'C':-3,'A':-3,'G':-4,'P':-4,'T':-2,'S':-3,'Y':3,'W':0,'Q':-4,'N':-4,'H':-2,'E':-4,'D':-4,'K':-4,'R':-4,'X':-1,'*':-6},
	'M':{'I':1,'L':2,'V':1,'F':0,'M':6,'C':-2,'A':-1,'G':-4,'P':-3,'T':-1,'S':-2,'Y':-2,'W':-2,'Q':0,'N':-3,'H':-2,'E':-2,'D':-4,'K':-2,'R':-2,'X':-1,'*':-6},
	'C':{'I':-2,'L':-2,'V':-1,'F':-3,'M':-2,'C':9,'A':-1,'G':-4,'P':-4,'T':-1,'S':-2,'Y':-3,'W':-3,'Q':-4,'N':-3,'H':-4,'E':-5,'D':-4,'K':-4,'R':-4,'X':-1,'*':-6},
	'A':{'I':-2,'L':-2,'V':0,'F':-3,'M':-1,'C':-1,'A':5,'G':0,'P':-1,'T':0,'S':1,'Y':-2,'W':-3,'Q':-1,'N':-2,'H':-2,'E':-1,'D':-2,'K':-1,'R':-2,'X':-1,'*':-6},
	'G':{'I':-5,'L':-4,'V':-4,'F':-4,'M':-4,'C':-4,'A':0,'G':6,'P':-3,'T':-2,'S':-1,'Y':-4,'W':-4,'Q':-2,'N':-1,'H':-3,'E':-3,'D':-2,'K':-2,'R':-3,'X':-1,'*':-6},
	'P':{'I':-4,'L':-3,'V':-3,'F':-4,'M':-3,'C':-4,'A':-1,'G':-3,'P':8,'T':-2,'S':-1,'Y':-4,'W':-5,'Q':-2,'N':-3,'H':-3,'E':-2,'D':-2,'K':-1,'R':-2,'X':-1,'*':-6},
	'T':{'I':-1,'L':-2,'V':0,'F':-2,'M':-1,'C':-1,'A':0,'G':-2,'P':-2,'T':5,'S':1,'Y':-2,'W':-4,'Q':-1,'N':0,'H':-2,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-6},
	'S':{'I':-3,'L':-3,'V':-2,'F':-3,'M':-2,'C':-2,'A':1,'G':-1,'P':-1,'T':1,'S':5,'Y':-2,'W':-4,'Q':0,'N':0,'H':-1,'E':0,'D':-1,'K':-1,'R':-1,'X':-1,'*':-6},
	'Y':{'I':-2,'L':-2,'V':-2,'F':3,'M':-2,'C':-3,'A':-2,'G':-4,'P':-4,'T':-2,'S':-2,'Y':7,'W':2,'Q':-2,'N':-3,'H':2,'E':-3,'D':-4,'K':-3,'R':-3,'X':-1,'*':-6},
	'W':{'I':-3,'L':-2,'V':-3,'F':0,'M':-2,'C':-3,'A':-3,'G':-4,'P':-5,'T':-4,'S':-4,'Y':2,'W':11,'Q':-3,'N':-4,'H':-3,'E':-4,'D':-6,'K':-4,'R':-4,'X':-1,'*':-6},
	'Q':{'I':-3,'L':-3,'V':-3,'F':-4,'M':0,'C':-4,'A':-1,'G':-2,'P':-2,'T':-1,'S':0,'Y':-2,'W':-3,'Q':6,'N':0,'H':1,'E':2,'D':-1,'K':1,'R':1,'X':-1,'*':-6},
	'N':{'I':-4,'L':-4,'V':-4,'F':-4,'M':-3,'C':-3,'A':-2,'G':-1,'P':-3,'T':0,'S':0,'Y':-3,'W':-4,'Q':0,'N':6,'H':0,'E':-1,'D':1,'K':0,'R':-1,'X':-1,'*':-6},
	'H':{'I':-4,'L':-3,'V':-4,'F':-2,'M':-2,'C':-4,'A':-2,'G':-3,'P':-3,'T':-2,'S':-1,'Y':2,'W':-3,'Q':1,'N':0,'H':8,'E':0,'D':-2,'K':-1,'R':0,'X':-1,'*':-6},
	'E':{'I':-4,'L':-4,'V':-3,'F':-4,'M':-2,'C':-5,'A':-1,'G':-3,'P':-2,'T':-1,'S':0,'Y':-3,'W':-4,'Q':2,'N':-1,'H':0,'E':6,'D':1,'K':1,'R':-1,'X':-1,'*':-6},
	'D':{'I':-4,'L':-5,'V':-4,'F':-4,'M':-4,'C':-4,'A':-2,'G':-2,'P':-2,'T':-1,'S':-1,'Y':-4,'W':-6,'Q':-1,'N':1,'H':-2,'E':1,'D':6,'K':-1,'R':-2,'X':-1,'*':-6},
	'K':{'I':-3,'L':-3,'V':-3,'F':-4,'M':-2,'C':-4,'A':-1,'G':-2,'P':-1,'T':-1,'S':-1,'Y':-3,'W':-4,'Q':1,'N':0,'H':-1,'E':1,'D':-1,'K':5,'R':2,'X':-1,'*':-6},
	'R':{'I':-3,'L':-3,'V':-3,'F':-4,'M':-2,'C':-4,'A':-2,'G':-3,'P':-2,'T':-1,'S':-1,'Y':-3,'W':-4,'Q':1,'N':-1,'H':0,'E':-1,'D':-2,'K':2,'R':6,'X':-1,'*':-6},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-1,'A':-1,'G':-1,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-6},
	'*':{'I':-6,'L':-6,'V':-6,'F':-6,'M':-6,'C':-6,'A':-6,'G':-6,'P':-6,'T':-6,'S':-6,'Y':-6,'W':-6,'Q':-6,'N':-6,'H':-6,'E':-6,'D':-6,'K':-6,'R':-6,'X':-6,'*':1}
	}
	#
	BLOSUM90 = {	# for alignments with 90% identity
	'I':{'I':5,'L':1,'V':3,'F':-1,'M':1,'C':-2,'A':-2,'G':-5,'P':-4,'T':-1,'S':-3,'Y':-2,'W':-4,'Q':-4,'N':-4,'H':-4,'E':-4,'D':-5,'K':-4,'R':-4,'X':-1,'*':-6},
	'L':{'I':1,'L':5,'V':0,'F':0,'M':2,'C':-2,'A':-2,'G':-5,'P':-4,'T':-2,'S':-3,'Y':-2,'W':-3,'Q':-3,'N':-4,'H':-4,'E':-4,'D':-5,'K':-3,'R':-3,'X':-1,'*':-6},
	'V':{'I':3,'L':0,'V':5,'F':-2,'M':0,'C':-2,'A':-1,'G':-5,'P':-3,'T':-1,'S':-2,'Y':-3,'W':-3,'Q':-3,'N':-4,'H':-4,'E':-3,'D':-5,'K':-3,'R':-3,'X':-1,'*':-6},
	'F':{'I':-1,'L':0,'V':-2,'F':7,'M':-1,'C':-3,'A':-3,'G':-5,'P':-4,'T':-3,'S':-3,'Y':3,'W':0,'Q':-4,'N':-4,'H':-2,'E':-5,'D':-5,'K':-4,'R':-4,'X':-1,'*':-6},
	'M':{'I':1,'L':2,'V':0,'F':-1,'M':7,'C':-2,'A':-2,'G':-4,'P':-3,'T':-1,'S':-2,'Y':-2,'W':-2,'Q':0,'N':-3,'H':-3,'E':-3,'D':-4,'K':-2,'R':-2,'X':-1,'*':-6},
	'C':{'I':-2,'L':-2,'V':-2,'F':-3,'M':-2,'C':9,'A':-1,'G':-4,'P':-4,'T':-2,'S':-2,'Y':-4,'W':-4,'Q':-4,'N':-4,'H':-5,'E':-6,'D':-5,'K':-4,'R':-5,'X':-1,'*':-6},
	'A':{'I':-2,'L':-2,'V':-1,'F':-3,'M':-2,'C':-1,'A':5,'G':0,'P':-1,'T':0,'S':1,'Y':-3,'W':-4,'Q':-1,'N':-2,'H':-2,'E':-1,'D':-3,'K':-1,'R':-2,'X':-1,'*':-6},
	'G':{'I':-5,'L':-5,'V':-5,'F':-5,'M':-4,'C':-4,'A':0,'G':6,'P':-3,'T':-3,'S':-1,'Y':-5,'W':-4,'Q':-3,'N':-1,'H':-3,'E':-3,'D':-2,'K':-2,'R':-3,'X':-1,'*':-6},
	'P':{'I':-4,'L':-4,'V':-3,'F':-4,'M':-3,'C':-4,'A':-1,'G':-3,'P':8,'T':-2,'S':-2,'Y':-4,'W':-5,'Q':-2,'N':-3,'H':-3,'E':-2,'D':-3,'K':-2,'R':-3,'X':-1,'*':-6},
	'T':{'I':-1,'L':-2,'V':-1,'F':-3,'M':-1,'C':-2,'A':0,'G':-3,'P':-2,'T':6,'S':1,'Y':-2,'W':-4,'Q':-1,'N':0,'H':-2,'E':-1,'D':-2,'K':-1,'R':-2,'X':-1,'*':-6},
	'S':{'I':-3,'L':-3,'V':-2,'F':-3,'M':-2,'C':-2,'A':1,'G':-1,'P':-2,'T':1,'S':5,'Y':-3,'W':-4,'Q':-1,'N':0,'H':-2,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-6},
	'Y':{'I':-2,'L':-2,'V':-3,'F':3,'M':-2,'C':-4,'A':-3,'G':-5,'P':-4,'T':-2,'S':-3,'Y':8,'W':2,'Q':-3,'N':-3,'H':1,'E':-4,'D':-4,'K':-3,'R':-3,'X':-1,'*':-6},
	'W':{'I':-4,'L':-3,'V':-3,'F':0,'M':-2,'C':-4,'A':-4,'G':-4,'P':-5,'T':-4,'S':-4,'Y':2,'W':11,'Q':-3,'N':-5,'H':-3,'E':-5,'D':-6,'K':-5,'R':-4,'X':-1,'*':-6},
	'Q':{'I':-4,'L':-3,'V':-3,'F':-4,'M':0,'C':-4,'A':-1,'G':-3,'P':-2,'T':-1,'S':-1,'Y':-3,'W':-3,'Q':7,'N':0,'H':1,'E':2,'D':-1,'K':1,'R':1,'X':-1,'*':-6},
	'N':{'I':-4,'L':-4,'V':-4,'F':-4,'M':-3,'C':-4,'A':-2,'G':-1,'P':-3,'T':0,'S':0,'Y':-3,'W':-5,'Q':0,'N':7,'H':0,'E':-1,'D':1,'K':0,'R':-1,'X':-1,'*':-6},
	'H':{'I':-4,'L':-4,'V':-4,'F':-2,'M':-3,'C':-5,'A':-2,'G':-3,'P':-3,'T':-2,'S':-2,'Y':1,'W':-3,'Q':1,'N':0,'H':8,'E':-1,'D':-2,'K':-1,'R':0,'X':-1,'*':-6},
	'E':{'I':-4,'L':-4,'V':-3,'F':-5,'M':-3,'C':-6,'A':-1,'G':-3,'P':-2,'T':-1,'S':-1,'Y':-4,'W':-5,'Q':2,'N':-1,'H':-1,'E':6,'D':1,'K':0,'R':-1,'X':-1,'*':-6},
	'D':{'I':-5,'L':-5,'V':-5,'F':-5,'M':-4,'C':-5,'A':-3,'G':-2,'P':-3,'T':-2,'S':-1,'Y':-4,'W':-6,'Q':-1,'N':1,'H':-2,'E':1,'D':7,'K':-1,'R':-3,'X':-1,'*':-6},
	'K':{'I':-4,'L':-3,'V':-3,'F':-4,'M':-2,'C':-4,'A':-1,'G':-2,'P':-2,'T':-1,'S':-1,'Y':-3,'W':-5,'Q':1,'N':0,'H':-1,'E':0,'D':-1,'K':6,'R':2,'X':-1,'*':-6},
	'R':{'I':-4,'L':-3,'V':-3,'F':-4,'M':-2,'C':-5,'A':-2,'G':-3,'P':-3,'T':-2,'S':-1,'Y':-3,'W':-4,'Q':1,'N':-1,'H':0,'E':-1,'D':-3,'K':2,'R':6,'X':-1,'*':-6},
	'X':{'I':-1,'L':-1,'V':-1,'F':-1,'M':-1,'C':-1,'A':-1,'G':-1,'P':-1,'T':-1,'S':-1,'Y':-1,'W':-1,'Q':-1,'N':-1,'H':-1,'E':-1,'D':-1,'K':-1,'R':-1,'X':-1,'*':-6},
	'*':{'I':-6,'L':-6,'V':-6,'F':-6,'M':-6,'C':-6,'A':-6,'G':-6,'P':-6,'T':-6,'S':-6,'Y':-6,'W':-6,'Q':-6,'N':-6,'H':-6,'E':-6,'D':-6,'K':-6,'R':-6,'X':-6,'*':1}
	}
	#
	## Function to calculate the sum of pairs - per aa position
	def sum_of_pairs(aa_list,sub_matrix,gap_pen): # Takes a list of aa characters for a particular alignment position, the substitution matrix and the gap penalty to use, and return the sum of pairs score
		pairwise_scores = []
		for i in range(len(aa_list)):
			for j in range(i + 1, len(aa_list)):
				#if aa_list[i] == '-' and aa_list[j] == '-':	# if both characters are gaps --> ignore...
				#	continue
				if aa_list[i] == '-' or aa_list[j] == '-':	# if one of them is a gap, penalize based on the defined gap penalty
				#elif aa_list[i] == '-' or aa_list[j] == '-':	# if one of them is a gap, penalize based on the defined gap penalty
					pairwise_scores.append(int(gap_pen))
				else:
					# Treat Selenocysteine (U) as Cysteine (C), and Pyrrolysine (O) as Any amino acid (X), as https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/multalignviewer/multalignviewer.html
					aa_i = aa_list[i]
					if aa_i == 'U':
						aa_i = 'C'
					elif aa_i == 'O':
						aa_i = 'X'
					aa_j = aa_list[j]
					if aa_j == 'U':
						aa_j = 'C'
					elif aa_j == 'O':
						aa_j = 'X'
					if sub_matrix == 'PAM1' or sub_matrix == 'PAM 1' or sub_matrix == 'pam1' or sub_matrix == 'pam 1' or sub_matrix == 'Pam1' or sub_matrix == 'Pam 1':
						pairwise_scores.append(PAM1[aa_i][aa_i])
					elif sub_matrix == 'PAM10' or sub_matrix == 'PAM 10' or sub_matrix == 'pam10' or sub_matrix == 'pam 10' or sub_matrix == 'Pam10' or sub_matrix == 'Pam 10':
						pairwise_scores.append(PAM10[aa_i][aa_j])
					elif sub_matrix == 'PAM30' or sub_matrix == 'PAM 30' or sub_matrix == 'pam30' or sub_matrix == 'pam 30' or sub_matrix == 'Pam30' or sub_matrix == 'Pam 30':
						pairwise_scores.append(PAM30[aa_i][aa_j])
					elif sub_matrix == 'PAM70' or sub_matrix == 'PAM 70' or sub_matrix == 'pam70' or sub_matrix == 'pam 70' or sub_matrix == 'Pam70' or sub_matrix == 'Pam 70':
						pairwise_scores.append(PAM70[aa_i][aa_j])
					elif sub_matrix == 'BLOSUM45' or sub_matrix == 'BLOSUM 45' or sub_matrix == 'blosum45' or sub_matrix == 'blosum 45' or sub_matrix == 'Blosum45' or sub_matrix == 'Blosum 45':
						pairwise_scores.append(BLOSUM45[aa_i][aa_j])
					elif sub_matrix == 'BLOSUM50' or sub_matrix == 'BLOSUM 50' or sub_matrix == 'blosum50' or sub_matrix == 'blosum 50' or sub_matrix == 'Blosum50' or sub_matrix == 'Blosum 50':
						pairwise_scores.append(BLOSUM50[aa_i][aa_j])
					elif sub_matrix == 'BLOSUM62' or sub_matrix == 'BLOSUM 62' or sub_matrix == 'blosum62' or sub_matrix == 'blosum 62' or sub_matrix == 'Blosum62' or sub_matrix == 'Blosum 62':
						pairwise_scores.append(BLOSUM62[aa_i][aa_j])
					elif sub_matrix == 'BLOSUM80' or sub_matrix == 'BLOSUM 80' or sub_matrix == 'blosum80' or sub_matrix == 'blosum 80' or sub_matrix == 'Blosum80' or sub_matrix == 'Blosum 80':
						pairwise_scores.append(BLOSUM80[aa_i][aa_j])
					elif sub_matrix == 'BLOSUM90' or sub_matrix == 'BLOSUM 90' or sub_matrix == 'blosum90' or sub_matrix == 'blosum 90' or sub_matrix == 'Blosum90' or sub_matrix == 'Blosum 90':
						pairwise_scores.append(BLOSUM90[aa_i][aa_j])
					else:
						sys.exit("Error: substitution matrix was not defined")
		return sum(pairwise_scores)
	#
	## Function to calculate the normalized sum of pairs (per aa position) - normalizing by the number of pairwise comparisons...
	def norm_sum_of_pairs(aa_list,sub_matrix,gap_pen): # Takes a list of aa characters for a particular alignment position, the substitution matrix and the gap penalty to use, and return the sum of pairs score
		if len(aa_list) > 1:
			pairwise_scores = []
			for i in range(len(aa_list)):
				for j in range(i + 1, len(aa_list)):
					#if aa_list[i] == '-' and aa_list[j] == '-':	# if both characters are gaps --> ignore...
					#	continue
					if aa_list[i] == '-' or aa_list[j] == '-':	# if one of them is a gap, penalize based on the defined gap penalty
					#elif aa_list[i] == '-' or aa_list[j] == '-':	# if one of them is a gap, penalize based on the defined gap penalty
						pairwise_scores.append(int(gap_pen))
					else:
						# Treat Selenocysteine (U) as Cysteine (C), and Pyrrolysine (O) as Any amino acid (X), as https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/multalignviewer/multalignviewer.html
						aa_i = aa_list[i]
						if aa_i == 'U':
							aa_i = 'C'
						elif aa_i == 'O':
							aa_i = 'X'
						aa_j = aa_list[j]
						if aa_j == 'U':
							aa_j = 'C'
						elif aa_j == 'O':
							aa_j = 'X'
						if sub_matrix == 'PAM1' or sub_matrix == 'PAM 1' or sub_matrix == 'pam1' or sub_matrix == 'pam 1' or sub_matrix == 'Pam1' or sub_matrix == 'Pam 1':
							pairwise_scores.append(PAM1[aa_i][aa_j])
						elif sub_matrix == 'PAM10' or sub_matrix == 'PAM 10' or sub_matrix == 'pam10' or sub_matrix == 'pam 10' or sub_matrix == 'Pam10' or sub_matrix == 'Pam 10':
							pairwise_scores.append(PAM10[aa_i][aa_j])
						elif sub_matrix == 'PAM30' or sub_matrix == 'PAM 30' or sub_matrix == 'pam30' or sub_matrix == 'pam 30' or sub_matrix == 'Pam30' or sub_matrix == 'Pam 30':
							pairwise_scores.append(PAM30[aa_i][aa_j])
						elif sub_matrix == 'PAM70' or sub_matrix == 'PAM 70' or sub_matrix == 'pam70' or sub_matrix == 'pam 70' or sub_matrix == 'Pam70' or sub_matrix == 'Pam 70':
							pairwise_scores.append(PAM70[aa_i][aa_j])
						elif sub_matrix == 'BLOSUM45' or sub_matrix == 'BLOSUM 45' or sub_matrix == 'blosum45' or sub_matrix == 'blosum 45' or sub_matrix == 'Blosum45' or sub_matrix == 'Blosum 45':
							pairwise_scores.append(BLOSUM45[aa_i][aa_j])
						elif sub_matrix == 'BLOSUM50' or sub_matrix == 'BLOSUM 50' or sub_matrix == 'blosum50' or sub_matrix == 'blosum 50' or sub_matrix == 'Blosum50' or sub_matrix == 'Blosum 50':
							pairwise_scores.append(BLOSUM50[aa_i][aa_j])
						elif sub_matrix == 'BLOSUM62' or sub_matrix == 'BLOSUM 62' or sub_matrix == 'blosum62' or sub_matrix == 'blosum 62' or sub_matrix == 'Blosum62' or sub_matrix == 'Blosum 62':
							pairwise_scores.append(BLOSUM62[aa_i][aa_j])
						elif sub_matrix == 'BLOSUM80' or sub_matrix == 'BLOSUM 80' or sub_matrix == 'blosum80' or sub_matrix == 'blosum 80' or sub_matrix == 'Blosum80' or sub_matrix == 'Blosum 80':
							pairwise_scores.append(BLOSUM80[aa_i][aa_j])
						elif sub_matrix == 'BLOSUM90' or sub_matrix == 'BLOSUM 90' or sub_matrix == 'blosum90' or sub_matrix == 'blosum 90' or sub_matrix == 'Blosum90' or sub_matrix == 'Blosum 90':
							pairwise_scores.append(BLOSUM90[aa_i][aa_j])
						else:
							sys.exit("Error: substitution matrix was not defined")
			if len(pairwise_scores)>0:
				return float(sum(pairwise_scores))/float(len(pairwise_scores))
			else:
				return 0
		else:
			return 0	# It should be 'NA', but put here 0, to be able to plot (it is 0 comparisons, therefore 0 score)
	#
	# Function to calculate the median value of a list of numeric values (int and/or float)
	def median(lst):
		n = len(lst)
		if n < 1:
			return None
		if n % 2 == 1:
			return sorted(lst)[n//2]
		else:
			return sum(sorted(lst)[n//2-1:n//2+1])/2.0
	#
	# Function to calculate the mean value of a list of numeric values (int and/or float)
	def mean(lst):
		n = len(lst)
		if n < 1:
			return None
		else:
			return float(sum(lst))/float(n)
	#
	#
	# Function to generate an list of tuples for a sliding window - i.e. input a list of objects, and it will generate a sliding window of them, e.g.: list(window([1,2,4,6,8,10],n=3)) will generate [(1, 2, 4), (2, 4, 6), (4, 6, 8), (6, 8, 10)]
	def window(seq, n=2):
		"Returns a sliding window (of width n) over data from the iterable"
		"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
		it = iter(seq)
		result = tuple(islice(it, n))
		if len(result) == n:
			yield result
		for elem in it:
			result = result[1:] + (elem,)
			yield result
	#
	#	
	# Function to average sliding window list of tuples. Will return a list of average values for each sliding window (ignores 'NA's, but if empty, will return 'NA')
	def ave_wind(wind_array):	# e.g. ave_wind(list(window([1,2,4,6,8,10],n=3)))
		ave = []
		for i in wind_array:
			new = [float(j) for j in list(i) if j != 'NA']	# get a list of floating points, and ignore 'NA' - check if getting an empty list (than return 'NA')
			if len(new) > 0:
				ave.append(mean(new))	# uses the mean function above...
			else:
				ave.append('NA')
		return ave
	#
	#	
	# Function to get the median sliding window list of tuples. Will return a list of median values for each sliding window (ignores 'NA's, but if empty, will return 'NA')
	def med_wind(wind_array):	# e.g. med_wind(list(window([1,2,4,6,8,10],n=3)))
		med = []
		for i in wind_array:
			new = [float(j) for j in list(i) if j != 'NA']	# get a list of floating points, and ignore 'NA' - check if getting an empty list (than return 'NA')
			if len(new) > 0:
				med.append(median(new))	# uses the median function above...
			else:
				med.append('NA')
		return med
	#		
	#
	## Sequence type:
	seq_type = str(args.sequence_type)	# understand if it is a protein or codon alignment sequence (codons need to be further translated)
	if seq_type == 'codon' or seq_type == 'Codon' or seq_type == 'CODON' or seq_type == 'codons' or seq_type == 'Codons' or seq_type == 'Codons':
		seq_type = 'codon'
	elif seq_type == 'protein' or seq_type == 'Protein' or seq_type == 'PROTEIN' or seq_type == 'Pep' or seq_type == 'PEP' or seq_type == 'aa' or seq_type == 'AA' or seq_type == 'p' or seq_type == 'P' or seq_type == 'a' or seq_type == 'A':
		seq_type = 'protein'
	else:
		sys.exit("Error: Please indicate sequence type, using the '-t' (or '--sequence-type') option. Should be either either 'protein' or 'codon' (default: 'protein')")
	#
	## Substitution matrix:
	sub_matrix = str(args.aa_sub_matrix)	# understand which substitution matrix to use for scoring
	if sub_matrix == 'PAM1' or sub_matrix == 'PAM 1' or sub_matrix == 'pam1' or sub_matrix == 'pam 1' or sub_matrix == 'Pam1' or sub_matrix == 'Pam 1':
		sub_matrix = 'PAM1'
	elif sub_matrix == 'PAM10' or sub_matrix == 'PAM 10' or sub_matrix == 'pam10' or sub_matrix == 'pam 10' or sub_matrix == 'Pam10' or sub_matrix == 'Pam 10':
		sub_matrix = 'PAM10'
	elif sub_matrix == 'PAM30' or sub_matrix == 'PAM 30' or sub_matrix == 'pam30' or sub_matrix == 'pam 30' or sub_matrix == 'Pam30' or sub_matrix == 'Pam 30':
		sub_matrix = 'PAM30'
	elif sub_matrix == 'PAM70' or sub_matrix == 'PAM 70' or sub_matrix == 'pam70' or sub_matrix == 'pam 70' or sub_matrix == 'Pam70' or sub_matrix == 'Pam 70':
		sub_matrix = 'PAM70'
	elif sub_matrix == 'BLOSUM45' or sub_matrix == 'BLOSUM 45' or sub_matrix == 'blosum45' or sub_matrix == 'blosum 45' or sub_matrix == 'Blosum45' or sub_matrix == 'Blosum 45':
		sub_matrix = 'BLOSUM45'
	elif sub_matrix == 'BLOSUM50' or sub_matrix == 'BLOSUM 50' or sub_matrix == 'blosum50' or sub_matrix == 'blosum 50' or sub_matrix == 'Blosum50' or sub_matrix == 'Blosum 50':
		sub_matrix = 'BLOSUM50'
	elif sub_matrix == 'BLOSUM62' or sub_matrix == 'BLOSUM 62' or sub_matrix == 'blosum62' or sub_matrix == 'blosum 62' or sub_matrix == 'Blosum62' or sub_matrix == 'Blosum 62':
		sub_matrix = 'BLOSUM62'
	elif sub_matrix == 'BLOSUM80' or sub_matrix == 'BLOSUM 80' or sub_matrix == 'blosum80' or sub_matrix == 'blosum 80' or sub_matrix == 'Blosum80' or sub_matrix == 'Blosum 80':
		sub_matrix = 'BLOSUM80'
	elif sub_matrix == 'BLOSUM90' or sub_matrix == 'BLOSUM 90' or sub_matrix == 'blosum90' or sub_matrix == 'blosum 90' or sub_matrix == 'Blosum90' or sub_matrix == 'Blosum 90':
		sub_matrix = 'BLOSUM90'
	else:
		sub_matrix = 'BLOSUM62'	
	#
	## Gap penalty:
	gap_pen = str(args.gap_penalty)	# understand which gap penalty to use for scoring
	if gap_pen.isalpha():
		sys.exit("Error: A non numeric string was provided for gap penalty. Please indicate only a numeric value for gap penalty using the '-gp' (or '--gap-penalty') option.\n Terminating the program!")
	#
	## Gap penalty:
	gap_pen = str(args.gap_penalty)	# understand which gap penalty to use for scoring
	if gap_pen.isalpha():
		sys.exit("Error: A non numeric string was provided for gap penalty. Please indicate only a numeric value for gap penalty using the '-gp' (or '--gap-penalty') option.\n Terminating the program!")
	#
	## Sliding window size:
	sw_size = str(args.sliding_window_size)	# understand which sliding window size to use for scoring and plotting (default 10 aa)
	if sw_size.isalpha():
		sys.exit("Error: A non numeric string was provided for sliding window size. Please indicate only a numeric value for sliding window size using the '-sw' (or '--sliding-window-size') option.\n Terminating the program!")
	#	
	## Find out which hydrophobicity scale to use - the default is 'Kyte-Doolittle'
	hydropho_scale = str(args.hydrophobicity_scale)
	if hydropho_scale == 'Kyte-Doolittle' or hydropho_scale == 'Kyte_Doolittle' or hydropho_scale == 'KD' or hydropho_scale == 'kd' or hydropho_scale == 'K-D' or hydropho_scale == 'K_D' or hydropho_scale == 'KYTE_DOOLITTLE'  or hydropho_scale == 'KYTE' or hydropho_scale == 'Kyte' or hydropho_scale == 'kyte' or hydropho_scale == 'kyte-doolittle' or hydropho_scale == 'kyte_doolittle':
		hydropho_scale = 'Hydrophobicity_Kyte-Doolittle'
	elif hydropho_scale == 'Hopp-Woods' or hydropho_scale == 'Hopp_Woods' or hydropho_scale == 'HW' or hydropho_scale == 'hw' or hydropho_scale == 'H-W' or hydropho_scale == 'H_W' or hydropho_scale == 'HOPP_WOODS'  or hydropho_scale == 'HOPP' or hydropho_scale == 'Hopp' or hydropho_scale == 'hopp' or hydropho_scale == 'hopp_woods' or hydropho_scale == 'hopp-woods':
		hydropho_scale = 'Hydrophobicity_Hopp-Woods'
	elif hydropho_scale == 'Cornette' or hydropho_scale == 'cornette' or hydropho_scale == 'C' or hydropho_scale == 'c' or hydropho_scale == 'CORNETTE':
		hydropho_scale = 'Hydrophobicity_Cornette'
	elif hydropho_scale == 'Eisenberg' or hydropho_scale == 'eisenberg' or hydropho_scale == 'E' or hydropho_scale == 'e' or hydropho_scale == 'EISENBERG':
		hydropho_scale = 'Hydrophobicity_Eisenberg'
	elif hydropho_scale == 'Rose' or hydropho_scale == 'rose' or hydropho_scale == 'R' or hydropho_scale == 'r' or hydropho_scale == 'ROSE':
		hydropho_scale = 'Hydrophobicity_Rose'
	elif hydropho_scale == 'Janin' or hydropho_scale == 'janin' or hydropho_scale == 'J' or hydropho_scale == 'j' or hydropho_scale == 'JANIN':
		hydropho_scale = 'Hydrophobicity_Janin'
	elif hydropho_scale == 'Engelman_GES' or hydropho_scale == 'Engelman-GES' or hydropho_scale == 'EGES' or hydropho_scale == 'eges' or hydropho_scale == 'EG' or hydropho_scale == 'eg' or hydropho_scale == 'E_G' or hydropho_scale == 'E-G' or hydropho_scale == 'ENGELMAN-GES'  or hydropho_scale == 'ENGELMAN' or hydropho_scale == 'Engelman' or hydropho_scale == 'engelman' or hydropho_scale == 'engelman_ges' or hydropho_scale == 'engelman-ges':
		hydropho_scale = 'Hydrophobicity_Engelman_GES'
	#
	## Find out which aa grouping classification to use - the default is 'Taylor'
	aa_grouping_method = str(args.aa_grouping)
	if aa_grouping_method == 'Katzir' or aa_grouping_method == 'KATZIR' or aa_grouping_method == 'K' or aa_grouping_method == 'k' or aa_grouping_method == 'KATZIR2006' or aa_grouping_method == 'Katzir2006' or aa_grouping_method == 'KATZIR_2006'  or aa_grouping_method == 'Katzir_2006' or aa_grouping_method == 'KATZIR-2006' or aa_grouping_method == 'Katzir-2006':
		aa_grouping_method = 'Katzir'
	elif aa_grouping_method == 'Murphy' or aa_grouping_method == 'MURPHY' or aa_grouping_method == 'M' or aa_grouping_method == 'm' or aa_grouping_method == 'Murphy2000' or aa_grouping_method == 'MURPHY2000' or aa_grouping_method == 'Murphy_2000'  or aa_grouping_method == 'MURPHY_2000' or aa_grouping_method == 'Murphy-2000' or aa_grouping_method == 'MURPHY-2000' or aa_grouping_method == 'Murphy_10' or aa_grouping_method == 'Murphy-10' or aa_grouping_method == 'Murphy_10_grouping':
		aa_grouping_method = 'Murphy'
	else:
		aa_grouping_method = 'Taylor'
	#
	## Get the list of taxa labels for each group:
	if args.species_groups_file is not None:
		with open(str(args.species_groups_file), 'r') as groupFile:
			species_lines = groupFile.read().splitlines()
		group_1 = []	# the first defined species is assigned to group1 - based on it's membership column, the script defines the groups 
		group_2 = []
		for s in species_lines:
			species = s.strip().split('\t')[0]
			group = s.strip().split('\t')[1]
			if group == species_lines[0].strip().split('\t')[1]:
				group_1.append(species)
			else:
				group_2.append(species)
	else:
		sys.exit("Error: Please provide a group membership file, using the '-g' (or '--species-groups-file') option. Should be a two lines file, each with a comma-sep list of taxa assosiated with either 'GROUP1=' or 'GROUP2='\nTerminating the program!")
	#
	species_label_to_name = {}
	if args.species_labels_names_file is not None:
		with open(str(args.species_labels_names_file), 'r') as spLabNamFile:
			for line in spLabNamFile:
				species_label_to_name[line.strip().split(str(args.species_labels_names_delimiter))[0]] = line.strip().split(str(args.species_labels_names_delimiter))[1]
	#
	#
	## Parse aa MSA fasta file:
	if args.fasta_msa_file is not None:
		with open(str(args.fasta_msa_file), 'r') as fasFile:
			fasta = fasFile.read().splitlines()
		if fasta[0][0] != ">": # check format
			sys.exit("Error: Please provide a multiple sequence alignmet in a fasta format, using the '-f' (or '--fasta-msa-file') option.\nTerminating the program!")
		seq_dict = {} # parse the input MSA fasta file into a dictionary {seqid:sequence}
		seqid = fasta[0].strip()[1:]
		seq = []
		for line in range(len(fasta)):
			if not fasta[line].startswith('>'):
				seq.append(fasta[line].strip())
				if line + 1 == len(fasta):
					seq_dict[seqid] = ''.join(seq)
			else:
				seq_dict[seqid] = ''.join(seq)
				seqid = fasta[line].strip()[1:]
				seq = []
	else:
		sys.exit("Error: Please provide a multiple sequence alignmet in a fasta format, using the '-f' (or '--fasta-msa-file') option.\nTerminating the program!")
	#
	## Translate codon file (if was specified) to protein sequence:
	new_seq_dict = {}
	if seq_type == 'codon':
		for s in seq_dict:
			split_to_codon_list = [seq_dict[s][i:i+3].upper() for i in range(0, len(seq_dict[s]), 3)]	# the .upper() makes sure that lowercase sequences will get converted to uppercase first...
			corresponding_aa_seq = ''.join([ codon.get(item,item) for item in split_to_codon_list])
			new_seq_dict[s] = corresponding_aa_seq
	else:
		new_seq_dict = seq_dict
	#
	## Separate the sequences of each group:
	group_1_seq_dict = {}
	for g1 in group_1:
		for id in new_seq_dict:
			if g1 in id:
				group_1_seq_dict[g1] = new_seq_dict[id]
	#
	group_2_seq_dict = {}
	for g2 in group_2:
		for id in new_seq_dict:
			if g2 in id:
				group_2_seq_dict[g2] = new_seq_dict[id]
	#
	#
	#
	## Traverse the alignment columns and compare the aa of each groups:
	criteria_satisfied_aa_res = ['AA_position\t' + 'Position_characterisitcs\t' + str(args.group1) + '_aa_(num_occurances)\t' + str(args.group2) + '_aa_(num_occurances)\t' + aa_grouping_method + '_physicochemical_properties_' + str(args.group1) + '_aa\t' + aa_grouping_method + '_physicochemical_properties_' + str(args.group2) + '_aa\t' + aa_grouping_method + '_physicochemical_properties_unique_' + str(args.group1) + '\t' + aa_grouping_method + '_physicochemical_properties_unique_' + str(args.group2)]	# store info for only aa that are the same in group1 and different (from group1 aa) in group2, start with the header
	all_aa_res = ['AA_position\t' + 'Position_characterisitcs\t' + str(args.group1) + '_aa_(num_occurances)\t' + str(args.group2) + '_aa_(num_occurances)\t' + aa_grouping_method + '_physicochemical_properties_' + str(args.group1) + '_aa\t'  + aa_grouping_method + '_physicochemical_properties_' + str(args.group2) + '_aa\t'  + aa_grouping_method + '_physicochemical_properties_unique_' + str(args.group1) + '\t'  + aa_grouping_method + '_physicochemical_properties_unique_' + str(args.group2)]	# store info for all aa positions. Start with the header - Currently not saving this (just used for debugging). If you want to save info for all aa positions, uncomment the saving block at the end
	all_aa_all_res = ['\t'.join(['AA_position' , 'Position_characterisitcs', str(args.group1) + '_species',str(args.group1) + '_aa_(num_occurances)',str(args.group2) + '_species',str(args.group2) + '_aa_(num_occurances)',aa_grouping_method + '_physicochemical_properties_' + str(args.group1) + '_aa',aa_grouping_method + '_physicochemical_properties_' + str(args.group2), aa_grouping_method + '_physicochemical_properties_unique_' + str(args.group1),aa_grouping_method + '_physicochemical_properties_unique_' + str(args.group2),str(args.group1) + '_aa_median_' + hydropho_scale ,str(args.group2) + '_aa_median_' + hydropho_scale, 'Total_sum_of_pairs_score_per_site', str(args.group1) + '_within-group_sum_of_pairs_score_per_site',str(args.group2) + '_within-group_sum_of_pairs_score_per_site','Between-groups_sum_of_pairs_score_per_site','Total_normalized_sum_of_pairs_score_per_site',str(args.group1) + '_within-group_normalized_sum_of_pairs_score_per_site',str(args.group2) + '_within-group_normalized_sum_of_pairs_score_per_site','Between-groups_normalized_sum_of_pairs_score_per_site'])]	# store all info for all aa positions, start with the header	
	# probably change the above output formats, to include the different aa groupings and hydrophobicity.
	#
	# This is an optional output (if the "a" option flag was indicated) of all the aa grouping comparisons:
	all_aa_grouping_res = ['\t'.join(['AA_position' , 'Position_characterisitcs', str(args.group1) + '_species',str(args.group1) + '_aa_(num_occurances)',str(args.group2) + '_species',str(args.group2) + '_aa_(num_occurances)', str(args.group1) + '_Katzir_2006_7_aa_groups', str(args.group2) + '_Katzir_2006_7_aa_groups', 'Common_' + str(args.group1) + '_' + str(args.group2) + '_Katzir_2006_7_aa_groups',str(args.group1) + '_Katzir_2006_radical_aa_replacement', str(args.group2) + '_Katzir_2006_radical_aa_replacement', 'Between_groups_Katzir_2006_radical_aa_replacement', str(args.group1) + '_Murphy_2000_10_aa_groups', str(args.group2) + '_Murphy_2000_10_aa_groups', 'Common_' + str(args.group1) + '_' + str(args.group2) + '_Murphy_2000_10_aa_groups',str(args.group1) + '_Murphy_2000_radical_aa_replacement', str(args.group2) + '_Murphy_2000_radical_aa_replacement', 'Between_groups_Murphy_2000_radical_aa_replacement', str(args.group1) + '_Zhang_2000_charge_classification', str(args.group2) + '_Zhang_2000_charge_classification', 'Common_' + str(args.group1) + '_' + str(args.group2) + '_Zhang_2000_charge_classification',str(args.group1) + '_Zhang_2000_charge_radical_aa_replacement', str(args.group2) + '_Zhang_2000_charge_radical_aa_replacement', 'Between_groups_Zhang_2000_charge_radical_aa_replacement', str(args.group1) + '_Zhang_2000_polarity_classification', str(args.group2) + '_Zhang_2000_polarity_classification', 'Common_' + str(args.group1) + '_' + str(args.group2) + '_Zhang_2000_polarity_classification',str(args.group1) + '_Zhang_2000_polarity_radical_aa_replacement', str(args.group2) + '_Zhang_2000_polarity_radical_aa_replacement', 'Between_groups_Zhang_2000_polarity_radical_aa_replacement', str(args.group1) + '_Zhang_2000_polarity_and_volumn_classification', str(args.group2) + '_Zhang_2000_polarity_and_volumn_classification', 'Common_' + str(args.group1) + '_' + str(args.group2) + '_Zhang_2000_polarity_and_volumn_classification',str(args.group1) + '_Zhang_2000_polarity_and_volumn_radical_aa_replacement', str(args.group2) + '_Zhang_2000_polarity_and_volumn_radical_aa_replacement', 'Between_groups_Zhang_2000_polarity_and_volumn_radical_aa_replacement', 	str(args.group1) + '_Betts_and_Russell_2003_less_exchangeable_aa', str(args.group2) + '_Betts_and_Russell_2003_less_exchangeable_aa', 'Common_' + str(args.group1) + '_' + str(args.group2) + '_Betts_and_Russell_2003_less_exchangeable_aa',str(args.group1) + '_Betts_and_Russell_2003_less_exchangeable_aa_radical_aa_replacement', str(args.group2) + '_Betts_and_Russell_2003_less_exchangeable_aa_radical_aa_replacement', 'Between_groups_Betts_and_Russell_2003_less_exchangeable_aa_radical_aa_replacement', str(args.group1) + '_Grantham_1974_physicochemical_max_distance',str(args.group2) + '_Grantham_1974_physicochemical_max_distance','Between_groups_Grantham_1974_physicochemical_max_distance',str(args.group1) + '_Grantham_1974_aa_radical_aa_replacement',str(args.group2) + '_Grantham_1974_aa_radical_aa_replacement','Between_groups_Grantham_1974_aa_radical_aa_replacement'])]
	# Between_groups_METHOD_X_radical_aa_replacement is if there are no common aa groups shared between the two groups (but only if there is at least one aa in each group...)
	# Grantham_1974_aa_radical_aa_replacement is when d >= 100 (based on Dagan et al., 2002)
	# Contains group comparisons of:
	# aa_physchem_prop_grouping_ERGPMM_book	(Katzir_2006_7_groups)
	# aa_physchem_prop_grouping_Murphy_2000	(Murphy_2000_10_groups)
	# physico_chemical_distance_Grantham (Grantham_1974_physicochemical_distance)
	# aa_charge_classification_Zhang (Zhang_2000_charge_classification)
	# aa_polarity_classification_Zhang	(Zhang_2000_polarity_classification)
	# aa_polarity_and_volumn_classification_Zhang	(Zhang_2000_polarity_and_volumn_classification)
	# Uniq_aa_Betts_and_Russell	(Betts_and_Russell_2003_less_exchangeable_aa)
	#
	#
	#
	Hydrophobicity_per_position_per_sp_dict = {}	# for hydrophobicity plot (the hydrophobicity value of each aa at each alignment position for each species (key))
	Hydrophobicity_per_position_dict = {}	# # store the mean and median (can add more estimates later...) hydrophobicity at each position (no groups, all the species) - use that afterward to also calculate the sliding window...
	Hydrophobicity_per_position_per_gr_dict_group1 = {}	# store the mean and median (can add more estimates later...) hydrophobicity for aa of group 1, at each position - use that afterward to also calculate the sliding window...
	Hydrophobicity_per_position_per_gr_dict_group2 = {}	# store the mean and median (can add more estimates later...) hydrophobicity for aa of group 2, at each position - use that afterward to also calculate the sliding window...
	sum_of_pairs_score_per_site_dict = {}	# Calculate the sum of pairs for each alignment position (User should indicate which substitution matrix to use, from a list, and the gap penalty (default = -11))
	sum_of_pairs_score_per_site_dict_group1 = {}	# Calculate the sum of pairs for each alignment position for group1 members
	sum_of_pairs_score_per_site_dict_group2 = {}	# Calculate the sum of pairs for each alignment position for group2 members
	sum_of_pairs_score_per_site_dict_between_groups = {}	# Calculate the sum of pairs between the two groups (not including within group pairs)
	norm_sum_of_pairs_score_per_site_dict = {}	# sum of pairs normalized by the number of pairwise comparisons
	norm_sum_of_pairs_score_per_site_dict_group1 = {}
	norm_sum_of_pairs_score_per_site_dict_group2 = {}
	norm_sum_of_pairs_score_per_site_dict_between_groups = {}
	sp_with_seq_per_site_dict_group1 = {}	# get the list of species with sequence per site for group 1
	sp_with_seq_per_site_dict_group2 = {}	# get the list of species with sequence per site for group 2
	#
	for c in range(len(new_seq_dict[list(new_seq_dict.keys())[0]])):	# For each aa position: 
		# For aa position use "c+1", as the iterator starts with 0...
		# Group1:
		group1_sp_list = []	# get the list of group 1 species with aa at position c
		group1_aa_list = []	# get the list of aa at position c for all group1 members 
		group1_aa_properties_Taylor = []	# get the list of aa properties from the Taylor classification (aa is assigned to multiple property groups - take the union of properties of all the aa of group1,  and intersect with the union of group2..., what they share and what's unique...)
		group1_aa_properties_Murphy_2000_10_groups = []	# get the list of aa groups (';' separate groups, ',' seperate aa within a classified group...) from the Murphy_2000_10_groups classification (aa is assigned to only one property group (set of amino acids...) - take the union of property groups of all group1 amino acids,  and intersect with the union of group2..., what they share and what's unique...)
		group1_aa_properties_ERGPMM_book_2006_7_groups = []	# get the list of aa groups (';' separate groups) from the ERGPMM_book_2006_7_groups classification (aa is assigned to only one property group - take the union of property groups of all group1 amino acids,  and intersect with the union of group2..., what they share and what's unique...)
		group1_aa_Hydrophobicity = []	# get the aa Hydrophobicity to calculate the median and mean (ignore gap and ambiguous characters - if all group members have a gap or ambiguous character, than median_Hydrophobicity or mean_Hydrophobicity will be equal to 'NA')
		#
		for gr1 in list(group_1_seq_dict):
			#print(gr1,c, len(group_1_seq_dict[gr1]))
			if group_1_seq_dict[gr1][c] not in ('-','*'):
				group1_sp_list.append(gr1)
			group1_aa_list.append(group_1_seq_dict[gr1][c])
			if group_1_seq_dict[gr1][c] not in ('-','X','*','U','O'): # For non-gap, non-ambiguous characters, or unusual amino acids (U: Selenocysteine, and O: Pyrrolysine) retrive the aa physiochemical properties
			#if (group_1_seq_dict[gr1][c] != '-') or (group_1_seq_dict[gr1][c] != 'X') or (group_1_seq_dict[gr1][c] != '*'):	# For non-gap, non-ambiguous characters, retrive the aa physiochemical properties
				# get the amino acid group
				#group1_aa_properties.extend(aa_physchem_prop[group_1_seq_dict[gr1][c]])
				group1_aa_properties_Taylor.extend(aa_physchem_prop_Taylor[group_1_seq_dict[gr1][c]])
				group1_aa_properties_Murphy_2000_10_groups.extend(aa_physchem_prop_grouping_Murphy_2000[group_1_seq_dict[gr1][c]])
				group1_aa_properties_ERGPMM_book_2006_7_groups.extend(aa_physchem_prop_grouping_ERGPMM_book[group_1_seq_dict[gr1][c]])
				#
				# get the Hydrophobicity values
				#group1_aa_Hydrophobicity.append(aa_Hydrophobicity[group_1_seq_dict[gr1][c]])
				if hydropho_scale == 'Hydrophobicity_Hopp-Woods':
					hydro = hydrophobicity_Hopp_Woods[group_1_seq_dict[gr1][c]]
				elif hydropho_scale == 'Hydrophobicity_Cornette':
					hydro = Hydrophobicity_Cornette[group_1_seq_dict[gr1][c]]
				elif hydropho_scale == 'Hydrophobicity_Eisenberg':
					hydro = Hydrophobicity_Eisenberg[group_1_seq_dict[gr1][c]]
				elif hydropho_scale == 'Hydrophobicity_Rose':
					hydro = Hydrophobicity_Rose[group_1_seq_dict[gr1][c]]
				elif hydropho_scale == 'Hydrophobicity_Janin':
					hydro = Hydrophobicity_Janin[group_1_seq_dict[gr1][c]]
				elif hydropho_scale == 'Hydrophobicity_Engelman_GES':
					hydro = Hydrophobicity_Engelman_GES[group_1_seq_dict[gr1][c]]
				else: # the default is Hydrophobicity_Kyte-Doolittle
					hydro = hydrophobicity_Kyte_Doolittle[group_1_seq_dict[gr1][c]]
				#
				group1_aa_Hydrophobicity.append(hydro)
				#
				if gr1 in Hydrophobicity_per_position_per_sp_dict:
					Hydrophobicity_per_position_per_sp_dict[gr1][str(c+1)] = hydro
				else:
					Hydrophobicity_per_position_per_sp_dict[gr1] = {str(c+1): hydro}
			else:
				if gr1 in Hydrophobicity_per_position_per_sp_dict:
					Hydrophobicity_per_position_per_sp_dict[gr1][str(c+1)] = 'NA'
				else:
					Hydrophobicity_per_position_per_sp_dict[gr1] = {str(c+1): 'NA'}
		#
		# Calculate sum_of_pairs scores per site group 1:
		#
		sum_of_pairs_score_per_site_group1 = sum_of_pairs(group1_aa_list,sub_matrix,gap_pen)
		sum_of_pairs_score_per_site_dict_group1[str(c+1)] = sum_of_pairs_score_per_site_group1
		norm_sum_of_pairs_score_per_site_group1 = norm_sum_of_pairs(group1_aa_list,sub_matrix,gap_pen)
		norm_sum_of_pairs_score_per_site_dict_group1[str(c+1)] = norm_sum_of_pairs_score_per_site_group1
		#
		# Calculate the mean and median aa hydrophobicity for group1:
		#
		if len(group1_aa_Hydrophobicity) > 0:	# calculate the median aa Hydrophobicity
			group1_median_aa_Hydrophobicity = str(median(group1_aa_Hydrophobicity))
			group1_mean_aa_Hydrophobicity = str(mean(group1_aa_Hydrophobicity))
		else:
			group1_median_aa_Hydrophobicity = 'NA'
			group1_mean_aa_Hydrophobicity = 'NA'
		#
		Hydrophobicity_per_position_per_gr_dict_group1[str(c+1)] = {'Median':group1_median_aa_Hydrophobicity,'Mean':group1_mean_aa_Hydrophobicity}	# store the mean and median (can add more estimates later...) hydrophobicity for aa of group 1, at each position - use that afterward to also calculate the sliding window...
		#
		sp_with_seq_per_site_dict_group1[c+1] = ';'.join(group1_sp_list)
		#
		# Group2:
		group2_sp_list = []	# get the list of group 2 species with aa at position c
		group2_aa_list = []	# get the list of aa at position c for all group2 members 
		group2_aa_properties_Taylor = []	# get the list of aa properties from the Taylor classification (aa is assigned to multiple property groups - take the union of properties of all the aa of group2,  and intersect with the union of group2..., what they share and what's unique...)
		group2_aa_properties_Murphy_2000_10_groups = []	# get the list of aa groups (';' separate groups, ',' seperate aa within a classified group...) from the Murphy_2000_10_groups classification (aa is assigned to only one property group (set of amino acids...) - take the union of property groups of all group2 amino acids,  and intersect with the union of group2..., what they share and what's unique...)
		group2_aa_properties_ERGPMM_book_2006_7_groups = []	# get the list of aa groups (';' separate groups) from the ERGPMM_book_2006_7_groups classification (aa is assigned to only one property group - take the union of property groups of all group2 amino acids,  and intersect with the union of group2..., what they share and what's unique...)
		group2_aa_Hydrophobicity = []	# get the aa Hydrophobicity to calculate the median and mean (ignore gap and ambiguous characters - if all group members have a gap or ambiguous character, than median_Hydrophobicity or mean_Hydrophobicity will be equal to 'NA')
		for gr2 in list(group_2_seq_dict):
			if group_2_seq_dict[gr2][c] not in ('-','*'):
				group2_sp_list.append(gr2)
			group2_aa_list.append(group_2_seq_dict[gr2][c])
			if group_2_seq_dict[gr2][c] not in ('-','X','*','U','O'): # For non-gap, non-ambiguous characters, or unusual amino acids (U: Selenocysteine, and O: Pyrrolysine) retrive the aa physiochemical properties
				# get the amino acid group
				group2_aa_properties_Taylor.extend(aa_physchem_prop_Taylor[group_2_seq_dict[gr2][c]])
				group2_aa_properties_Murphy_2000_10_groups.extend(aa_physchem_prop_grouping_Murphy_2000[group_2_seq_dict[gr2][c]])
				group2_aa_properties_ERGPMM_book_2006_7_groups.extend(aa_physchem_prop_grouping_ERGPMM_book[group_2_seq_dict[gr2][c]])
				# get the Hydrophobicity values
				if hydropho_scale == 'Hydrophobicity_Hopp-Woods':
					hydro = hydrophobicity_Hopp_Woods[group_2_seq_dict[gr2][c]]
				elif hydropho_scale == 'Hydrophobicity_Cornette':
					hydro = Hydrophobicity_Cornette[group_2_seq_dict[gr2][c]]
				elif hydropho_scale == 'Hydrophobicity_Eisenberg':
					hydro = Hydrophobicity_Eisenberg[group_2_seq_dict[gr2][c]]
				elif hydropho_scale == 'Hydrophobicity_Rose':
					hydro = Hydrophobicity_Rose[group_2_seq_dict[gr2][c]]
				elif hydropho_scale == 'Hydrophobicity_Janin':
					hydro = Hydrophobicity_Janin[group_2_seq_dict[gr2][c]]
				elif hydropho_scale == 'Hydrophobicity_Engelman_GES':
					hydro = Hydrophobicity_Engelman_GES[group_2_seq_dict[gr2][c]]
				else: # the default is Hydrophobicity_Kyte-Doolittle
					hydro = hydrophobicity_Kyte_Doolittle[group_2_seq_dict[gr2][c]]
				#
				group2_aa_Hydrophobicity.append(hydro)
				#
				if gr2 in Hydrophobicity_per_position_per_sp_dict:
					Hydrophobicity_per_position_per_sp_dict[gr2][str(c+1)] = hydro
				else:
					Hydrophobicity_per_position_per_sp_dict[gr2] = {str(c+1): hydro}
			else:
				if gr2 in Hydrophobicity_per_position_per_sp_dict:
					Hydrophobicity_per_position_per_sp_dict[gr2][str(c+1)] = 'NA'
				else:
					Hydrophobicity_per_position_per_sp_dict[gr2] = {str(c+1): 'NA'}
		#
		# Calculate sum_of_pairs scores per site group 2:
		#
		sum_of_pairs_score_per_site_group2 = sum_of_pairs(group2_aa_list,sub_matrix,gap_pen)
		sum_of_pairs_score_per_site_dict_group2[str(c+1)] = sum_of_pairs_score_per_site_group2
		norm_sum_of_pairs_score_per_site_group2 = norm_sum_of_pairs(group2_aa_list,sub_matrix,gap_pen)
		norm_sum_of_pairs_score_per_site_dict_group2[str(c+1)] = norm_sum_of_pairs_score_per_site_group2
		#
		# Calculate the mean and median aa hydrophobicity for group2:
		#
		if len(group2_aa_Hydrophobicity) > 0:	# calculate the median aa Hydrophobicity
			group2_median_aa_Hydrophobicity = str(median(group2_aa_Hydrophobicity))
			group2_mean_aa_Hydrophobicity = str(mean(group2_aa_Hydrophobicity))
		else:
			group2_median_aa_Hydrophobicity = 'NA'
			group2_mean_aa_Hydrophobicity = 'NA'
		#
		Hydrophobicity_per_position_per_gr_dict_group2[str(c+1)] = {'Median':group2_median_aa_Hydrophobicity,'Mean':group2_mean_aa_Hydrophobicity}	# store the mean and median (can add more estimates later...) hydrophobicity for aa of group 2, at each position - use that afterward to also calculate the sliding window...
		#
		sp_with_seq_per_site_dict_group2[c+1] = ';'.join(group2_sp_list)
		#
		#
		#
		# Both groups calculations:
		#
		# Calculate the mean and median aa hydrophobicity for all species together:
		aa_Hydrophobicity = group1_aa_Hydrophobicity + group2_aa_Hydrophobicity
		if len(aa_Hydrophobicity) > 0:
			median_aa_Hydrophobicity = str(median(aa_Hydrophobicity))
			mean_aa_Hydrophobicity = str(mean(aa_Hydrophobicity))
		else:
			median_aa_Hydrophobicity = 'NA'
			mean_aa_Hydrophobicity = 'NA'
		#
		Hydrophobicity_per_position_dict[str(c+1)] = {'Median':median_aa_Hydrophobicity,'Mean':mean_aa_Hydrophobicity}	# store the mean and median (can add more estimates later...) hydrophobicity for aa of group 1, at each position - use that afterward to also calculate the sliding window...
		#
		# Calculate sum_of_pairs scores per site both and between groups:
		#
		both_group_aa_list = group1_aa_list + group2_aa_list
		sum_of_pairs_score_per_site = sum_of_pairs(both_group_aa_list,sub_matrix,gap_pen)
		sum_of_pairs_score_per_site_dict[str(c+1)] = sum_of_pairs_score_per_site
		norm_sum_of_pairs_score_per_site = norm_sum_of_pairs(both_group_aa_list,sub_matrix,gap_pen)
		norm_sum_of_pairs_score_per_site_dict[str(c+1)] = norm_sum_of_pairs_score_per_site
		# Calculate the score between groups (one from each lists as pairs - not like sum_of_pairs_score_per_site_dict which contain all comparisons...)
		between_scores = []	# store all between group scores (sum them up)
		for i in group1_aa_list:
			for j in group2_aa_list:
				between_scores.append(sum_of_pairs([i,j],sub_matrix,gap_pen))
		sum_of_pairs_score_per_site_between_groups = sum(between_scores)
		sum_of_pairs_score_per_site_dict_between_groups[str(c+1)] = sum_of_pairs_score_per_site_between_groups
		if len(between_scores) != 0:
			norm_sum_of_pairs_score_per_site_between_groups = float(sum(between_scores))/float(len(between_scores))
		else:
			norm_sum_of_pairs_score_per_site_between_groups = float(0)
		norm_sum_of_pairs_score_per_site_dict_between_groups[str(c+1)] = norm_sum_of_pairs_score_per_site_between_groups
		#
		#
		# Get aa frequencies:
		aa_freq_group1 = dict((x,group1_aa_list.count(x)) for x in set(group1_aa_list))	# get the counts (number of occurances) per aa in the group1 list
		group1_aa_list_with_counts = []
		for a in list(aa_freq_group1.keys()):
			group1_aa_list_with_counts.append(a + '(' + str(aa_freq_group1[a]) + ')')
		aa_freq_group2 = dict((x,group2_aa_list.count(x)) for x in set(group2_aa_list))	# get the counts (number of occurances) per aa in the group2 list
		group2_aa_list_with_counts = []
		for a in list(aa_freq_group2.keys()):
			group2_aa_list_with_counts.append(a + '(' + str(aa_freq_group2[a]) + ')')
		#
		## For the aa_grouping output:
		group1_Katzir_2006_7_groups = []
		group2_Katzir_2006_7_groups = []
		group1_Murphy_2000_10_groups = []
		group2_Murphy_2000_10_groups = []
		group1_Zhang_2000_charge_classification = []
		group2_Zhang_2000_charge_classification = []
		group1_Zhang_2000_polarity_classification = []
		group2_Zhang_2000_polarity_classification = []
		group1_Zhang_2000_polarity_and_volumn_classification = []
		group2_Zhang_2000_polarity_and_volumn_classification = []
		group1_Betts_and_Russell_2003_less_exchangeable_aa = []
		group2_Betts_and_Russell_2003_less_exchangeable_aa = []
		within_group1_Grantham_1974_physicochemical_max_distance = 0
		within_group2_Grantham_1974_physicochemical_max_distance = 0
		between_group1_group2_Grantham_1974_physicochemical_max_distance = 0
		#
		group1_aa_list_uniq_and_valid = list(set([i for i in group1_aa_list if i not in ('-','X','*','U','O')]))	# Ignore gaps, ambiguous characters or unusual amino acids (U: Selenocysteine, and O: Pyrrolysine), and keep a unique set
		group2_aa_list_uniq_and_valid = list(set([j for j in group2_aa_list if j not in ('-','X','*','U','O')]))	# Ignore gaps, ambiguous characters or unusual amino acids (U: Selenocysteine, and O: Pyrrolysine), and keep a unique set
		# Within group1:
		if len(group1_aa_list_uniq_and_valid) > 0:
			within_group1_Grantham_1974_physicochemical_max_distance = max([physico_chemical_distance_Grantham[x][y] for x in group1_aa_list_uniq_and_valid for y in group1_aa_list_uniq_and_valid])
			if within_group1_Grantham_1974_physicochemical_max_distance >= 100:
				group1_Grantham_radical_aa_replacement = 'Yes'
			else:
				group1_Grantham_radical_aa_replacement = 'No'
			for i in group1_aa_list_uniq_and_valid:
				group1_Katzir_2006_7_groups.append(aa_physchem_prop_grouping_ERGPMM_book[i])
				group1_Murphy_2000_10_groups.append(aa_physchem_prop_grouping_Murphy_2000[i])
				group1_Zhang_2000_charge_classification.append(aa_charge_classification_Zhang[i])
				group1_Zhang_2000_polarity_classification.append(aa_polarity_classification_Zhang[i])
				group1_Zhang_2000_polarity_and_volumn_classification.append(aa_polarity_and_volumn_classification_Zhang[i])
				if i in Uniq_aa_Betts_and_Russell:
					group1_Betts_and_Russell_2003_less_exchangeable_aa.append(i)
				else:
					group1_Betts_and_Russell_2003_less_exchangeable_aa.append('other')
			if len(set(group1_Katzir_2006_7_groups)) == 0:
				group1_Katzir = 'NA'
				group1_Katzir_radical = 'NA'
			elif len(set(group1_Katzir_2006_7_groups)) == 1:
				group1_Katzir = ';'.join(list(set(group1_Katzir_2006_7_groups)))
				group1_Katzir_radical = 'No'
			else:
				group1_Katzir = ';'.join(list(set(group1_Katzir_2006_7_groups)))
				group1_Katzir_radical = 'Yes'

			if len(set(group1_Murphy_2000_10_groups)) == 0:
				group1_Murphy = 'NA'
				group1_Murphy_radical = 'NA'
			elif len(set(group1_Murphy_2000_10_groups)) == 1:
				group1_Murphy = ';'.join(list(set(group1_Murphy_2000_10_groups)))
				group1_Murphy_radical = 'No'
			else:
				group1_Murphy = ';'.join(list(set(group1_Murphy_2000_10_groups)))
				group1_Murphy_radical = 'Yes'

			if len(set(group1_Zhang_2000_charge_classification)) == 0:
				group1_Zhang_charge = 'NA'
				group1_Zhang_charge_radical = 'NA'
			elif len(set(group1_Zhang_2000_charge_classification)) == 1:
				group1_Zhang_charge = ';'.join(list(set(group1_Zhang_2000_charge_classification)))
				group1_Zhang_charge_radical = 'No'
			else:
				group1_Zhang_charge = ';'.join(list(set(group1_Zhang_2000_charge_classification)))
				group1_Zhang_charge_radical = 'Yes'		

			if len(set(group1_Zhang_2000_polarity_classification)) == 0:
				group1_Zhang_polarity = 'NA'
				group1_Zhang_polarity_radical = 'NA'
			elif len(set(group1_Zhang_2000_polarity_classification)) == 1:
				group1_Zhang_polarity = ';'.join(list(set(group1_Zhang_2000_polarity_classification)))
				group1_Zhang_polarity_radical = 'No'
			else:
				group1_Zhang_polarity = ';'.join(list(set(group1_Zhang_2000_polarity_classification)))
				group1_Zhang_polarity_radical = 'Yes'

			if len(set(group1_Zhang_2000_polarity_and_volumn_classification)) == 0:
				group1_Zhang_polarity_and_volumn = 'NA'
				group1_Zhang_polarity_and_volumn_radical = 'NA'
			elif len(set(group1_Zhang_2000_polarity_and_volumn_classification)) == 1:
				group1_Zhang_polarity_and_volumn = ';'.join(list(set(group1_Zhang_2000_polarity_and_volumn_classification)))
				group1_Zhang_polarity_and_volumn_radical = 'No'
			else:
				group1_Zhang_polarity_and_volumn = ';'.join(list(set(group1_Zhang_2000_polarity_and_volumn_classification)))
				group1_Zhang_polarity_and_volumn_radical = 'Yes'
			if len(set(group1_Betts_and_Russell_2003_less_exchangeable_aa)) == 0:
				group1_Betts_and_Russell_less_exchangeable_aa = 'NA'
				group1_Betts_and_Russell_less_exchangeable_aa_radical = 'NA'
			elif len(set(group1_Betts_and_Russell_2003_less_exchangeable_aa)) == 1:
				group1_Betts_and_Russell_less_exchangeable_aa = ';'.join(list(set(group1_Betts_and_Russell_2003_less_exchangeable_aa)))
				group1_Betts_and_Russell_less_exchangeable_aa_radical = 'No'
			else:
				group1_Betts_and_Russell_less_exchangeable_aa = ';'.join(list(set(group1_Betts_and_Russell_2003_less_exchangeable_aa)))
				group1_Betts_and_Russell_less_exchangeable_aa_radical = 'Yes'
		else:
			group1_Katzir = 'NA'
			group1_Katzir_radical = 'NA'
			group1_Murphy = 'NA'
			group1_Murphy_radical = 'NA'
			group1_Zhang_charge = 'NA'
			group1_Zhang_charge_radical = 'NA'
			group1_Zhang_polarity = 'NA'
			group1_Zhang_polarity_radical = 'NA'
			group1_Zhang_polarity_and_volumn = 'NA'
			group1_Zhang_polarity_and_volumn_radical = 'NA'
			group1_Betts_and_Russell_less_exchangeable_aa = 'NA'
			group1_Betts_and_Russell_less_exchangeable_aa_radical = 'NA'
			group1_Grantham_radical_aa_replacement = 'NA'
		#
		# Within group2:
		if len(group2_aa_list_uniq_and_valid) > 0:
			within_group2_Grantham_1974_physicochemical_max_distance = max([physico_chemical_distance_Grantham[x][y] for x in group2_aa_list_uniq_and_valid for y in group2_aa_list_uniq_and_valid])
			if within_group2_Grantham_1974_physicochemical_max_distance >= 100:
				group2_Grantham_radical_aa_replacement = 'Yes'
			else:
				group2_Grantham_radical_aa_replacement = 'No'
			for j in group2_aa_list_uniq_and_valid:
				group2_Katzir_2006_7_groups.append(aa_physchem_prop_grouping_ERGPMM_book[j])
				group2_Murphy_2000_10_groups.append(aa_physchem_prop_grouping_Murphy_2000[j])
				group2_Zhang_2000_charge_classification.append(aa_charge_classification_Zhang[j])
				group2_Zhang_2000_polarity_classification.append(aa_polarity_classification_Zhang[j])
				group2_Zhang_2000_polarity_and_volumn_classification.append(aa_polarity_and_volumn_classification_Zhang[j])
				if j in Uniq_aa_Betts_and_Russell:
					group2_Betts_and_Russell_2003_less_exchangeable_aa.append(j)
				else:
					group2_Betts_and_Russell_2003_less_exchangeable_aa.append('other')
			if len(set(group2_Katzir_2006_7_groups)) == 0:
				group2_Katzir = 'NA'
				group2_Katzir_radical = 'NA'
			elif len(set(group2_Katzir_2006_7_groups)) == 1:
				group2_Katzir = ';'.join(list(set(group2_Katzir_2006_7_groups)))
				group2_Katzir_radical = 'No'
			else:
				group2_Katzir = ';'.join(list(set(group2_Katzir_2006_7_groups)))
				group2_Katzir_radical = 'Yes'

			if len(set(group2_Murphy_2000_10_groups)) == 0:
				group2_Murphy = 'NA'
				group2_Murphy_radical = 'NA'
			elif len(set(group2_Murphy_2000_10_groups)) == 1:
				group2_Murphy = ';'.join(list(set(group2_Murphy_2000_10_groups)))
				group2_Murphy_radical = 'No'
			else:
				group2_Murphy = ';'.join(list(set(group2_Murphy_2000_10_groups)))
				group2_Murphy_radical = 'Yes'

			if len(set(group2_Zhang_2000_charge_classification)) == 0:
				group2_Zhang_charge = 'NA'
				group2_Zhang_charge_radical = 'NA'
			elif len(set(group2_Zhang_2000_charge_classification)) == 1:
				group2_Zhang_charge = ';'.join(list(set(group2_Zhang_2000_charge_classification)))
				group2_Zhang_charge_radical = 'No'
			else:
				group2_Zhang_charge = ';'.join(list(set(group2_Zhang_2000_charge_classification)))
				group2_Zhang_charge_radical = 'Yes'		

			if len(set(group2_Zhang_2000_polarity_classification)) == 0:
				group2_Zhang_polarity = 'NA'
				group2_Zhang_polarity_radical = 'NA'
			elif len(set(group2_Zhang_2000_polarity_classification)) == 1:
				group2_Zhang_polarity = ';'.join(list(set(group2_Zhang_2000_polarity_classification)))
				group2_Zhang_polarity_radical = 'No'
			else:
				group2_Zhang_polarity = ';'.join(list(set(group2_Zhang_2000_polarity_classification)))
				group2_Zhang_polarity_radical = 'Yes'

			if len(set(group2_Zhang_2000_polarity_and_volumn_classification)) == 0:
				group2_Zhang_polarity_and_volumn = 'NA'
				group2_Zhang_polarity_and_volumn_radical = 'NA'
			elif len(set(group2_Zhang_2000_polarity_and_volumn_classification)) == 1:
				group2_Zhang_polarity_and_volumn = ';'.join(list(set(group2_Zhang_2000_polarity_and_volumn_classification)))
				group2_Zhang_polarity_and_volumn_radical = 'No'
			else:
				group2_Zhang_polarity_and_volumn = ';'.join(list(set(group2_Zhang_2000_polarity_and_volumn_classification)))
				group2_Zhang_polarity_and_volumn_radical = 'Yes'
			if len(set(group2_Betts_and_Russell_2003_less_exchangeable_aa)) == 0:
				group2_Betts_and_Russell_less_exchangeable_aa = 'NA'
				group2_Betts_and_Russell_less_exchangeable_aa_radical = 'NA'
			elif len(set(group2_Betts_and_Russell_2003_less_exchangeable_aa)) == 1:
				group2_Betts_and_Russell_less_exchangeable_aa = ';'.join(list(set(group2_Betts_and_Russell_2003_less_exchangeable_aa)))
				group2_Betts_and_Russell_less_exchangeable_aa_radical = 'No'
			else:
				group2_Betts_and_Russell_less_exchangeable_aa = ';'.join(list(set(group2_Betts_and_Russell_2003_less_exchangeable_aa)))
				group2_Betts_and_Russell_less_exchangeable_aa_radical = 'Yes'
		else:
			group2_Katzir = 'NA'
			group2_Katzir_radical = 'NA'
			group2_Murphy = 'NA'
			group2_Murphy_radical = 'NA'
			group2_Zhang_charge = 'NA'
			group2_Zhang_charge_radical = 'NA'
			group2_Zhang_polarity = 'NA'
			group2_Zhang_polarity_radical = 'NA'
			group2_Zhang_polarity_and_volumn = 'NA'
			group2_Zhang_polarity_and_volumn_radical = 'NA'
			group2_Betts_and_Russell_less_exchangeable_aa = 'NA'
			group2_Betts_and_Russell_less_exchangeable_aa_radical = 'NA'
			group2_Grantham_radical_aa_replacement = 'NA'
		#
		# Between group1 and group2:
		# Between_groups_XXX_radical_aa_replacement set as 'Yes' if there are no common aa groups shared between the two groups (but only if there is at least one aa in each group...)
		if len(group1_aa_list_uniq_and_valid) > 0 and len(group2_aa_list_uniq_and_valid) > 0:
			#Grantham:
			between_group1_group2_Grantham_1974_physicochemical_max_distance = max([physico_chemical_distance_Grantham[x][y] for x in group1_aa_list_uniq_and_valid for y in group2_aa_list_uniq_and_valid])
			if between_group1_group2_Grantham_1974_physicochemical_max_distance >= 100:
				between_group1_group2_Grantham_radical_aa_replacement = 'Yes'
			else:
				between_group1_group2_Grantham_radical_aa_replacement = 'No'
			#
			# Katzir:
			group1_group2_common_aa_groups_Katzir = list(set(group1_Katzir_2006_7_groups) & set(group2_Katzir_2006_7_groups))
			if len(group1_group2_common_aa_groups_Katzir) == 0:
				group1_group2_common_aa_groups_Katzir = 'NA'
				between_groups_radiacal_replacement_Katzir = 'Yes'
			else:
				group1_group2_common_aa_groups_Katzir = ';'.join(group1_group2_common_aa_groups_Katzir)
				between_groups_radiacal_replacement_Katzir = 'No'
			# Murphy:
			group1_group2_common_aa_groups_Murphy = list(set(group1_Murphy_2000_10_groups) & set(group2_Murphy_2000_10_groups))
			if len(group1_group2_common_aa_groups_Murphy) == 0:
				group1_group2_common_aa_groups_Murphy = 'NA'
				between_groups_radiacal_replacement_Murphy = 'Yes'
			else:
				group1_group2_common_aa_groups_Murphy = ';'.join(group1_group2_common_aa_groups_Murphy)
				between_groups_radiacal_replacement_Murphy = 'No'
			# Zhang_2000_charge:
			group1_group2_common_aa_groups_Zhang_charge = list(set(group1_Zhang_2000_charge_classification) & set(group2_Zhang_2000_charge_classification))
			if len(group1_group2_common_aa_groups_Zhang_charge) == 0:
				group1_group2_common_aa_groups_Zhang_charge = 'NA'
				between_groups_radiacal_replacement_Zhang_charge = 'Yes'
			else:
				group1_group2_common_aa_groups_Zhang_charge = ';'.join(group1_group2_common_aa_groups_Zhang_charge)
				between_groups_radiacal_replacement_Zhang_charge = 'No'
			# Zhang_2000_polarity:
			group1_group2_common_aa_groups_Zhang_polarity = list(set(group1_Zhang_2000_polarity_classification) & set(group2_Zhang_2000_polarity_classification))
			if len(group1_group2_common_aa_groups_Zhang_polarity) == 0:
				group1_group2_common_aa_groups_Zhang_polarity = 'NA'
				between_groups_radiacal_replacement_Zhang_polarity = 'Yes'
			else:
				group1_group2_common_aa_groups_Zhang_polarity = ';'.join(group1_group2_common_aa_groups_Zhang_polarity)
				between_groups_radiacal_replacement_Zhang_polarity = 'No'
			# Zhang_2000_polarity_and_volumn:
			group1_group2_common_aa_groups_Zhang_polarity_and_volumn = list(set(group1_Zhang_2000_polarity_and_volumn_classification) & set(group2_Zhang_2000_polarity_and_volumn_classification))
			if len(group1_group2_common_aa_groups_Zhang_polarity_and_volumn) == 0:
				group1_group2_common_aa_groups_Zhang_polarity_and_volumn = 'NA'
				between_groups_radiacal_replacement_Zhang_polarity_and_volumn = 'Yes'
			else:
				group1_group2_common_aa_groups_Zhang_polarity_and_volumn = ';'.join(group1_group2_common_aa_groups_Zhang_polarity_and_volumn)
				between_groups_radiacal_replacement_Zhang_polarity_and_volumn = 'No'
			# Betts_and_Russell_2003
			group1_group2_common_aa_groups_Betts_and_Russell = list(set(group1_Betts_and_Russell_2003_less_exchangeable_aa) & set(group2_Betts_and_Russell_2003_less_exchangeable_aa))
			if len(group1_group2_common_aa_groups_Betts_and_Russell) == 0:
				group1_group2_common_aa_groups_Betts_and_Russell = 'NA'
				between_groups_radiacal_replacement_Betts_and_Russell = 'Yes'
			else:
				group1_group2_common_aa_groups_Betts_and_Russell = ';'.join(group1_group2_common_aa_groups_Betts_and_Russell)
				between_groups_radiacal_replacement_Betts_and_Russell = 'No'
		else:
			group1_group2_common_aa_groups_Katzir = 'NA'
			between_groups_radiacal_replacement_Katzir = 'NA'		
			group1_group2_common_aa_groups_Murphy = 'NA'
			between_groups_radiacal_replacement_Murphy = 'NA'
			group1_group2_common_aa_groups_Zhang_charge = 'NA'
			between_groups_radiacal_replacement_Zhang_charge = 'NA'
			group1_group2_common_aa_groups_Zhang_polarity = 'NA'
			between_groups_radiacal_replacement_Zhang_polarity = 'NA'
			group1_group2_common_aa_groups_Zhang_polarity_and_volumn = 'NA'
			between_groups_radiacal_replacement_Zhang_polarity_and_volumn = 'NA'
			group1_group2_common_aa_groups_Betts_and_Russell = 'NA'
			between_groups_radiacal_replacement_Betts_and_Russell = 'NA'
			between_group1_group2_Grantham_radical_aa_replacement = 'NA'
		#
		#
		## Define the group comparisons classes:
		# Relationships:
		# All_gaps (all the species have a gap symbol for that position - probably the alignment was extracted from a larger alignment, and all gap positions were not removed...)
		# Exist_only_in_group1_conserved (only group 1 members have an aa, and it is conserved, group 2 members all have gaps..)
		# Exist_only_in_group1_polymorphic (only group 1 members have an aa, and it is polymorphic, group 2 members all have gaps..)
		# Exist_only_in_group2_conserved (only group 2 members have an aa, and it is conserved, group 1 members all have gaps..)
		# Exist_only_in_group2_polymorphic (only group 2 members have an aa, and it is polymorphic, group 1 members all have gaps..)
		# Conserved_in_all_species (same one amino acid in all species)
		# Diverged_in_all_species (no common amino acid between any of the species)
		# Conserved_within_groups_diverged_between_groups (same one amino acid in each group, that is different across groups)
		# Polymorphic_within_groups_diverged_between_groups (Multiple amino acids in each group, that are different across groups - not shared among groups...)
		# Polymorphic_within_groups_shared_between_groups (Multiple amino acids in each group, but the same ones in both groups (common polymorphism))
		# Polymorphic_within_groups_shared_and_uniq_between_groups (Multiple amino acids in each group, where some of the amino acids are common and some that are unique to a particular group)
		# Polymorphic_within_group1_diverged_between_groups_conserved_within_group2 (Multiple amino acids in group 1 (polymorphic), that don't exist in group 2 (therefore 'diverged'), where within group2, the same one amino acid exist in all members (conserved))
		# Polymorphic_within_group1_shared_between_groups_conserved_within_group2 (Multiple amino acids in group 1 (polymorphic), that some are also in group 2 (therefore 'shared'), where within group2, the same one amino acid exist in all members (conserved))
		# Conserved_within_group1_diverged_between_groups_polymorphic_within_group2 (same one amino acid in group 1 (conserved), that don't exist in group 2 (therefore 'diverged'), where within group2, multiple amino acids exist (polymorphic))
		# Conserved_within_group1_shared_between_groups_polymorphic_within_group2 (same one amino acid in group 1 (conserved), that also exist in group 2 (therefore 'shared'), Multiple amino acids in group 1 (polymorphic), that some are also in group 2 (therefore 'shared'), where within group2, multiple amino acids exist (polymorphic))
		# 
		# Of interest, are the ones that diverge between the groups (or exist/missing only in the group of interest)...:
		# Conserved_within_groups_diverged_between_groups, Polymorphic_within_groups_diverged_between_groups, Polymorphic_within_group1_diverged_between_groups_conserved_within_group2, Conserved_within_group1_diverged_between_groups_polymorphic_within_group2, Exist_only_in_group1_conserved, Exist_only_in_group2_conserved
		#
		#
		# Define gap status for each group (missing, only one species, or exist (more than one species))
		group1_existance_status = ''
		group1_poly_cons_status = ''
		group1_all_diverged = ''
		if len(aa_freq_group1) > 0:
			sp_without_gaps = [aa_freq_group1[i] for i in aa_freq_group1 if i != '-']
			if sum(sp_without_gaps) == 1:
				group1_existance_status = 'exist_only_in_one_species'
			elif sum(sp_without_gaps) > 1:
				group1_existance_status = 'exist'
				if len(sp_without_gaps) == 1:
					group1_poly_cons_status = 'conserved'
				else:
					group1_poly_cons_status = 'polymorphic'
					if sum(sp_without_gaps) == len(sp_without_gaps):
						group1_all_diverged = 'Yes'
			else:
				group1_existance_status = 'missing'
		else:
			 group1_existance_status = 'missing'
		#
		#
		group2_existance_status = ''
		group2_poly_cons_status = ''
		group2_all_diverged = ''
		if len(aa_freq_group2) > 0:
			sp_without_gaps = [aa_freq_group2[i] for i in aa_freq_group2 if i != '-']
			if sum(sp_without_gaps) == 1:
				group2_existance_status = 'exist_only_in_one_species'
			elif sum(sp_without_gaps) > 1:
				group2_existance_status = 'exist'
				if len(sp_without_gaps) == 1:
					group2_poly_cons_status = 'conserved'
				else:
					group2_poly_cons_status = 'polymorphic'
					if sum(sp_without_gaps) == len(sp_without_gaps):
						group2_all_diverged = 'Yes'
			else:
				group2_existance_status = 'missing'
		else:
			 group2_existance_status = 'missing'
		#
		#
		# Compare the two groups characters and define if they share or completely diverged
		share_diverged_status = ''
		common_char = list(set(group1_aa_list) & set(group2_aa_list))
		common_non_gap_char = [j for j in common_char if j != '-']
		if len(common_non_gap_char) > 0:
			share_diverged_status = 'share'
		else:
			share_diverged_status = 'diverged'
		#
		#
		# Now based on combinations of the above flags, can describe the overall status of a site:
		site_between_group_status = ''
		if group1_existance_status == 'missing' and group2_existance_status == 'missing':
			site_between_group_status = 'All_gaps'
		if group1_existance_status == 'exist_only_in_one_species' and group2_existance_status == 'exist_only_in_one_species':
			site_between_group_status = 'Exist_only_in_one_species_per_group'
		if group1_existance_status == 'missing' and group2_existance_status == 'exist_only_in_one_species':
			site_between_group_status = 'Exist_only_in_one_species_of_' + str(args.group2)
		if group2_existance_status == 'missing' and group1_existance_status == 'exist_only_in_one_species':
			site_between_group_status = 'Exist_only_in_one_species_of_' + str(args.group1)
		if group1_existance_status == 'missing' and group2_existance_status == 'exist':
			if group2_poly_cons_status == 'conserved':
				site_between_group_status = 'Exist_only_in_'  + str(args.group2) + '_conserved'
			elif group2_poly_cons_status == 'polymorphic':
				site_between_group_status = 'Exist_only_in_' + str(args.group2) + ' _polymorphic'
		if group2_existance_status == 'missing' and group1_existance_status == 'exist':
			if group1_poly_cons_status == 'conserved':
				site_between_group_status = 'Exist_only_in_' + str(args.group1) + '_conserved'
			elif group1_poly_cons_status == 'polymorphic':
				site_between_group_status = 'Exist_only_in_' + str(args.group1) + '_polymorphic'
		if group1_existance_status == 'exist' and group2_existance_status == 'exist_only_in_one_species':
			if group1_poly_cons_status == 'conserved':
				if share_diverged_status == 'share':
					site_between_group_status = 'Conserved_within_' + str(args.group1) + '_shared_between_groups_Exist_only_in_one_species_of_'  + str(args.group2)
				else:
					site_between_group_status = 'Conserved_within_' + str(args.group1) + '_diverged_between_groups_Exist_only_in_one_species_of_' + str(args.group2)
			elif group1_poly_cons_status == 'polymorphic':
				if share_diverged_status == 'share':
					site_between_group_status = 'Polymorphic_within_' + str(args.group1) + '_shared_between_groups_Exist_only_in_one_species_of_' + str(args.group2)
				else:
					site_between_group_status = 'Polymorphic_within_' + str(args.group1) + '_diverged_between_groups_Exist_only_in_one_species_of_' + str(args.group2)
		if group2_existance_status == 'exist' and group1_existance_status == 'exist_only_in_one_species':
			if group2_poly_cons_status == 'conserved':
				if share_diverged_status == 'share':
					site_between_group_status = 'Exist_only_in_one_species_of_' + str(args.group1) + '_shared_between_groups_Conserved_within_' + str(args.group2)
				else:
					site_between_group_status = 'Exist_only_in_one_species_of_' + str(args.group1) + '_diverged_between_groups_Conserved_within_' + str(args.group2)
			elif group2_poly_cons_status == 'polymorphic':
				if share_diverged_status == 'share':
					site_between_group_status = 'Exist_only_in_one_species_of_' + str(args.group1) + '_shared_between_groups_Polymorphic_within_' + str(args.group2)
				else:
					site_between_group_status = 'Exist_only_in_one_species_of_' + str(args.group1) + '_diverged_between_groups_Polymorphic_within_' + str(args.group2)
		if group1_existance_status == 'exist' and group2_existance_status == 'exist':
			if group1_poly_cons_status == 'conserved' and group2_poly_cons_status == 'conserved':
				if share_diverged_status == 'share':
					if '-' in list(set(group1_aa_list + group2_aa_list)):	# not all species have sequence
						site_between_group_status = 'Conserved_within_groups_and_between_groups'
					else:
						site_between_group_status = 'Conserved_in_all_species'
				else:
					site_between_group_status = 'Conserved_within_groups_diverged_between_groups'
			if group1_poly_cons_status == 'polymorphic' and group2_poly_cons_status == 'polymorphic':
				if share_diverged_status == 'share':
					site_between_group_status = 'Polymorphic_within_groups_shared_between_groups'
				else:
					if group1_all_diverged == 'Yes' and group2_all_diverged == 'Yes':
						site_between_group_status = 'Diverged_in_all_species'
					else:
						site_between_group_status = 'Polymorphic_within_groups_diverged_between_groups'
			if group1_poly_cons_status == 'conserved' and group2_poly_cons_status == 'polymorphic':
				if share_diverged_status == 'share':			
					site_between_group_status = 'Conserved_within_' + str(args.group1) + '_shared_between_groups_Polymorphic_within_' + str(args.group2)
				else:
					site_between_group_status = 'Conserved_within_' + str(args.group1) + '_diverged_between_groups_Polymorphic_within_' + str(args.group2)
			if group1_poly_cons_status == 'polymorphic' and group2_poly_cons_status == 'conserved':
				if share_diverged_status == 'share':			
					site_between_group_status = 'Polymorphic_within_' + str(args.group1) + '_shared_between_groups_Conserved_within_' + str(args.group2)
				else:
					site_between_group_status = 'Polymorphic_within_' + str(args.group1) + '_diverged_between_groups_Conserved_within_' + str(args.group2)
		#
		#
		## Get unique aa properties of each species group:
		if aa_grouping_method == 'Katzir':
			aa_prop_group1 = ';'.join(list(set(group1_aa_properties_ERGPMM_book_2006_7_groups)))
			aa_prop_group2 = ';'.join(list(set(group2_aa_properties_ERGPMM_book_2006_7_groups)))			
			aa_prop_uniq_group1 = list(set(group1_aa_properties_ERGPMM_book_2006_7_groups) - set(group2_aa_properties_ERGPMM_book_2006_7_groups))
			if len(aa_prop_uniq_group1) == 0:
				aa_prop_uniq_group1 = 'NA'
			else:
				aa_prop_uniq_group1 = ';'.join(aa_prop_uniq_group1)
			aa_prop_uniq_group2 = list(set(group2_aa_properties_ERGPMM_book_2006_7_groups) - set(group1_aa_properties_ERGPMM_book_2006_7_groups))
			if len(aa_prop_uniq_group2) == 0:
				aa_prop_uniq_group2 = 'NA'
			else:
				aa_prop_uniq_group2 = ';'.join(aa_prop_uniq_group2)
		elif aa_grouping_method == 'Murphy':
			aa_prop_group1 = ';'.join(list(set(group1_aa_properties_Murphy_2000_10_groups)))
			aa_prop_group2 = ';'.join(list(set(group2_aa_properties_Murphy_2000_10_groups)))
			aa_prop_uniq_group1 = list(set(group1_aa_properties_Murphy_2000_10_groups) - set(group2_aa_properties_Murphy_2000_10_groups))
			if len(aa_prop_uniq_group1) == 0:
				aa_prop_uniq_group1 = 'NA'
			else:
				aa_prop_uniq_group1 = ';'.join(aa_prop_uniq_group1)
			aa_prop_uniq_group2 = list(set(group2_aa_properties_Murphy_2000_10_groups) - set(group1_aa_properties_Murphy_2000_10_groups))
			if len(aa_prop_uniq_group2) == 0:
				aa_prop_uniq_group2 = 'NA'
			else:
				aa_prop_uniq_group2 = ';'.join(aa_prop_uniq_group2)
		else:
			aa_prop_group1 = ';'.join(list(set(group1_aa_properties_Taylor)))
			aa_prop_group2 = ';'.join(list(set(group2_aa_properties_Taylor)))
			aa_prop_uniq_group1 = list(set(group1_aa_properties_Taylor) - set(group2_aa_properties_Taylor))
			if len(aa_prop_uniq_group1) == 0:
				aa_prop_uniq_group1 = 'NA'
			else:
				aa_prop_uniq_group1 = ';'.join(aa_prop_uniq_group1)
			aa_prop_uniq_group2 = list(set(group2_aa_properties_Taylor) - set(group1_aa_properties_Taylor))
			if len(aa_prop_uniq_group2) == 0:
				aa_prop_uniq_group2 = 'NA'
			else:
				aa_prop_uniq_group2 = ';'.join(aa_prop_uniq_group2)		
		#
		## Before, adding the info to the output, replace empty variables with 'NA': (can happen if for instance there are only members of one group in the alignment...)
		if len(group1_sp_list) == 0:
			group1_sp_list = ['NA']
		if len(group2_sp_list) == 0:
			group2_sp_list = ['NA']
		if len(group1_aa_list_with_counts) == 0:
			group1_aa_list_with_counts = ['NA']
		if len(group2_aa_list_with_counts) == 0:
			group2_aa_list_with_counts = ['NA']
		if len(aa_prop_group1) == 0:
			aa_prop_group1 = 'NA'
		if len(aa_prop_group2) == 0:
			aa_prop_group2 = 'NA'
		#
		## Add all the info to the outputs:
		all_aa_res.append('\t'.join([str(c+1), site_between_group_status, ';'.join(group1_aa_list_with_counts), ';'.join(group2_aa_list_with_counts),aa_prop_group1, aa_prop_group2, aa_prop_uniq_group1, aa_prop_uniq_group2]))
		all_aa_all_res.append('\t'.join([str(c+1), site_between_group_status, ';'.join(group1_sp_list), ';'.join(group1_aa_list_with_counts), ';'.join(group2_sp_list), ';'.join(group2_aa_list_with_counts), aa_prop_group1, aa_prop_group2, aa_prop_uniq_group1, aa_prop_uniq_group2, str(group1_median_aa_Hydrophobicity), str(group2_median_aa_Hydrophobicity), str(sum_of_pairs_score_per_site), str(sum_of_pairs_score_per_site_group1), str(sum_of_pairs_score_per_site_group2), str(sum_of_pairs_score_per_site_between_groups), str(norm_sum_of_pairs_score_per_site), str(norm_sum_of_pairs_score_per_site_group1), str(norm_sum_of_pairs_score_per_site_group2), str(norm_sum_of_pairs_score_per_site_between_groups)]))
		# Check relationship status, and save the ones of interest in a separate file:
		if site_between_group_status == 'Conserved_within_groups_diverged_between_groups' or site_between_group_status == 'Polymorphic_within_groups_diverged_between_groups' or site_between_group_status == 'Polymorphic_within_' + str(args.group1) + '_diverged_between_groups_conserved_within_' + str(args.group2) or site_between_group_status == 'Conserved_within_' + str(args.group1) + '_diverged_between_groups_polymorphic_within_' + str(args.group2) or site_between_group_status == 'Exist_only_in_group1_conserved' or site_between_group_status == 'Exist_only_in_' + str(args.group2) + '_conserved':
			criteria_satisfied_aa_res.append('\t'.join([str(c+1), site_between_group_status, ';'.join(group1_aa_list_with_counts), ';'.join(group2_aa_list_with_counts), aa_prop_group1, aa_prop_group2, aa_prop_uniq_group1, aa_prop_uniq_group2]))
		# Add to all_aa_grouping output:
		all_aa_grouping_res.append('\t'.join([str(c+1), site_between_group_status,';'.join(group1_sp_list), ';'.join(group1_aa_list_with_counts), ';'.join(group2_sp_list), ';'.join(group2_aa_list_with_counts), group1_Katzir,group2_Katzir,group1_group2_common_aa_groups_Katzir,group1_Katzir_radical,group2_Katzir_radical,between_groups_radiacal_replacement_Katzir,	group1_Murphy,group2_Murphy,group1_group2_common_aa_groups_Murphy,group1_Murphy_radical,group2_Murphy_radical,between_groups_radiacal_replacement_Murphy, group1_Zhang_charge,group2_Zhang_charge,group1_group2_common_aa_groups_Zhang_charge,group1_Zhang_charge_radical,group2_Zhang_charge_radical,between_groups_radiacal_replacement_Zhang_charge,	group1_Zhang_polarity,group2_Zhang_polarity,group1_group2_common_aa_groups_Zhang_polarity,group1_Zhang_polarity_radical,group2_Zhang_polarity_radical,between_groups_radiacal_replacement_Zhang_polarity, group1_Zhang_polarity_and_volumn,group2_Zhang_polarity_and_volumn,group1_group2_common_aa_groups_Zhang_polarity_and_volumn,group1_Zhang_polarity_and_volumn_radical,group2_Zhang_polarity_and_volumn_radical,between_groups_radiacal_replacement_Zhang_polarity_and_volumn, group1_Betts_and_Russell_less_exchangeable_aa,group2_Betts_and_Russell_less_exchangeable_aa,group1_group2_common_aa_groups_Betts_and_Russell,group1_Betts_and_Russell_less_exchangeable_aa_radical,group2_Betts_and_Russell_less_exchangeable_aa_radical,between_groups_radiacal_replacement_Betts_and_Russell, str(within_group1_Grantham_1974_physicochemical_max_distance), str(within_group2_Grantham_1974_physicochemical_max_distance), str(between_group1_group2_Grantham_1974_physicochemical_max_distance),group1_Grantham_radical_aa_replacement,group2_Grantham_radical_aa_replacement,between_groups_radiacal_replacement_Katzir]))
		# For Grantham_1974_aa_radical_aa_replacement, considers d<100 as a conservative amino acid replacement, and d>=100 as radical amino acid replacement.		
		#
	#
	#
	#
	# Add a readme file with an explanation about the output:
	#readme = ''
	#now = datetime.datetime.now()	# get the current time and date (using the 'datetime' module...)
	#date_and_time_of_analysis = now.strftime("%Y-%m-%d %H:%M")	# get it in a simple format: e.g. '2019-07-17 23:56'
	# continue to build the readme....
	#
	## Get sorted position list for plotting:
	aa_positions = list(sum_of_pairs_score_per_site_dict.keys())	# get list of aa positions
	aa_positions.sort(key=int)	# order numerically...
	#
	## Get the various lists for plotting:
	sum_of_pairs_score_per_site_list = []
	sum_of_pairs_score_per_site_group1_list = []
	sum_of_pairs_score_per_site_group2_list = []
	sum_of_pairs_score_per_site_between_groups_list = []
	norm_sum_of_pairs_score_per_site_list = []
	norm_sum_of_pairs_score_per_site_group1_list = []
	norm_sum_of_pairs_score_per_site_group2_list = []
	norm_sum_of_pairs_score_per_site_between_groups_list = []
	Hydrophobicity_per_position_median_list = []
	Hydrophobicity_per_position_average_list = []
	Hydrophobicity_per_position_group1_median_list = []
	Hydrophobicity_per_position_group1_average_list = []
	Hydrophobicity_per_position_group2_median_list = []
	Hydrophobicity_per_position_group2_average_list = []
	#
	for p in aa_positions:
		if len(sum_of_pairs_score_per_site_dict) != 0:
			sum_of_pairs_score_per_site_list.append(sum_of_pairs_score_per_site_dict[p])
		if len(sum_of_pairs_score_per_site_dict_group1) != 0:
			sum_of_pairs_score_per_site_group1_list.append(sum_of_pairs_score_per_site_dict_group1[p])
		if len(sum_of_pairs_score_per_site_dict_group2) != 0:
			sum_of_pairs_score_per_site_group2_list.append(sum_of_pairs_score_per_site_dict_group2[p])
		if len(sum_of_pairs_score_per_site_dict_between_groups) != 0:
			sum_of_pairs_score_per_site_between_groups_list.append(sum_of_pairs_score_per_site_dict_between_groups[p])
		if len(norm_sum_of_pairs_score_per_site_dict) != 0:
			norm_sum_of_pairs_score_per_site_list.append(norm_sum_of_pairs_score_per_site_dict[p])
		if len(norm_sum_of_pairs_score_per_site_dict_group1) != 0:
			norm_sum_of_pairs_score_per_site_group1_list.append(norm_sum_of_pairs_score_per_site_dict_group1[p])
		if len(norm_sum_of_pairs_score_per_site_dict_group2) != 0:
			norm_sum_of_pairs_score_per_site_group2_list.append(norm_sum_of_pairs_score_per_site_dict_group2[p])
		if len(norm_sum_of_pairs_score_per_site_dict_between_groups) != 0:
			norm_sum_of_pairs_score_per_site_between_groups_list.append(norm_sum_of_pairs_score_per_site_dict_between_groups[p])	
		if len(Hydrophobicity_per_position_dict) != 0:
			Hydrophobicity_per_position_median_list.append(Hydrophobicity_per_position_dict[p]['Median'])
			Hydrophobicity_per_position_average_list.append(Hydrophobicity_per_position_dict[p]['Mean'])
		if len(Hydrophobicity_per_position_per_gr_dict_group1) != 0:
			Hydrophobicity_per_position_group1_median_list.append(Hydrophobicity_per_position_per_gr_dict_group1[p]['Median'])
			Hydrophobicity_per_position_group1_average_list.append(Hydrophobicity_per_position_per_gr_dict_group1[p]['Mean'])
		if len(Hydrophobicity_per_position_per_gr_dict_group2) != 0:
			Hydrophobicity_per_position_group2_median_list.append(Hydrophobicity_per_position_per_gr_dict_group2[p]['Median'])
			Hydrophobicity_per_position_group2_average_list.append(Hydrophobicity_per_position_per_gr_dict_group2[p]['Mean'])
	#
	#
	# Calculate the sliding window average and median for Hydrophobicity and sum of pairs:
	aa_positions_int = [int(i) for i in aa_positions]
	sw_positions = [min(i) for i in list(window(aa_positions_int,n=int(sw_size)))]
	sw_positions_ranges = [str(min(i)) + ':' + str(max(i)) for i in list(window(aa_positions_int,n=int(sw_size)))]
	sum_of_pairs_score_per_site_sw_median = med_wind(list(window(sum_of_pairs_score_per_site_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_sw_averages = ave_wind(list(window(sum_of_pairs_score_per_site_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_group1_sw_median = med_wind(list(window(sum_of_pairs_score_per_site_group1_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_group1_sw_averages = ave_wind(list(window(sum_of_pairs_score_per_site_group1_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_group2_sw_median = med_wind(list(window(sum_of_pairs_score_per_site_group2_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_group2_sw_averages = ave_wind(list(window(sum_of_pairs_score_per_site_group2_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_between_groups_sw_median = med_wind(list(window(sum_of_pairs_score_per_site_between_groups_list,n=int(sw_size))))
	sum_of_pairs_score_per_site_between_groups_sw_averages = ave_wind(list(window(sum_of_pairs_score_per_site_between_groups_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_sw_median = med_wind(list(window(norm_sum_of_pairs_score_per_site_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_sw_averages = ave_wind(list(window(norm_sum_of_pairs_score_per_site_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_group1_sw_median = med_wind(list(window(norm_sum_of_pairs_score_per_site_group1_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_group1_sw_averages = ave_wind(list(window(norm_sum_of_pairs_score_per_site_group1_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_group2_sw_median = med_wind(list(window(norm_sum_of_pairs_score_per_site_group2_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_group2_sw_averages = ave_wind(list(window(norm_sum_of_pairs_score_per_site_group2_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_between_groups_sw_median = med_wind(list(window(norm_sum_of_pairs_score_per_site_between_groups_list,n=int(sw_size))))
	norm_sum_of_pairs_score_per_site_between_groups_sw_averages = ave_wind(list(window(norm_sum_of_pairs_score_per_site_between_groups_list,n=int(sw_size))))
	Hydrophobicity_sw_median = med_wind(list(window(Hydrophobicity_per_position_median_list,n=int(sw_size))))
	Hydrophobicity_sw_averages = ave_wind(list(window(Hydrophobicity_per_position_average_list,n=int(sw_size))))
	Hydrophobicity_group1_sw_median = med_wind(list(window(Hydrophobicity_per_position_group1_median_list,n=int(sw_size))))
	Hydrophobicity_group1_sw_averages = ave_wind(list(window(Hydrophobicity_per_position_group1_average_list,n=int(sw_size))))
	Hydrophobicity_group2_sw_median = med_wind(list(window(Hydrophobicity_per_position_group2_median_list,n=int(sw_size))))
	Hydrophobicity_group2_sw_averages = ave_wind(list(window(Hydrophobicity_per_position_group2_average_list,n=int(sw_size))))
	#
	#
	## Generate an output for the sliding window values:
	sw_out = ['\t'.join(['Sliding_window' , 'All_species_sliding_window_median_' + hydropho_scale, 'All_species_sliding_window_average_' + hydropho_scale, str(args.group1) + '_sliding_window_median_' + hydropho_scale, str(args.group1) + '_sliding_window_average_' + hydropho_scale, str(args.group2) + '_sliding_window_median_' + hydropho_scale, str(args.group2) + '_sliding_window_average_' + hydropho_scale, 'Total_sliding_window_median_sum_of_pairs_score_per_site', 'Total_sliding_window_average_sum_of_pairs_score_per_site', str(args.group1) + '_within-group_sliding_window_median_sum_of_pairs_score_per_site', str(args.group1) + '_within-group_sliding_window_average_sum_of_pairs_score_per_site', str(args.group2) + '_within-group_sliding_window_median_sum_of_pairs_score_per_site', str(args.group2) + '_within-group_sliding_window_average_sum_of_pairs_score_per_site', 'Between-groups_sliding_window_median_sum_of_pairs_score_per_site', 'Between-groups_sliding_window_average_sum_of_pairs_score_per_site', 'Total_sliding_window_median_normalized_sum_of_pairs_score_per_site', 'Total_sliding_window_average_normalized_sum_of_pairs_score_per_site', str(args.group1) + '_within-group_sliding_window_median_normalized_sum_of_pairs_score_per_site', str(args.group1) + '_within-group_sliding_window_average_normalized_sum_of_pairs_score_per_site', str(args.group2) + '_within-group_sliding_window_median_normalized_sum_of_pairs_score_per_site', str(args.group2) + '_within-group_sliding_window_average_normalized_sum_of_pairs_score_per_site', 'Between-groups_sliding_window_median_normalized_sum_of_pairs_score_per_site', 'Between-groups_sliding_window_average_normalized_sum_of_pairs_score_per_site'])]
	for s in range(len(sw_positions)):
		sw_out.append('\t'.join([str(sw_positions_ranges[s]),str(Hydrophobicity_sw_median[s]),str(Hydrophobicity_sw_averages[s]),str(Hydrophobicity_group1_sw_median[s]),str(Hydrophobicity_group1_sw_averages[s]),str(Hydrophobicity_group2_sw_median[s]),str(Hydrophobicity_group2_sw_averages[s]),str(sum_of_pairs_score_per_site_sw_median[s]),str(sum_of_pairs_score_per_site_sw_averages[s]),str(sum_of_pairs_score_per_site_group1_sw_median[s]),str(sum_of_pairs_score_per_site_group1_sw_averages[s]),str(sum_of_pairs_score_per_site_group2_sw_median[s]),str(sum_of_pairs_score_per_site_group2_sw_averages[s]),str(sum_of_pairs_score_per_site_between_groups_sw_median[s]),str(sum_of_pairs_score_per_site_between_groups_sw_averages[s]),str(norm_sum_of_pairs_score_per_site_sw_median[s]),str(norm_sum_of_pairs_score_per_site_sw_averages[s]),str(norm_sum_of_pairs_score_per_site_group1_sw_median[s]),str(norm_sum_of_pairs_score_per_site_group1_sw_averages[s]),str(norm_sum_of_pairs_score_per_site_group2_sw_median[s]),str(norm_sum_of_pairs_score_per_site_group2_sw_averages[s]),str(norm_sum_of_pairs_score_per_site_between_groups_sw_median[s]),str(norm_sum_of_pairs_score_per_site_between_groups_sw_averages[s])]))
	#
	#
	## If the --plot flag was indicated, plot the figures:
	# Plot the different scores (within, between, within_group1, within_group2), and think how to detect diverging sites...    play with the sum of pairs between and within groups? See which function make sense (within_group1 + within_group2 > between_group1_group2)
	# Plot values along the alignment positions (Needs that matplotlib will be installed and importable - better to use matplotlib v2.0.2 that works):
	if str(args.plot) == 'True':
		# sum_of_pairs_score:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(aa_positions_int,sum_of_pairs_score_per_site_list, 'm^', label = 'Total sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(aa_positions_int,sum_of_pairs_score_per_site_group1_list, 'co', label = str(args.group1) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(aa_positions_int,sum_of_pairs_score_per_site_group2_list, 'gs', label = str(args.group2) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(aa_positions_int,sum_of_pairs_score_per_site_between_groups_list, 'yv', label = 'Between-groups sum of pairs score per site', ms = 2,alpha=0.90)
			plt.ylabel('Sum of pairs score', fontsize=10)
			plt.xlabel('Alignment position', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sum_of_pairs_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sum of pairs plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		# norm_sum_of_pairs_score:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(aa_positions_int,norm_sum_of_pairs_score_per_site_list, 'm^', label = 'Total normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(aa_positions_int,norm_sum_of_pairs_score_per_site_group1_list, 'co', label = str(args.group1) + ' within-group normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(aa_positions_int,norm_sum_of_pairs_score_per_site_group2_list, 'gs', label = str(args.group2) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(aa_positions_int,norm_sum_of_pairs_score_per_site_between_groups_list, 'yv', label = 'Between-groups normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.ylabel('Normalized sum of pairs score', fontsize=10)
			plt.xlabel('Alignment position', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_norm_sum_of_pairs_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate normalized sum of pairs plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		# Hydrophobicity:
		try:
			import matplotlib.pyplot as plt
			fig_size = [8, 6]
			fig, axs = plt.subplots(len(Hydrophobicity_per_position_per_sp_dict),1,sharex=True,sharey=True,figsize=fig_size)
			plt.style.use('ggplot')
			fig.text(0.04, 0.5, hydropho_scale.replace('_',' ') + ' level', va='center',rotation='vertical')
			fig.text(0.5, 0.04, 'Alignment position', ha='center')
			count = 0
			for s in Hydrophobicity_per_position_per_sp_dict:
				sp_name = s
				if s in species_label_to_name.keys():
					sp_name = species_label_to_name[s]
				else:
					sp_name = s
				aa_pos = []
				hydropat = []
				for a in Hydrophobicity_per_position_per_sp_dict[s]:
					if Hydrophobicity_per_position_per_sp_dict[s][a] != 'NA':
						aa_pos.append(a)
						hydropat.append(Hydrophobicity_per_position_per_sp_dict[s][a])
				if s in group_1:
					axs[count].plot(aa_pos, hydropat, linestyle='None', marker='o', color='b' , label = sp_name, ms = 2)
					axs[count].legend(loc='center left', bbox_to_anchor=(1, 0.5))
				elif s in group_2:
					axs[count].plot(aa_pos, hydropat, linestyle='None', marker='o', color='g' , label = sp_name, ms = 2)
					axs[count].legend(loc='center left', bbox_to_anchor=(1, 0.5))
				else:
					axs[count].plot(aa_pos, hydropat, linestyle='None', marker='o', color='r', label = sp_name, ms = 2)
					axs[count].legend(loc='center left', bbox_to_anchor=(1, 0.5))
				count = count + 1
			#
			axs[0].set_title('Hydrophobicity plot for ' + str(args.fasta_msa_file).rsplit('.',1)[0],fontsize=12)
			fig.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_Hydrophobicity_plot.pdf", bbox_inches='tight')	
			fig.clear()
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate Hydrophobicity plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		#
		#
		## Plot sliding window sum of pairs (average or median) - Needs that matplotlib will be installed and importable - better to use matplotlib v2.0.2 that works:
		# Sliding window median sum_of_pairs_score:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(sw_positions,sum_of_pairs_score_per_site_sw_median, 'm^', label = 'Total sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,sum_of_pairs_score_per_site_group1_sw_median, 'co', label = str(args.group1) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,sum_of_pairs_score_per_site_group2_sw_median, 'gs', label = str(args.group2) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,sum_of_pairs_score_per_site_between_groups_sw_median, 'yv', label = 'Between-groups sum of pairs score per site', ms = 2,alpha=0.90)
			plt.ylabel('Sliding window median sum of pairs score', fontsize=10)
			plt.xlabel('Sliding window', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sliding_window_median_sum_of_pairs_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sliding window median sum of pairs plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		#
		# Sliding window average sum_of_pairs_score:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(sw_positions,sum_of_pairs_score_per_site_sw_averages, 'm^', label = 'Total sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,sum_of_pairs_score_per_site_group1_sw_averages, 'co', label = str(args.group1) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,sum_of_pairs_score_per_site_group2_sw_averages, 'gs', label = str(args.group2) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,sum_of_pairs_score_per_site_between_groups_sw_averages, 'yv', label = 'Between-groups sum of pairs score per site', ms = 2,alpha=0.90)
			plt.ylabel('Sliding window average sum of pairs score', fontsize=10)
			plt.xlabel('Sliding window', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sliding_window_average_sum_of_pairs_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sliding window average sum of pairs plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		#
		# Sliding window median norm_sum_of_pairs_score:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_sw_median, 'm^', label = 'Total normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_group1_sw_median, 'co', label = str(args.group1) + ' within-group normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_group2_sw_median, 'gs', label = str(args.group2) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_between_groups_sw_median, 'yv', label = 'Between-groups normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.ylabel('Sliding window median normalized sum of pairs score', fontsize=10)
			plt.xlabel('Sliding window', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sliding_window_median_norm_sum_of_pairs_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sliding window median normalized sum of pairs plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		#
		# Sliding window average norm_sum_of_pairs_score:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_sw_averages, 'm^', label = 'Total normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_group1_sw_averages, 'co', label = str(args.group1) + ' within-group normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_group2_sw_averages, 'gs', label = str(args.group2) + ' within-group sum of pairs score per site', ms = 2,alpha=0.90)
			plt.plot(sw_positions,norm_sum_of_pairs_score_per_site_between_groups_sw_averages, 'yv', label = 'Between-groups normalized sum of pairs score per site', ms = 2,alpha=0.90)
			plt.ylabel('Sliding window average normalized sum of pairs score', fontsize=10)
			plt.xlabel('Sliding window', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sliding_window_average_norm_sum_of_pairs_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sliding window average normalized sum of pairs plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		#
		#
		## Plot sliding window Hydrophobicity (Here not plotting per species sliding window, but per group - can add per species later...):
		# Median Hydrophobicity:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(sw_positions,Hydrophobicity_sw_median, 'm^', label = 'All species', ms = 2,alpha=0.90)
			plt.plot(sw_positions,Hydrophobicity_group1_sw_median, 'co', label = str(args.group1), ms = 2,alpha=0.90)
			plt.plot(sw_positions,Hydrophobicity_group2_sw_median, 'gs', label = str(args.group2), ms = 2,alpha=0.90)
			plt.ylabel('Median sliding window ' + hydropho_scale.replace('_',' '), fontsize=10)
			plt.xlabel('Sliding window', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sliding_window_median_hydrophobicity_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sliding window median hydrophobicity plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
		#
		#
		# Average Hydrophobicity:
		try:
			import matplotlib.pyplot as plt
			f = plt.figure()
			plt.style.use('ggplot')
			plt.plot(sw_positions,Hydrophobicity_sw_averages, 'm^', label = 'All species', ms = 2,alpha=0.90)
			plt.plot(sw_positions,Hydrophobicity_group1_sw_averages, 'co', label = str(args.group1), ms = 2,alpha=0.90)
			plt.plot(sw_positions,Hydrophobicity_group2_sw_averages, 'gs', label = str(args.group2), ms = 2,alpha=0.90)
			plt.ylabel('Average sliding window ' + hydropho_scale.replace('_',' '), fontsize=10)
			plt.xlabel('Sliding window', fontsize=10)
			plt.legend(fontsize=8)
			f.savefig(str(args.fasta_msa_file).rsplit('.',1)[0] + "_sliding_window_average_hydrophobicity_" + sub_matrix + "_" + str(gap_pen) + ".pdf", bbox_inches='tight')
		except:
			print("Couldn't import matplotlib, therefore, couldn't generate sliding window average hydrophobicity plot figures. If you wish output plots, please install matplotlib first: e.g. by typing in your terminal 'pip install matplotlib==2.0.2'. Another posibility - you may want to try exacuting the script using 'pythonw tgaa.py' instead of 'python tgaa.py'")
	#
	#
	#
	## Save results:
	with open(str(args.output), 'w') as outFile:
		outFile.write('\n'.join(criteria_satisfied_aa_res) + '\n')
	with open(str(args.output).rsplit('.',1)[0] + '_all_aa.' + str(args.output).rsplit('.',1)[1], 'w') as outFile:
		outFile.write('\n'.join(all_aa_res) + '\n')
	with open(str(args.output).rsplit('.',1)[0] + '_all_aa_all_info.' + str(args.output).rsplit('.',1)[1], 'w') as outFile:
		outFile.write('\n'.join(all_aa_all_res) + '\n')
	with open(str(args.output).rsplit('.',1)[0] + '_sliding_window.' + str(args.output).rsplit('.',1)[1], 'w') as outFile:
		outFile.write('\n'.join(sw_out) + '\n')
	if str(args.output_all_aa_grouping_file) == 'True':	# if the flag was called, need to output all the aa-grouping comparisons to an additional file...
		with open(str(args.output).rsplit('.',1)[0] + '_all_aa_grouping_comp.' + str(args.output).rsplit('.',1)[1], 'w') as outFile:
			outFile.write('\n'.join(all_aa_grouping_res) + '\n')

if __name__ == '__main__':
	main()

