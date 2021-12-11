'''
VERSION
	0.1
	Python3.8
AUTHOR
	Fuentes Mendez David Gregorio
DESCRIPTION
	Simple SNP management module
CATEGORY
	Package
USAGE
	Helps in the analysis of SNP variants
	Has useful methods that can tell you:
		-if your snp is inside your gbk
		-
ARGUMNETS

SEE ALSO
	My repository in GitHub: https://github.com/DavidGreeg
NOTE
	I thought it would be nice to do my first program comlletely in english
	Coming soon, taco version, I mean, spanish version.	
'''


import re
import Bio
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *


#====================================ARGUMENTS========================================#
#-----------------------------------------------------------------------ARGUMENT PARSER
# Start ArgumentParser()
args = argparse.ArgumentParser(description = "Simple SNP management program, intended \
							for helping in the analysis of other files where we might \
							wonder about SNPs location or things like presence in the \
							coding sequence, or even what would frameshift / missense \
							mutations look like there")

# ----------------------------------------------------------------------------ARGUMENTS
# Define the arguments
args.add_argument("-s","--snp_file",
	type = str,
	help = "Path to csv file containing positions, alleles and traits of the SNPs, in \
				that order so the appropiate SNP class objects can be made",
	required = False)
args.add_argument("-p","--pos_file",
	type = str,
	help = "Path of txt file containing the positions of each SNP you'll use",
	required = False)
args.add_argument("-a","--allele_file",
	type = str,
	help = "Path of txt file containing the alleles of each SNP you want to use",
	required = False)
args.add_argument("-t","--traits_file",
	type = str,
	help = "Path of txt file containing the traits associated to each SNP you'll use",
	required = False)
args.add_argument("-g","--gbk_file",
	type = str,
	help = "Path to gbk file containing a gene or genes that where we may check SNPs",
	required = False)

argx = args.parse_args() 
#-------------------------------------------------------------------ARGUMENT PROCESSING
# Check if there is no arguments or illogical options selected


no_args = True
if not (argx.snp_file or argx.pos_file or argx.allele_file or argx.traits_file):
	print('No arguments were selected.')
	if (argx.snp_file) and (argx.pos_file or argx.allele_file or argx.traits_file):
		print("Illogical argument selection. Both -s and either -p or -a or -t, can't \
		be selected al the same time")
else:
	no_args = False
	if argx.snp_file:
		s_file = argx.snp_file
	else:
		p_file = argx.pos_file
		a_file = argx.allele_file
		t_file = argx.traits_file

#-------------------------------------------------------------------FUNCTION DEFINITION
# Define the funcitons this module will contain


pra