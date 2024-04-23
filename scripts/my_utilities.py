'''
VERSION
	0.1
	Python3.8
AUTHOR
	Fuentes Mendez David Gregorio
DESCRIPTION
	This module comprises some functions I thougth might
	be useful later, or moved here from some code outside
	because I wanted to modularize that part.
CATEGORY
	Package/Module
USAGE
	Move this file to the same directory as the python
	file where you'd like to use it. Then you call it
	with:		
		´import my_utilities as my_ut´
	
	And call any of the functions below as:
		´my_ut.funtion(arg1, arg2, ...)´
ARGUMENTS
	[Will later be added properly]
SEE ALSO
	Check out my repository on Github:
		https://github.com/DavidGreeg
'''

def check_ext(file:str, ext:str):
	ext_dict = {'Text':['txt'],'Word':['doc','docm','docx'],'PDF':['pdf'], 'HTML':['html'],'GenBank':['gb','gbk'],'FASTA':['fasta','fas','fna','faa','fa'],'CSV':['csv'], 'TSV':['tsv']}
	
	if ext in ext_dict.keys():
		str_splt = file.split(".")
		str_extn = str_splt[-1] 
		if (len(str_extn) > 1) and (str_extn in ext_dict[ext]):
			return True
		else:
			print('File extension (%s) and given option (%s) are not the same in %s' %(str_extn, ext, file))
			return False
	else:
		print('Invalid file extension')
		return False

def get_lastword(string):
	ls = string.split()
	string = ls[-1]
	return string

def hgdfhgdf(string):
	string = string.lower()
	types=['dna','rna','prot','protein']
	if string in types:
		return True
	else:
		print('Invalid sequence type')
		return False
