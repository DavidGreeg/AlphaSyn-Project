'''
VERSION
	0.1
	Python3.8
AUTHOR
	Fuentes Mendez David Gregorio
DESCRIPTION
	This module comprises some functions that might come in handy for
	new users who may just want to do a quick look to their gbk files.
	Useful 
CATEGORY
	Package/Module
USAGE
	
ARGUMENTS
SEE ALSO
'''




from Bio import SeqIO

for rec in SeqIO.parse("SNCA_ReferenceGRCh38_HomoSapiens.gb", "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "CDS":
                print("Protein: %s; ID: %s;" %(feature.qualifiers["product"][0],feature.qualifiers["protein_id"][0]))
                print(feature.location)
                i=0
                for loc in feature.location.parts:
                    i+=1
                    print("Location exon %s: %s:%s" %(i,loc.start,loc.end))
                print(feature.location.extract(rec).seq)
                print("------------------------------------------------------------")
#               print(feature.location.parts)
#               print(dir(feature.location.parts[0]))
#               print(str(feature.location.parts[0]))
#               print(feature.location.parts[3].start)
