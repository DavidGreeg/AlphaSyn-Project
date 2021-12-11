#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Al principio vamos a abrir una terminal para realizar observaciones rápidas
get_ipython().run_line_magic('qtconsole', '')


# ### Pruebas de funciones

# In[2]:


from Bio.Seq import Seq

cds_ref = Seq("ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA")
cds_loc = Seq("ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA")

print(len(cds_ref),len(cds_loc))


# ### GenBank records

# In[37]:


#to: gbk_utilities.py
from Bio import SeqIO
import my_utilities as my_ut

file1="SNCA_ReferenceGRCh38_HomoSapiens.gb"
file2="SNCA_RefSeq_HomoSapiens.gb"

def report_cds(file):
    if my_ut.check_ext(file,"GenBank"):
        for rec in SeqIO.parse(file, "genbank"):
            if rec.features:
                cds=False
                for feature in rec.features:
                    if feature.type == "CDS":
                        cds=True
                        print("CDS Product: %s; ID: %s;" %(feature.qualifiers["product"][0],feature.qualifiers["protein_id"][0]))
                        i=0
                        for loc in feature.location.parts:
                            i+=1
                            print("Location exon %s: %s:%s" %(i,loc.start,loc.end))
                        print("\nSequence DNA: \n%s"%(feature.location.extract(rec).seq))
                        print("\nSequence Prot: \n%s"%(feature.location.extract(rec).seq.translate()))                    
                        print("------------------------------------------------------------")
                if not cds:
                    print("No CDS sequences were found")
        return

#Couldn't make this s**t work, spent more than 2 hours trying to figure out why it doesn't let me import my 
#modules. Keeps printing nameError and I honestly can't figure out why
#
#It seems like the problem is some incompatibility between jupyter and module importation
def get_seqs_cds(file, seq_type:str, fasta=False, out_file=""):
    if my_ut.hgdfhgdf(seq_type):
        if my_ut.check_ext(file,"GenBank"):
            if fasta:
                if my_ut.check_ext(out_file,"FASTA"):
                    o_file=open(out_file, "w")

            for rec in SeqIO.parse(file, "genbank"):
                if rec.features:
                    cds=False
                    for feature in rec.features:

                        if feature.type == "CDS" and fasta:
                            cds=True
                            o_file.write(">%s %s" %(feature.qualifiers["product"][0],feature.qualifiers["protein_id"][0]))
                            if seq_type == "rna":
                                o_file.write("%s"%(feature.location.extract(rec).seq.transcribe()))
                            elif seq_type == ("prot" or "protein"):
                                o_file.write("%s"%(feature.location.extract(rec).seq.translate()))
                            else:
                                o_file.write("%s"%(feature.location.extract(rec).seq))

                        if feature.type == "CDS" and not fasta:
                            cds=True
                            print(">%s %s" %(feature.qualifiers["product"][0],feature.qualifiers["protein_id"][0]))
                            if seq_type == "rna":
                                print("%s"%(feature.location.extract(rec).seq.transcribe()))
                            elif seq_type == ("prot" or "protein"):
                                print("%s"%(feature.location.extract(rec).seq.translate()))
                            else:
                                print("%s" %(feature.location.extract(rec).seq))

                    if fasta:
                    	o_file.close()

                    if not cds:
                        print("No CDS sequences were found")
                    return
####################################################################################################################


report_cds(file1)


# In[5]:


#to: gbk_utilities.py
import re
import numpy as np
from Bio import SeqIO

file = "SNCA_ReferenceGRCh38_HomoSapiens.gb"

def location_ori(file):
#This function will return us the location of our ORIGIN sequence of this GenBank file.
#Will work properly if there's only one ORIGIN (SeqRecord) per file
    for rec in SeqIO.parse(file, "genbank"):
        if rec.annotations:
            if rec.annotations['accessions']:
                locs = rec.annotations['accessions'][2]
    index = [1,3]
    #This index vector is created because the following 're.split()' will output a list
    #of strings with the structure: ['complement' '89724099' '' '89838324' '']
    #We also cast it to np.int_ afterwards
    locs = np.array(re.split("[().]", locs))[index]
    locs = list(np.int_(locs))
    return locs

print(location_ori(file))
## NOTE: Make location_oris() later


# In[6]:


#to: gbk_utilities.py

file1 = "SNCA_ReferenceGRCh38_HomoSapiens.gb"

def get_cds_indx(file):
    if my_ut.check_ext(file,"GenBank"):        
        for rec in SeqIO.parse(file, "genbank"):            
            if rec.features:                
                cds=False
                cind_ls = list() #cds_index_list
                eind_ls = list() #exon_index_list
                i=0
                for feature in rec.features:                    
                    
                    if feature.type == "CDS":                        
                        cds=True                    
                        
                        #You might wanna look at the intervals each exon has with the omitted commands
                        #print("CDS Product: %s; ID: %s;" %(feature.qualifiers["product"][0],feature.qualifiers["protein_id"][0]))
                        
                        i=0
                        for loc in feature.location.parts:
                            i+=1
                            eind_ls = eind_ls + list(range(loc.start,loc.end))
                            
                            #print("Location exon %s: %s:%s" %(i,loc.start,loc.end))
                            #print("------------------------------------------------------------")
                        
                        cind_ls.append(eind_ls)                        
                if not cds:
                    print("No CDS sequences were found")
                    
        return cind_ls

print(get_cds_indx(file1))
###This funtion returns a list of lists, where each element is the location of each nucleotide in each
###exon in each cds-variant/isoform of the gene. This will come in handy in order to make logic vectors.


# In[42]:


#to: snp_utilities.py
import numpy as np
#import gbk_utilities as gbk_u

#Please make sure your SNPs are from the same chromosome as the GenBank file you're compairing them with
file_snp = "SNP_Locations.txt"
file_gbk = "SNCA_ReferenceGRCh38_HomoSapiens.gb"
start_gb = location_ori(file_gbk)[0]
#The location of our GenBank file inside the chromosome
#Given that 'start_gb' is the first base inside our CDS sequence, then we would like to get the following
#result when snp were equal to start_gb: (snp=start_gb; snp-start_gb=1)
#Then, and for the sake of legibility we're just gonna rest 1 to start_gb
start_gb = start_gb - 1

def snps_list_parse(file, start_gb:int=0):
    locs = open(file)
    locl = str(locs.read()).split()

    snps = np.int_(locl)    #Here we change our file datatype to 'int' as well as we make it a numpy array
                            #which will make it easier to get their relative positions
    snps = snps - start_gb
    snps = snps[snps >= 1]  #Here we select only the SNPs whose location is at or after the start position
                            #in our GenBank file
    snps = snps.tolist()    #Then we cast the array back to list (optional)
    return snps

print(snps_list_parse(file_snp, start_gb))
len(snps_list_parse(file_snp, start_gb))


# ### Check if SNPs are inside

# In[10]:


#to: snp_utilities.py 
import my_utilities as my_ut

file_gb = "SNCA_ReferenceGRCh38_HomoSapiens.gb"

file_snp = "SNP_Locations.txt"
start_gb = location_ori(file_gb)[0]
snp_ls = snps_list_parse(file_snp, start_gb)

def report_snp_in_gb(file,snps):
    if my_ut.check_ext(file, "GenBank"):
        for snp in snps:
            for rec in SeqIO.parse(file, "genbank"):
                if rec.features:
                    for feature in rec.features:
                        if feature.type == "CDS":
                            if snp in feature.location:
                                print("%s is in %s, ID: %s" %(snp,feature.qualifiers["product"][0],feature.qualifiers["protein_id"][0]))
    return


report_snp_in_gb(file_gb,snp_ls)


# In[43]:


import my_utilities as my_ut

file_gb = "SNCA_ReferenceGRCh38_HomoSapiens.gb"

file_snp = "SNP_Locations.txt"
start_gb = location_ori(file_gb)[0]
snp_ls = snps_list_parse(file_snp, start_gb)


def get_snp_in_gb(file,snps):
    snps_in_gb = []
    if my_ut.check_ext(file, "GenBank"):
        for rec in SeqIO.parse(file, "genbank"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        for snp in snps:
                            if snp in feature.location:
                                if not snp in snps_in_gb:
                                    snps_in_gb.append(snp)
                                    continue
        return snps_in_gb

    
print(get_snp_in_gb(file_gb,snp_ls))
len(get_snp_in_gb(file_gb,snp_ls))
snp_in_cds = get_snp_in_gb(file_gb,snp_ls)


# We've passed from 16k+ SNPs inside the CDS, to only 150 inside the coding region.

# In[ ]:


#to snp_utilities
#assign values to each snp based on the number of alleles [1:4] 

file_snp = "SNP_Locations.txt"
file_alleles = "SNP_Alleles.txt"
snp_in_cds = w


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




