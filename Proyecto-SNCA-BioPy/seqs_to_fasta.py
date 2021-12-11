from Bio import SeqIO
import my_utilities as my_ut

file1="SNCA_ReferenceGRCh38_HomoSapiens.gb"
file2="SNCA_RefSeq_HomoSapiens.gb"

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


get_seqs_cds(file1,"dna")