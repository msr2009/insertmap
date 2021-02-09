"""
convert_to_fasta.py

stupid wrapper for biopython SeqIO.convert
to convert genbank/ape/embl format files to fasta

Matt Rich, 1/2021
"""
import sys
from Bio import SeqIO
from re import sub

def remove_characters(s, characters="\(\)\{\}\.\*\:"):
	return sub("["+characters+"]", "", s)


for f in sys.argv[1:]:
		#grab format from filename
		seqformat = f.split(".")[-1]
		path = "/".join(f.split("/")[0:-1])


		#create output file name, same as chromosome name
		#filename is first 16 characters of chromosome, which is how idxstats
		#outputs the name (this is used late in the pipeline). Unfortunately
		#this makes these filenames a bit less readable... 
		out_stmt = "Converting {} from ".format(f)

		record = None
		if seqformat in ["fa", "fsa", "fasta"]: #== "fa" or seqformat == "fsa" or seqformat == "fasta":
				record = SeqIO.read(f, "fasta")
				out_stmt += "fasta to "
		elif seqformat in ["gb", "gbk", "ape"]: # == "gb" or seqformat == "gbk" or seqformat == "ape":
				record = SeqIO.read(f, "genbank")
				out_stmt += "genbank to "
		elif seqformat == "embl":
				record = SeqIO.read(f, "embl")
				out_stmt += "embl to "
		else:
				print("{} has unknown format. Requires .fa, .fsa, .fasta, .gb, .gbk, .ape, or .embl as file type.".format(f), 
						file=sys.stderr)
				continue

		#write to new file, 
		#change chromosome name to remove any special characters
		new_chr = remove_characters(record.id)
		record.id = new_chr

		#output new fasta file
		output_file = path + "/" + new_chr[0:16]+".fa"
		out_stmt = out_stmt + "fasta (" + output_file + ")"
		print(out_stmt, file=sys.stderr) 
		
		SeqIO.write(record, output_file, "fasta")
		

