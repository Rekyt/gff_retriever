#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that from a gene id of GFF file and specified length, gives you the number to extract the upstream sequence of this gene

### IMPORTS ###
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

def load_gff(gff_file):
	"""Returns a list of parsed gff."""
	with open(gff_file,"r") as f:
		print "Parsing file {}...".format(gff_file)
		rec = []
		for line in GFF.parse(f):
			rec.append(line)

	return rec

def load_fasta(fasta_file):
	"""Returns a list of parsed fasta."""
	with open(fasta_file,"r") as f:
		seqs = []
		for seq in SeqIO.parse(f,"fasta"):
			seqs.append(seq)
	return seqs

def fasta_len(fasta_file,seqid):
	"""Return length of sequence seqid in fasta_file."""

	length = 0

	# Find the length in fasta_file of seqid
	record = load_fasta(fasta_file)		
	for r in record:
		if record.id == seqid:
			length = len(record.seq) #length of the sequence

	if length == 0:
		print "Sequence {} was not found in {}".format(seqid,fasta_file)
	return length


def retrieve_seq(fasta_file,uplist):
	"""Take retrieve_up return and return the same list with sequences appended."""

	# Puts fasta file in memory and parses it
	records = load_fasta(fasta_file) # list of sequences


	# retrieves upstream sequence of each gene
	for u in uplist:
		seqid = u[-1]
		for rec in records:
			if rec.id == seqid:
				seq = rec.seq[u[0]:u[1]] # extract sequence, beware of indexes, as BioPython indexes from 0
				u.append(seq)

	return uplist

def write_fasta(file_name,upseqs):
	"""Write a fasta_file from upstream sequences list returned by retrieve_up."""

	records = []

	for u in upseqs:
		strand = u[2]
		gene = u[3] # gene name
		seqid = u[4] # name of scaffold from which the gene is extracted
		seq = u[-1] # Seq object

		
		if strand == 1:
			strand = "+"
		elif strand == -1:
			strand = "-"
		else:
			strand = "?"

		ident = seqid+"|"+gene+"|"+strand

		rec = SeqRecord(seq,id=ident,name=gene,description="")

		records.append(rec)

	SeqIO.write(records,file_name,"fasta")

def extract_cds(fasta_file,gff_file,T=None,table=None):
	"""Returns a list of Coding Sequences extracted from gff_file and fasta_file."""
	if T == None:
		T = False
	if table == None:
		table = 1
	# Load gff file in memory
	rec = load_gff(gff_file)

	fasta = load_fasta(fasta_file)

	# Create a list
	cds = retrieve_pos("CDS",rec) # list of CDS

	# Append sequences at the end of each entry in cds
	for c in cds:
		start = c[0]
		end = c[1]
		seq_id = c[-1]
		seq = seq_extract(seq_id,start,end,fasta)
		c.append(seq)

	
	# This loop assemble translated genes
	genes = []
	i = 0
	while i < len(cds)-1:
		start = cds[i][0]
		strand = cds[i][2]
		phase = cds[i][3]
		gene_id = cds[i][4]
		seq_id = cds[i][5]
		dna_seq = Seq("")
		dna = []
		dna.append(cds[i][-1])
		while cds[i+1][4] == gene_id and i+1 < len(cds)-1:
			i += 1
			dna.append(cds[i][-1])

		if strand == 1:
			for d in dna:
				dna_seq += d
		elif strand == -1:
			dna.reverse()
			for d in dna:
				dna_seq += d.reverse_complement()
		
		if T == True:

			prot_seq = dna_seq.translate(table=table)
			genes.append([start,phase,strand,gene_id,seq_id,prot_seq])

		else:
			genes.append([start,phase,strand,gene_id,seq_id,dna_seq])
		i += 1

	return genes

def seq_extract(seq_name,start,end,fasta_rec):
	"""Retrieves the DNA sequence in seq_name with positions start and end. fasta_rec is the results of load_fasta"""
	seq = Seq("")
	for f in fasta_rec:
		if f.id == seq_name:
			seq = f.seq[start:end]
	return seq

def retrieve_pos(seq_type,gff_rec):
	"""Return a list of positions sequences of given type using a parsed gff."""

	positions = []
	for r in gff_rec:
		seq_name = r.id # Scaffold name
		feat = r.features # list of all genes in scaffold
		for f in feat:
			gene_name = f.id
			subfeat = f.sub_features # list of all mRNA from this gene
			for s in subfeat:
				subsubfeat = s.sub_features # list of all exons, CDSs and introns
				for ss in subsubfeat:
					if ss.type == seq_type:
						start = ss.location.start.position # Beware it uses position in sequence using pythonic indexes
						end = ss.location.end.position
						strand = ss.location.strand
						phase = int("".join(ss.qualifiers["phase"]))

						posinfo = [start,end,strand,phase,gene_name,seq_name]
						positions.append(posinfo)
	return positions

def main():

	parser = argparse.ArgumentParser(description="Program to parse GFF files, extract all CDS and write a fasta file out of it.")
	parser.add_argument("gff",help="a .gff file, see specification on http://www.sequenceontology.org/gff3.shtml")
	parser.add_argument("fasta",help="a Fasta file")
	parser.add_argument("output",help="the name of the output file")
	parser.add_argument("-v","--verbose",help="verbose, add outputs in Terminal.",action="count",default=0)
	parser.add_argument("-t","--translate",help="return translated proteins if precised",action="store_true")
	choices = range(1,7) + range(9,17) + range(21,26)
	parser.add_argument("-tab","--table",help="NCBI codons table to use for translation. See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",type=int,choices=choices,default=1)

	args = parser.parse_args()

	if args.verbose >= 1:
		print "GFF file: {}\nFasta file: {}\nOutput file: {}".format(args.gff,args.fasta,args.output)
		print "Extracting CDS..."
		genes = extract_cds(args.fasta,args.gff,args.translate,args.table)
		print "Extracted {} CDSs".format(len(genes))
		print "Writing Fasta file..."
		write_fasta(args.output,genes)
		print "Fasta file {} written!".format(args.output)
	else:
		genes = extract_cds(args.fasta,args.gff,args.translate,args.table)
		write_fasta(args.output,genes)

if __name__ == "__main__":

	main()

