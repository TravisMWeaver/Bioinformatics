#!/usr/bin/python
#
# Author: Travis Weaver
# Date: November 11, 2016
# Last Updated: November 11, 2016
#
# Program to compare given protein sequences found in humans to protein sequences 
# found in common house mice. See BLASTinfo.txt.

from Bio import SeqIO
from Bio import Entrez
Entrez.email = "t.weaver3@umiami.edu"
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys, string, textwrap

proteinNameList = []
proteinDict = {}
is_verbose = 0

resultsFileWriter = open("huMouseProteinCmp.txt", 'a')

for seq_record in SeqIO.parse("human.fa", "fasta"):
	resultsFileWriter.write("-----Initial Input-----\n")
	resultsFileWriter.write("Human Protein ID: ")
	resultsFileWriter.write(seq_record.id)
	resultsFileWriter.write("\nProtein Sequence: ")
	resultsFileWriter.write(repr(seq_record.seq))
	resultsFileWriter.write("\nLength: ")
	resultsFileWriter.write(str(len(seq_record)))

	if is_verbose == 1:
		print("-----Initial Input-----")
		print("Human Protein ID: ", seq_record.id)
		print("Protein Sequence: ", repr(seq_record.seq))
		print("Length: ", len(seq_record))

	proteinNameList.append(seq_record.id)
	proteinDict[seq_record.id] = seq_record.seq
	result_handle = NCBIWWW.qblast("blastp", "nr", seq_record.seq, entrez_query="mouse[orgn]", matrix_name="BLOSUM62", hitlist_size=1)

	blast_record = NCBIXML.read(result_handle)

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < 0.001:
				resultsFileWriter.write("\n\n")
				resultsFileWriter.write("*****Alignment Starts Here*****\n")
				resultsFileWriter.write("\nSequence: ")
				resultsFileWriter.write(alignment.title)
				resultsFileWriter.write("\nLength: ")
				resultsFileWriter.write(str(alignment.length))
				resultsFileWriter.write("\nE Value: ")
				resultsFileWriter.write(str(hsp.expect))
				resultsFileWriter.write("\nRaw Score: ")
				resultsFileWriter.write(str(hsp.score))
				resultsFileWriter.write("\nBitscore: ")
				resultsFileWriter.write(str(hsp.bits))
				resultsFileWriter.write("\nQuery Sequence: ")
				resultsFileWriter.write(hsp.query[0:75])
				resultsFileWriter.write("...")
				resultsFileWriter.write("\nMatch Sequence: ")
				resultsFileWriter.write(hsp.match[0:75])
				resultsFileWriter.write("...")
				resultsFileWriter.write("\nSubject Sequence: ")
				resultsFileWriter.write(hsp.sbjct[0:75])
				resultsFileWriter.write("...")
				resultsFileWriter.write("\n\n")

				if is_verbose == 1:
					print("\n")
					print("*****Alignment Starts Here*****")
					print("Sequence: ", alignment.title)
					print("Length: ", alignment.length)
					print("E Value: ", hsp.expect)
					print("Raw Score: ", hsp.score)
					print("Bitscore: ", hsp.bits)
					print("Query Sequence: ", hsp.query[0:75] + "...")
					print("Match Sequence: ", hsp.match[0:75] + "...")
					print("Subject Sequence: ", hsp.sbjct[0:75] + "...")
					print("\n")

resultsFileWriter.close()

#fasta_string = open("human.fa").read()
#result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string, hitlist_size=1)
#result_handle = NCBIWWW.qblast("blastn", "nt", "8332116", hitlist_size=1)

#save_file = open("my_blast.xml", "w")
#save_file.write(result_handle.read())
#save_file.close()
#result_handle.close()

#
#
# End of code