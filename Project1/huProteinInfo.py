#!/usr/bin/python
#
# Author: Travis Weaver
# Date: November 11, 2016
# Last Updated: November 11, 2016
#
# Program to pull additional information about protein coding genes and export the
# information to a fasta file. Additional information includes approved names,
# gene locations, and synonyms.

from Bio import SeqIO
from Bio import Entrez
import sys, string, textwrap

proteinNameList = []
proteinDict = {}
is_verbose = 0
iterator = 2

for seq_record in SeqIO.parse("human.fa", "fasta"):
	if is_verbose == 1:
		print(seq_record.id)
		print(repr(seq_record.seq))
		print(len(seq_record))

	proteinNameList.append(seq_record.id)
	proteinDict[seq_record.id] = seq_record.seq

proteinFile = open("protein-coding_gene.txt", 'r')
proteinFileContents = proteinFile.readlines()
proteinFile.close()

if is_verbose == 1:
	for proteinName, proteinSeq in proteinDict.iteritems():
		print proteinName, proteinSeq

fastaFileWriter = open("humanProteinInfo.fa", 'a')

for i in proteinFileContents:
	proteinFileLines = string.split(i)

	for j in proteinNameList:
		if j == proteinFileLines[1]:
			fastaFileWriter.write(">")
			fastaFileWriter.write(j)
			fastaFileWriter.write(" | ")

			while proteinFileLines[iterator] != "Approved":
				fastaFileWriter.write(proteinFileLines[iterator])
				fastaFileWriter.write(" ")
				iterator += 1

			iterator += 1

			while iterator < (len(proteinFileLines) - 2):
				fastaFileWriter.write(proteinFileLines[iterator])
				fastaFileWriter.write(" ")
				iterator += 1

			fastaFileWriter.write("| ")
			fastaFileWriter.write(proteinFileLines[iterator])
			fastaFileWriter.write("\n")

			wrappedSeq = textwrap.fill(str(proteinDict[j]), 60)

			fastaFileWriter.write(wrappedSeq)
			fastaFileWriter.write("\n")

		iterator = 2

fastaFileWriter.close()

# End of code