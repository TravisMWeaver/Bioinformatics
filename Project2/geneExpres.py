#!/usr/bin/python
#
# Author: Travis Weaver
# Date: December 9, 2016
# Last Updated: December 9, 2016
#

from scipy import stats
import sys, string, math
import numpy as np

geneFile = open("Influenza_inf_expr.txt", 'r')
infGeneData = geneFile.readlines()
geneFile.close()

geneFile = open("Influenza_con_expr.txt", 'r')
conGeneData = geneFile.readlines()
geneFile.close()

iterator = 0

geneFileWriter = open("SigExpGenes.txt", 'a')

for i, j in zip(infGeneData, conGeneData):
	infGeneSplit = string.split(i)
	conGeneSplit = string.split(j)

	geneName = infGeneSplit[0]
	print(geneName)

	infVals = infGeneSplit[1:]
	conVals = conGeneSplit[1:]

	if iterator != 0:
		infVals = map(float, infVals)
		conVals = map(float, conVals)

		t_test = stats.ttest_ind(infVals, conVals)
		print "The t-statistic is %.3f and the p-value is %.3f." % t_test

		if t_test[1] < 0.01:
			geneFileWriter.write("Gene: ")
			geneFileWriter.write(geneName)
			geneFileWriter.write("\n")
			geneFileWriter.write("P-Value: ")
			geneFileWriter.write(str(t_test[1]))
			geneFileWriter.write("\n")

	iterator += 1

#
# End of file