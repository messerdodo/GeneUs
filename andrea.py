import csv as tsv
import string;

###############################################################################
## This method reads a GTF file and returns a structured object.              #
##
###############################################################################
def GTFParsing(gtfFile) :

	#Exons containter
	exons = [];#5
	cds = []; #7

	#Opens the GTF file
	with open(gtfFile) as gtf:
		for line in csv.reader(gtf, delimiter = '\t'):
			if line[2] == 'exon':
				#New exon data
				exon = [''] * 5;
				exon[0] = int(line[3]);
				exon[1] = int(line[4]);
				exon[2] = line[6];
				transcriptId = line[8].split(';')[0];
				exon[3] = transcriptId[transcriptId.find('"') + 1 : len(transcriptId) - 1];
				exons = exons + [exon];
			elif line[2] == 'CDS':
				print 'Ciao';

	return exons;
				

