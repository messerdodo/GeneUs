#Paroni Pellegrini Previtali

import csv as tsv
import string;

###############################################################################
## This method reads a GTF file and returns two structured object, one for   ##
## the exons and one for the CDS.                                            ##
## The exons object is a list where each item is an array with the following ##
## structure: begin_position, end_position, strand (+/-), transcript_id.     ##
## The CDS object is structured as follow:                                   ##
## begin_position, end_position, strand (+/-), transcript_id.                ##
###############################################################################
def GTFParsing(gtfFile) :

	#Exons containter
	exons = [];#4
	cds = []; #7

	#Opens the GTF file
	with open(gtfFile) as gtf:
		for line in tsv.reader(gtf, delimiter = '\t'):
			if line[2] == 'exon':
				#New exon data
				exon = [''] * 4;
				#Exon begin
				exon[0] = int(line[3]);
				#Exon end
				exon[1] = int(line[4]);
				#Strand (+/-)
				exon[2] = line[6];
				#Transcript Id
				transcriptId = line[8].split(';')[0];
				exon[3] = transcriptId[transcriptId.find('"') + 1 : len(transcriptId) - 1];
				#Adds the exons to the list
				exons = exons + [exon];
			elif line[2] == 'CDS':
				#New CDS data
				singleCDS = [''] * 4;
				#CDS begin
				singleCDS[0] = int(line[3]);
				#CDS end
				singleCDS[1] =  int(line[4]);
				#Strand
				singleCDS[2] = line[6];
				annotations = line[8].split(';');
				#Transcript Id
				transcriptId = annotations[0];
				singleCDS[3] = transcriptId[transcriptId.find('"') + 1 : len(transcriptId) - 1];
				#Adds the CDS to the collection
				cds = cds + [singleCDS];

	return exons, cds;



#import fasta file (returns single string)
def getFasta(path):
    filefasta = open(path, 'r')
    lines = filefasta.readlines()
    filefasta.close()
    #trim first line
    fasta = lines[1:len(lines)]
    #remove \n
    for i in range(len(fasta)):
        fasta[i] = fasta[i].replace('\n', '')
    #characters per line
    cpl = len(fasta[1])
    return fasta, cpl

def demo():
    fasta, cpl = getFasta('ENm006.fa')
    exons, cds = GTFParsing('GAB3_annot.gtf');
    print exons;
    print cds;
    print fasta, cpl

if __name__ == '__main__':
    demo()
