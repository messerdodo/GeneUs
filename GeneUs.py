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



#imports fasta file (returns single string)
def getFasta(path):
    #opens file
    filefasta = open(path, 'r')
    lines = filefasta.readlines()
    filefasta.close()
    #trims first line
    fasta = lines[1:len(lines)]
    #removes '\n'
    for i in range(len(fasta)):
        fasta[i] = fasta[i].replace('\n', '')
    #characters per line
    cpl = len(fasta[1])
    return fasta, cpl

#returns a substring of the fasta
#fasta is used as a matrix (rows, character position)
#unique characters-per-line value is supposed
def getFastaString(begin, end, fasta, cpl):
    fastastring = ''
    #character position in the row is determined by its offset
    offset_begin = begin % cpl
    offset_end = end % cpl
    #calculates the row that containes the character
    index_begin = (begin - (offset_begin))/cpl
    index_end = (end - (offset_end))/cpl
    #temp vector of fasta
    temp = fasta[index_begin:index_end]
    #trims the offsets of the first and the last row
    temp[0] = temp[0][offset_begin:]
    temp[len(temp)-1] = temp[len(temp)-1][offset_end:]
    #updates fastastring
    for i in range(0,len(temp)):
        fastastring = fastastring + temp[i]
    return fastastring

def demo():
    fasta, cpl = getFasta('ENm006.fa')
    exons, cds = GTFParsing('GAB3_annot.gtf');
    print 'exons:\n', exons;
    print 'cds:\n', cds;
    print 'fasta:\n', fasta, '\ncpl:\n', cpl
    fastastring = getFastaString(1650, 2040, fasta, cpl)
    print 'fastastring:\n', fastastring

if __name__ == '__main__':
    demo()
