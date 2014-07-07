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


###############################################################################
# This method returns an array of lines of a FASTA file                       #
# Variable @cpl is the number of characters per line of the fasta file. It is #
#  assumed that this value is unique for all the FASTA file                   #
###############################################################################
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

##############################################################################
# This method returns a substring of FASTA                                   #
# Variable @fasta is used as a matrix (lines, character offset)              #
# Variable @cpl (characters per line) is the number of nucleic bases on a    #
#  single line of FASTA file                                                 #
##############################################################################
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


###############################################################################
## This method returns the complements of the passed nucleic acid.           ##
###############################################################################
def complements(base):
	if base.upper() == 'A':
		return 'T';
	elif base.upper() == 'T':
		return 'A';
	elif base.upper() == 'C':
		return 'G';
	else:
		return 'C';

###############################################################################
## This method returns the passed sequence reverted and complemented.        ##
###############################################################################
def reverseAndComplement(sequence):
	newSequence = '';
	lenght = len(sequence);
	for i in range(lenght):
		#Reverts and complements
		newSequence = newSequence + complements(sequence[lenght - 1 - i]);
	return newSequence;

##############################################################################
# This method sorts the given annotation list by the begin value.            #
# Variable @increasing is set to True (default value) if increasing sort     #
#  is needed; it can be set to False if decreasing sort is needed.           #
##############################################################################
def sortByBegin(couples, increasing = True ):
	#couples[0] contains begin position
	#couples[1] contains end position
	if increasing:
	    #sorts couples in increasing order (default)
            sorted_couples = sorted(couples, key=lambda couples: couples[0])
	else:
	    #sorts couples in decreasing order
	    sorted_couples = sorted(couples, key=lambda couples: couples[0], reverse = True)
	return sorted_couples
	
##############################################################################
# This method returns the introns list, given the exons annotiations.        #
##############################################################################
def getIntrons(exons, fasta, cpl):
	introns = [];
	#Sorts the exons by the begin attribute.
	exons = sortByBegin(exons);
	#Checks if the gene is encoded in the watson or crick strand
	crickStrand = exons[0][2] == '-';
	#Obtains the not overlapped exons
	preprocessedExons = [(exons[0][0], exons[0][1])];
	j = 0;
	for i in range(1, len(exons)):
		#The begin of the i-th exon is before the end of the previous exon (overlapping)
		if exons[i][0] <= preprocessedExons[j][1]:
			preprocessedExons[j] = (preprocessedExons[j][0], max(preprocessedExons[j][1], exons[i][0]));
		else:
			#No overlapping
			preprocessedExons = preprocessedExons + [(exons[i][0], exons[i][1])];
			j = j + 1;
	
	#Obtains all the introns (as gaps between the exons)
	for i in range(len(preprocessedExons) - 1):
		#Gets the intron begin and end.
		intronBegin = preprocessedExons[i][1] + 1;
		intronEnd = preprocessedExons[i + 1][0] - 1;
		#Gets the corresponding intron pattern
		intronPattern = getFastaString(intronBegin, intronEnd, fasta, cpl);
		#If the gene is in the crick strand executes the reverse and complement task.
		if crickStrand:
			intronPattern = reverseAndComplement(intronPattern);
		introns = introns + [(intronBegin, intronEnd, intronPattern)];

	return introns;

###############################################################################
## This method returns a structure containing the trascripts id and their    ##
## exons.                                                                    ##
###############################################################################
def getsExonsGrouppedByTranscriptId(exons):
	transcriptIds = [];
	grouppedExons = [];
	#Analizes all the exons and classifies them by the transcript id
	for exon in exons:
		#The transcript is already seen
		if exon[3] in transcriptIds:
			group = transcriptIds.index(exon[3]);
			grouppedExons[group] = grouppedExons[group] + [(exon[0], exon[1])];
		else:
			#New transcript id.
			grouppedExons = grouppedExons + [[(exon[0], exon[1])]];
			transcriptIds = transcriptIds + [exon[3]];

	return transcriptIds, grouppedExons;

###############################################################################
## This method returns a structures containing each transcript in a FASTA    ##
## format.                                                                   ##
############################################################################### 
def getTranscripts(exons, fasta, cpl):
	transcriptsIds, grouppedExons, strands = getsExonsGrouppedByTranscriptId(exons);
	i = 0;
	transcripts = [];
	#Analizes each transcipt
	for exonsInTranscript in grouppedExons:
		transcript = '';
		#Crick strand: A reverse and complement task is required
		if strands[i] == '-':
			#Reverse way.
			exonsInTranscript = sortByBegin(exonsInTranscript, increasing = False);
			#Gets the sequence for each exon and do the reverse and complement task
			for exon in exonsInTranscript: 
				transcript = transcript + reverseAndComplement(getFastaString(exon[0], exon[1], fasta, cpl));	
		else: #Watson strand
			#Standard way.
			exonsInTranscript = sortByBegin(exonsInTranscript);
			#Gets the sequence for each exon and produces the transcript sequence.
			for exon in exonsInTranscript: 
				transcript = transcript + getFastaString(exon[0], exon[1], fasta, cpl);
		transcripts = transcripts + [transcript];
		i = i + 1;
	return transcriptsIds, transcripts;


def demo():
    fasta, cpl = getFasta('ENm006.fa')
    exons, cds = GTFParsing('GAB3_annot.gtf');
    print 'exons:\n', exons;
    print 'cds:\n', cds;
    #tests
    #fastastring = getFastaString(1650, 2040, fasta, cpl)
    #print fastastring
    #rec_fastastring = reverseAndComplement(fastastring)
    #print rec_fastastring
    #print sortByBegin(exons)
    #print sortByBegin(exons, False)
    introns = getIntrons(exons, fasta, cpl);
    print 'Introns:'
    for intron in introns:
    	print intron;
    print "Transcripts:";
    transcriptsIds, transcripts = getTranscripts(exons, fasta, cpl);
    for i in range(len(transcriptsIds)):
    	print 'Transcript: ' + transcriptsIds[i];
    	print transcripts[i];

if __name__ == '__main__':
    demo()
