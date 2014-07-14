#Paroni Andrea
#Pellegrini Simone Maria
#Previtali Giorgia

import csv as tsv
import string;
import os
import textwrap
import sys;
import getopt;

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
	cds = []; #4

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
## This method returns an array of lines of a FASTA file                     ##
## Variable @cpl is the number of characters per line of the fasta file. It  ##
##  is assumed that this value is unique for all the FASTA file              ##
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
        fasta[i] = fasta[i].replace('\r', '')
    #characters per line
    cpl = len(fasta[1])
    return fasta, cpl

###############################################################################
## This method returns a substring of the given FASTA file (variables @fasta ##
##  and @cpl). Length of the substring is given by variables @begin and @end ##
## Variable @fasta is used as a matrix (lines, character offset)             ##
## Variable @cpl (characters per line) is the number of nucleic bases on a   ##
##  single line of FASTA file                                                ##
## The substring is stored in the output variable @fastastring.              ##
###############################################################################
def getFastaString(begin, end, fasta, cpl):
    fastastring = ''
    begin = begin - 1; #zero based
    end = end - 1; #zero based

    #character position in the row is determined by its offset
    offset_begin = begin % cpl
    offset_end = end % cpl
    #calculates the row that containes the character
    index_begin = begin/cpl
    index_end = end/cpl
    #Begin and end in the same row.
    if index_begin == index_end :
    	fastastring = fasta[index_begin][offset_begin : offset_end + 1];
    else:
		#temp vector of fasta
    	temp = [];
    	if index_begin + 1 <= index_end - 1:
    		temp = fasta[index_begin + 1 : index_end]
    	#trims the offsets of the first and the last row
    	fastastring = fasta[index_begin][offset_begin:]; 
    	#updates fastastring
    	for i in range(len(temp)):
        	fastastring = fastastring + temp[i];
    	fastastring = fastastring + fasta[index_end][:(offset_end + 1)];
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

###############################################################################
## This method sorts the given annotations list by the begin value.          ##
## The input variable @couples is a list of couples (begin, end).            ##
## Variable @increasing is set to True (default value) if increasing sort    ##
##  is needed; it can be set to False if decreasing sort is needed.          ##
## Sorted annotations are stored in the output variable @sorted_couples.     ##
###############################################################################
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
	
###############################################################################
## This method returns the introns list, given the exons annotiations.       ##
## The introns and the exons are taken from the given fasta file (variables  ##
##  @fasta and @cpl).                                                        ##
## The introns list is stored in the output variable @introns                ##
###############################################################################
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
## annotations.                                                              ##
###############################################################################
def getsAnnotationsGrouppedByTranscriptId(annotations):
	transcriptIds = [];
	grouppedAnnotations = [];
	strands = [];
	#Analizes all the annotations and classifies them by the transcript id
	for annotation in annotations:
		#The transcript is already seen
		if annotation[3] in transcriptIds:
			group = transcriptIds.index(annotation[3]);
			grouppedAnnotations[group] = grouppedAnnotations[group] + [(annotation[0], annotation[1])];
		else:
			#New transcript id.
			grouppedAnnotations = grouppedAnnotations + [[(annotation[0], annotation[1])]];
			transcriptIds = transcriptIds + [annotation[3]];
			strands = strands + [annotation[2]];

	return transcriptIds, grouppedAnnotations, strands;

###############################################################################
## This method returns all the transcripts relative to the given annotations ##
##  (variable @annotations). The transcripts identifiers are stored in the   ##
##  output variable @transcriptsIds. The the transcripts' nucleic bases are  ##
##  stored in the output variable @transcripts.                              ##
## Transcripts are taken from the given genomic (variables @fasta and @cpl). ##
###############################################################################
def getTranscripts(annotations, fasta, cpl):
	transcriptsIds, grouppedannotations, strands = getsAnnotationsGrouppedByTranscriptId(annotations);
	i = 0;
	transcripts = [];
	#Analizes each transcipt
	for annotationsInTranscript in grouppedannotations:
		transcript = '';
		#Crick strand: A reverse and complement task is required
		if strands[i] == '-':
			#Reverse way.
			annotationsInTranscript = sortByBegin(annotationsInTranscript, increasing = False);
			#Gets the sequence for each exon and do the reverse and complement task
			for annotation in annotationsInTranscript:
				transcript = transcript + reverseAndComplement(getFastaString(annotation[0], annotation[1], fasta, cpl));	
		else: #Watson strand
			#Standard way.
			annotationsInTranscript = sortByBegin(annotationsInTranscript);
			#Gets the sequence for each exon and produces the transcript sequence.
			for annotation in annotationsInTranscript:
				transcript = transcript + getFastaString(annotation[0], annotation[1], fasta, cpl);
		transcripts = transcripts + [transcript];
		i = i + 1;
	return transcriptsIds, transcripts;

###############################################################################
## This method returns all the Coding Sequences connected to the given       ##
## annotations (variable @annotations).                                      ##
## CDSs are taken from the given genomic (variables @fasta and @cpl).        ##
###############################################################################
def getCDS(annotations, fasta, cpl):
	transcriptsIds, grouppedAnnotations , strands = getsAnnotationsGrouppedByTranscriptId(annotations);
	cdsSequences = [];
	i = 0;
	for annotationsInTranscript in grouppedAnnotations :
		singleCDS = '';
		#Crick strand: A reverse and complement task is required
		if strands[i] == '-':
			#Reverse way.
			annotationsInTranscript = sortByBegin(annotationsInTranscript, increasing = False);
			#Gets the sequence for each exon and do the reverse and complement task
			for annotation in annotationsInTranscript:
				singleCDS = singleCDS + reverseAndComplement(getFastaString(annotation[0], annotation[1], fasta, cpl));
		else: #Watson strand
			#Standard way.
			annotationsInTranscript = sortByBegin(annotationsInTranscript);
			#Gets the sequence for each exon and produces the transcript sequence.
			for annotation in annotationsInTranscript:
				singleCDS = singleCDS + getFastaString(annotation[0], annotation[1], fasta, cpl);
		cdsSequences = cdsSequences + [singleCDS];
		i = i + 1;
	return cdsSequences, transcriptsIds;

###############################################################################
# This method writes as many FASTA file as variable @ids                      #
# Each file is filled with @contents[i] and stored in a directory named after #
#  @directory. The width of each line is limited by @cpl.                     #
###############################################################################
def fastaExport(ids, contents, cpl, directory = '.', fileName = 'file'):
    #checks weather the directory already exists or not
    if not os.path.exists(directory):
        #no directory, so creates it
        os.makedirs(directory)
    path = directory + "/" + fileName + '.fa';
  	#Creates a new FASTA file.
    with open(path, 'w') as fasta:
 		for i in range(len(ids)):
			#for every id, saves the trasncript in a fasta file
			current_id = ids[i]
			heading = '>' + current_id
			#writes head line of fasta
			fasta.write(heading + '\n')
			fasta.flush();
			#wraps the text width
			fasta.write(textwrap.fill(contents[i], width = cpl));
			fasta.flush();
			fasta.write('\n');
			fasta.flush();
    return

def main(argv):
	fastaFile = '';
	gtfFile = '';
	outputFolder = '';
	try:
		opts, args = getopt.getopt(argv, 'f:g:o:h', ['fasta=', 'gtf=', 'ofolder=', 'help'])
	except getopt.GetoptError:
		print 'GeneUs.py -f <fastafile> -g <gtffile> -o <outputfolder>';
		sys.exit(2);

	
	for opt, arg in opts:
		if opt == '-h':
			print 'GeneUs.py -f <fastafile> -g <gtffile> -o <outputfolder>';
			sys.exit()
		elif opt in ("-f", "--fasta"):
			fastaFile = arg;
		elif opt in ("-g", "--gtf"):
			gtfFile = arg;
		elif opt in ("-o", '--ofolder'):
			outputFolder = arg;

	print "GeneUs initializing...";

	#Reads the input files
	fasta, cpl = getFasta(fastaFile);
	exons, cds = GTFParsing(gtfFile);
	#Retrieves the introns and stores them
	introns = getIntrons(exons, fasta, cpl);
	intronsHeaders = [];
	intronsSeqs = [];
	for i in range(len(introns)):
		intronsHeaders = intronsHeaders + ['Intron {num}|Begin: {begin}|End: {end}'.format(num = i + 1, begin = introns[i][0], end = introns[i][1])]; 
		intronsSeqs = intronsSeqs + [introns[i][2]];

	print "Saving introns at " + outputFolder + "/introns.fa ...";
	fastaExport(intronsHeaders, intronsSeqs, cpl, outputFolder, fileName = "introns");
	print "Introns saved at " + outputFolder + "/introns.fa";

	#Retrieves the transcripts and saves them
	transcriptsIds, transcripts = getTranscripts(exons, fasta, cpl);
	print "Saving transcripts at " + outputFolder + "/transcripts.fa ...";
	fastaExport(transcriptsIds, transcripts, cpl, outputFolder, fileName = "transcripts");
	print "Transcripts saved at " + outputFolder + "/transcripts.fa";

	#Retrieves the cds and saves them
	cdsSeqs, transcriptsIds = getCDS(cds, fasta, cpl);
	cdsHeaders = [];
	#Prepares the FASTA header
	for i in range(len(cdsSeqs)):
		cdsHeaders = cdsHeaders + ['CDS %d' %(i + 1) + "|Transcript id: " + transcriptsIds[i]];
	print "Saving coding sequences at " + outputFolder + "/CDS.fa ...";
	fastaExport(cdsHeaders, cdsSeqs, cpl, outputFolder, fileName = "CDS");
	print "Coding Sequences saved at " + outputFolder + "/CDS.fa";

	print "Task completed!"


if __name__ == '__main__':
	main(sys.argv[1:]);