#Paroni Pellegrini Previtali

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
    return fasta

def demo():
    fasta = getFasta('ENm006.fa')
    print fasta

if __name__ == '__main__':
    demo()