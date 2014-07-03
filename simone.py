#temporanea simone

def getFasta(path):
    filefasta = open(path, 'r')
    lines = filefasta.readlines()
    filefasta.close()
    #trim first line
    fastatemp = lines[1:len(lines)]
    fasta = ''
    #remove \n
    for i in range(len(fastatemp)):
        fastatemp[i] = fastatemp[i].replace('\n', '')
        fasta = fasta + fastatemp[i]
    return fasta