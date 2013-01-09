import csv,sys,random

afmFileIn = sys.argv[1]
afmFileOut = sys.argv[3]
pMissing = float(sys.argv[2])

assert afmFileIn != afmFileOut
assert 0 < pMissing < 1

afmReader = csv.reader(open(afmFileIn,'r'),delimiter='\t')
afmWriter = csv.writer(open(afmFileOut,'w'),delimiter='\t')

afmWriter.writerow(afmReader.next())

for inputLine in afmReader:
    
    afmWriter.writerow( [inputLine[0]] + [ "NA" if random.uniform(0,1) < pMissing else x for x in inputLine[1:] ] )
