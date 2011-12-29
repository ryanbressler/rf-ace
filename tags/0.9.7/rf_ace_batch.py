import sys,os

assert sys.argv[1] != sys.argv[4]
assert int(sys.argv[2]) <= int(sys.argv[3])

xlist = xrange(int(sys.argv[2]),int(sys.argv[3])+1)

for targetidx in xlist:
    os.system('bin/rf_ace --traindata '+sys.argv[1]+' --target '+str(targetidx)+' --associations '+sys.argv[4]+'_'+str(targetidx))

os.system('cat '+sys.argv[4]+'_* > '+sys.argv[4])
os.system('rm '+sys.argv[4]+'_*')
