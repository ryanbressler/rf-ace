import sys,os

assert sys.argv[1] != sys.argv[4]
assert int(sys.argv[2]) <= int(sys.argv[3])

xlist = xrange(int(sys.argv[2]),int(sys.argv[3])+1)

for targetidx in xlist:
    os.system('bin/rf_ace --input '+sys.argv[1]+' --targetidx '+str(targetidx)+' -n 500 -t 1 --output '+sys.argv[4]+'_'+str(targetidx)+'.out')

os.system('cat '+sys.argv[4]+'_*.out > '+sys.argv[4]+'.out')
os.system('rm '+sys.argv[4]+'_*.out')
