"""
% Caculate the decrease using the variance formula (slow+unstable)
if isNumerical
    DI = deltaImpurity_var_regr(x_left,x_right);

    %Calculate the decrease using the mean formulat (fast+stable)
    DI_test = deltaImpurity_mean_regr(x_left,x_right);

    %Make sure the two measures agree
    assert( abs(DI - DI_test ) < 1e-5, 'error: impurity functions disagree');
else

    DI = deltaImpurity_class(x_left,x_right);

end

function DI = deltaImpurity_mean_regr(x_left,x_right)

x = [x_left(:);x_right(:)];

mu = mean(x);
n = length(x);
muL = mean(x_left);
nL = length(x_left);
muR = mean(x_right);
nR = length(x_right);

DI = -mu^2 + nL/n*muL^2 + nR/n*muR^2;

function DI = deltaImpurity_var_regr(x_left,x_right)

x = [x_left(:);x_right(:)];
n = length(x);
nL = length(x_left);
nR = length(x_right);
DI = var(x,1) - nL/n*var(x_left,1) - nR/n*var(x_right,1);

function DI = deltaImpurity_class(x_left,x_right)
x = [x_left(:);x_right(:)];
n = length(x);
nL = length(x_left);
nR = length(x_right);
DI = giniIndex(x) - nL/n*giniIndex(x_left) - nR/n*giniIndex(x_right);

function GI = giniIndex(x)
GI = hist(x,unique(x))/length(x);
if ~isempty(GI)
    GI = 1 - sum(GI.^2);
else
    GI = 0;
end                                   
"""

import sys
import getopt
import numpy

def myHist(list):
	dic = {}
	for l in list:
		print l
		if (dic.get(l)):
			dic[l] = dic[l] + 1
		else:
			dic[l] = 1
	return dic

def giniIndex(x):
	print "Begin giniIndex"
	print x
	L = len(x)
	sorted_x = sorted(x)
	hist = myHist(sorted_x)
	"""
	numeric_sx = []
	for v in sorted_x:
		numeric_sx.append(float(v))
	print numeric_sx
	myset = set(x)
	y = numpy.cumsum(numeric_sx)
	B = sum(y) / (y[-1] * L)
	return 1 + 1./L - 2*B
	"""
	print hist.keys()
	GI = 0.0	
	for v in hist.values():
		GI = GI + pow(float(v)/float(L), 2)
	return 1 - GI;	

"""
function DI = deltaImpurity_mean_regr(x_left,x_right)

x = [x_left(:);x_right(:)];

mu = mean(x);
n = length(x);
muL = mean(x_left);
nL = length(x_left);
muR = mean(x_right);
nR = length(x_right);

DI = -mu^2 + nL/n*muL^2 + nR/n*muR^2;
"""

def diIndex(x,y):
	x = [float(v) for v in x]
	y = [float(v) for v in y]
	w = x + y
		
	L = len(w)
	xL = len(x)
	yL = len(y)
	mw = numpy.mean(w)
	mx = numpy.mean(x)
	my = numpy.mean(y)
	return (-1.0)*pow(mw,2) + float(xL)/L*(pow(mx,2)) + float(yL)/L*(pow(my,2))
	
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
    #for arg in args:
        #process(arg) # process() is defined elsewhere
    l = args[0].split(",")
    r = args[1].split(",")
    w = l + r
    left = giniIndex(l)
    right = giniIndex(r)
    print left
    print right
    wgi = giniIndex(w) - (float(len(l))/float(len(w))*left + float(len(r))/float(len(w))*right )
    print wgi

    print diIndex(l,r)	
    	
if __name__ == "__main__":
    main()

