import numpy as np
from scipy import special
import matplotlib.pyplot as plt


wx=0.001
P=1
def Power(x):
    return 0.5*P*(1-special.erf(np.sqrt(2)*x/wx))


import scipy
def c(x,wx,P):
    return 0.5*P*(1-scipy.special.erf(np.sqrt(2)*x/wx))

def calcx(a,b):
    return 0.5*a+0.02*b


xlist2=[(4,19),(4,11),(4,7.5),(4,5),(4,2.5),(3,24.5),(3,21),(3,9.5)]
xlist2=[calcx(a[0],a[1]) for a in xlist2]
xlist2=[x-min(xlist2) for x in xlist2]

Plist2=[0.03,0.520,1,1.5,2,2.53,3.01,3.41]

w0_2,P_2=scipy.optimize.curve_fit(Power,xlist2,Plist2)[0]
