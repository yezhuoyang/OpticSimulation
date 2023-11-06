import numpy as np
a=0.2
waveLength=633*10**(-9)
#The waist of z==0(Not the smallest waist)
w0=0.0016
zfo=a+0.415+0.155
fo=0.20 
fe=0.05
zfe=fo+fe+zfo



x1=a+0.058
x2=a+0.13
x3=a+0.3
x4=a+0.415+0.085
x5=a+0.415+0.085+0.3
w1=1.25
w2=1.5
w3=1.5
w4=1.75
w5=1.9
xlist=[x1,x2,x3,x4,x5]
wlist=[w1,w2,w3,w4,w5]


q0=1j*np.pi*w0**2/waveLength



def ABCD_free(x):
    return np.array([[1,x],[0,1]])

def ABCD_lens(f):
    return np.array([[1,0],[-1/f,1]])


def q_evolve(q,ABCD):
    A=ABCD[0][0]
    B=ABCD[0][1]
    C=ABCD[1][0]
    D=ABCD[1][1]
    return (A*q+B)/(C*q+D)



def calc_q(z):
    L=waveLength
    if z==0:
        return q0  
    if z<=zfo:
        ABCD=ABCD_free(z)
        return q_evolve(q0,ABCD)
    if z<=zfe:
        ABCD1=ABCD_free(zfo)
        ABCD2=ABCD_lens(fo)
        ABCD3=ABCD_free(z-zfo)
        ABCD=np.matmul(ABCD3,ABCD2)
        ABCD=np.matmul(ABCD,ABCD1)
        return q_evolve(q0,ABCD)
    ABCD1=ABCD_free(zfo)
    ABCD2=ABCD_lens(fo)
    ABCD3=ABCD_free(zfe-zfo)   
    ABCD4=ABCD_lens(fe)   
    ABCD5=ABCD_free(z-zfe) 
    ABCD=np.matmul(ABCD5,ABCD4)
    ABCD=np.matmul(ABCD,ABCD3) 
    ABCD=np.matmul(ABCD,ABCD2)  
    ABCD=np.matmul(ABCD,ABCD1)  
    return q_evolve(q0,ABCD)



def calc_wasit_0(z):
    ABCD=ABCD_free(z)
    q=q_evolve(q0,ABCD)
    return calc_waist_by_q(q)


def calc_waist_by_q(q):
    L=waveLength
    return np.sqrt(abs(-1*L/np.pi/np.imag(1/q)))



def calc_R_by_q(q):
    return 1/(np.real(1/q))


def calc_waist(z):
    return calc_waist_by_q(calc_q(z))

def calc_R(z):
    return calc_R_by_q((calc_q(z)))

import matplotlib.pyplot as plt

x = np.linspace(0.01, 2)
result=[]
for element in x:
    result.append(calc_waist(element)*1000)
plt.plot(x, result,label="The waist of laser beam")
plt.scatter(xlist,wlist)
plt.xlabel("z/m")
plt.ylabel("Waist/m")
plt.axvline(x=zfo,color='red',linestyle='--',label="The position of Objective Lens")
plt.axvline(x=zfe,color='green',linestyle='--',label="The position of eyepiece Lens")
plt.axhline(y=0.43,color='yellow',linestyle='--',label="The position of eyepiece Lens")


plt.title('The beamwaist of laser beam during propagation')
plt.legend()


