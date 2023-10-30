waveLength=633*10**(-9)
#The waist of z==0(Not the smallest waist)
wini=0.00305
zfo=0.06
fo=0.20 
fe=0.05
zfe=fo+fe+zfo


def calc_q(z):
    L=waveLength
    if z==0:
        return -1j*L/np.pi/wini**2    
    if z<=zfo:
        return z+calc_q(0)
    if z<=zfe:
        qfo_before=calc_q(zfo)
        qfo_after=qfo_before-1/fo
        return 1/(1/(qfo_after)+z-zfo)
    qfe_before=calc_q(zfe)
    qfe_after=qfe_before-1/fe
    return 1/(1/(qfe_after)+z-zfe)

def calc_q_transe(ABCD,q1inv):
    return (ABCD[1][0]+ABCD[1][1]*q1inv)/(ABCD[0][0]+ABCD[0][1]*q1inv)

def calc_waist_by_q(qinv):
    L=waveLength
    return np.sqrt(abs(-1*L/np.pi/np.imag(qinv)))

def calc_R_by_q(qinv):
    return 1/(np.real(qinv))

def calc_ABCD_lens(f):
    return np.array([[1,0],[-1/f,1]])

def calc_ABCD_freeprop(L):
    return np.array([[1,L],[0,1]])

def calc_q_ABCD(z):
    L=waveLength
    if z==0:
        return -1j*L/np.pi/wini**2    
    if z<=zfo:
        ABCD=calc_ABCD_freeprop(z)
        return calc_q_transe(ABCD,calc_q_ABCD(0))
    if z<=zfe:
        ABCD_fo=calc_ABCD_lens(fo)
        ABCD_freeprop=calc_ABCD_freeprop(z-zfo)
        return calc_q_trans(ABCD_freeprop,calc_q_trans(ABCD_fo,calc_q_ABCD(zfo)))
    ABCD_fe=calc_ABCD_lens(fe)
    ABCD_freeprop=calc_ABCD_freeprop(z-zfe)    
    return calc_q_trans(ABCD_freeprop,calc_q_trans(ABCD_fe,calc_q_ABCD(zfe)))

def calc_waist(z):
    return calc_waist_by_q(calc_q(z))

def calc_R(z):
    return calc_R_by_q((calc_q(z)))

import matplotlib.pyplot as plt

x = np.linspace(0.001, 0.3)
result=[]
for element in x:
    result.append(calc_R(element))
plt.plot(x, result,label="The curvature of laser beam")
plt.xlabel("z/m")
plt.ylabel("Curvature/m")
plt.axvline(x=0.1,color='red',linestyle='--',label="The position of Objective Lens")
plt.axvline(x=0.175,color='green',linestyle='--',label="The position of eyepiece Lens")
plt.title('The curvature of laser beam during propagation')
plt.legend()


x = np.linspace(0.0,1 ,1000)
result=[]
for element in x:
    result.append(calc_waist(element)*1000)
plt.plot(x, result,label="The waist of laser beam")
plt.xlabel("z/m")
plt.ylabel("Waist/mm")
plt.axvline(x=zfo,color='red',linestyle='--',label="The position of Objective Lens")
plt.axvline(x=zfe,color='green',linestyle='--',label="The position of eyepiece Lens")
plt.axvline(x=zfo+0.52+0.126,color='yellow',linestyle='--',label="The position of waist measurement1")
plt.scatter(zfo+0.52+0.126,0.933,marker='*',s=80,c='yellow')
plt.axvline(x=zfo+0.52+0.156,color='purple',linestyle='--',label="The position of waist measurement2")
plt.scatter(zfo+0.52+0.156,1.0946,marker='*',s=80,c='purple')
plt.axvline(x=zfo+0.52+0.18,color='black',linestyle='--',label="The position of waist measurement2")
plt.scatter(zfo+0.52+0.18,1.1314,marker='*',s=80,c='black')
plt.title('The waist of laser beam during propagation')
plt.legend()