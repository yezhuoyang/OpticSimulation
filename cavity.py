import numpy as np

def ABCD_free(L):
    return np.array([[1,L],[0,1]])

def ACBD_con_mir(R):
    return np.array([[1,0],[-2/R,1]])

def ABCD_cavity_circle(a,b,R):
    ABCD_free_1=ABCD_free(2*b+np.sqrt(a**2+b**2))
    ABCD_con=ACBD_con_mir(R) 
    ABCD=np.matmul(ABCD_con,ABCD_free_1)
    ABCD_free_2=ABCD_free(np.sqrt(a**2+b**2))
    return  np.matmul(ABCD_free_2,ABCD)

def ABCD_cavity(a,b,R,z):
    ABCD_once=ABCD_cavity_circle(a,b,R)
    circle=2*b+2*np.sqrt(a**2+b**2)
    period=int(z/circle)
    matrix=ABCD_once
    for i in range(0,period):
        matrix=np.matmul(ABCD_once,matrix)
    ABCD=matrix
    remain=z-period*circle
    if remain<(2*b+np.sqrt(a**2+b**2)):
        ABCD=np.matmul(ABCD_free(remain),ABCD)
        return ABCD
    else:
        ABCD=np.matmul(ABCD_free(2*b+np.sqrt(a**2+b**2)),ABCD)
        ABCD=np.matmul(ACBD_con_mir(R),ABCD)
        ABCD=np.matmul(ABCD_free(remain-(2*b+np.sqrt(a**2+b**2))),ABCD)
        return ABCD
    
def q_evolve(q,ABCD):
    A=ABCD[0][0]
    B=ABCD[0][1]
    C=ABCD[1][0]
    D=ABCD[1][1]
    return (A*q+B)/(C*q+D)

def calc_waist_by_q(q):
    L=waveLength
    return np.sqrt(abs(-1*L/np.pi/np.imag(1/q)))
def calc_R_by_q(q):
    return 1/(np.real(1/q))

def Stable(ABCD):
    A=ABCD[0][0]
    D=ABCD[1][1]
    return (A+D)**2<=4


def calc_waist_cavity(Z):
    ABCD=ABCD_cavity(a,b,R,Z)
    q=q_evolve(q0,ABCD)
    return calc_waist_by_q(q)


def min_pos(Z_list,waist_list):
    minwaist=1000
    minindex=-1
    for i in range(0,len(Z_list)):
        if waist_list[i]<minwaist:
            minwaist=waist_list[i]
            minindex=i
    return minindex


waveLength=633*(10**-9)
w0=0.001
z0=0
R=0.36
q0=1j*np.pi*w0**2/waveLength
q0=q_evolve(q0,ABCD_free(1))
a=0.1
b=0.2
circle=2*b+2*(a**2+b**2)**0.5
Z_list=np.linspace(1,circle,1000)
waist_list=[]
for z in Z_list:
    waist_list.append(calc_waist_cavity(z))
    
print(f"The minimus position is at {Z_list[min_pos(Z_list,waist_list)]}")
plt.plot(Z_list,waist_list)
plt.xlim(0.86,0.88)
plt.ylim(0,0.0002)
    
