
import math

R=3
lambdaa=1

Tol=0.1
e = 2**(0.5)
mp = 0.2**0.5
me = mp / 4
ma = mp
#ma = me / (pow(3, 2.0 / 3) * pow(pi, 4.0 / 3))
G = 1
h = 1
eps = (8 * me * e * e) / (h * h * 3**(2.0 / 3) * math.pi **(1.0 / 3))
eta = 4 * math.pi * e * e * 2 * ma / (h * h)

#first we solve the 4th order problem. This will be our initial guess. 
def f(t, ue, ua,  sa1, sa2, sa3): 
    if (t == 0):
        return -sa2 * (-sa2 / ua) - eta * (4 - G * ma*ma / e**2) * ua * ua * ua + eta * (2 + G * me * ma / e**2) * ue**(3.0 / 2) * ua
    elif (ua == 0):
        return 0
    else:
        return -sa3 * (4.0 / t - 2 * sa1 / (ua)) - sa2 * (-8.0 * sa1 / (t * ua) - sa2 / (ua)+2 * sa1 * sa1 / (ua * ua)) - sa1 * (4 * sa1 * sa1 / (t * ua * ua)) - eta * (4 - G * ma*ma / (e*e)) * ua * ua * ua + eta * (2 + G * me * ma / (e*e)) * ue**(3.0 / 2) * ua
        

def g(t,ue, ua, se):
    if (ue == 0):
        return 0
    elif (t == 0): 
        return -(2 + G * me * ma / (e*e)) * eps * ua * ua + eps * (1 - G * me*me / (e*e)) * ue**(3.0 / 2)
    else:
        return -(2.0 / t) * se - (2 + G * me * ma / (e*e)) * eps * ua * ua + eps * (1 - G * me*me / (e*e)) * (ue**(3.0 / 2)) 

def F( r1, r2, n):
    if (r1 < r2):
        return r1**(n + 2) / r2**(n + 1)
    elif (r2 < r1):
        return r2**n / r1**(n - 1)
    else:
        return r1

#so the initial list will have rstep of rmax/2048
#then we iterate again through this

t = 0
ue = 0.17**(2.0/3)
ua = 1
step =0.001
se = 0
sa =0
i = 0
aList=[1]
eList=[0.17**(2.0/3)]
#this first pass is to find the size of the initial guess
while (i < 6001):

    if (ue + step * se <= 0):
        ue = 0
        se = 0

    elif (ua + step * sa1 <= 0):
        ua = 0
        sa1 = 0

    nsa3 = sa3 + step * f(t, ue, ua, sa1, sa2, sa3)
    nsa2 = sa2 + step * sa3
    nsa1 = sa1 + step * sa2
    nse = se + step * g(t, ue, ua, se)

    nua = ua + step * sa1
    nue = ue + step * se
    sa3 = nsa3
    sa2 = nsa2
    sa1 = nsa1
    ua = nua
    se = nse
    ue = nue

    if (ua <= 0):
        ua = 0
        sa1 = 0
        sa2 = 0
        sa3 = 0
    
    i += 1
    aList.append(ua)
    eList.append(ue)
    t = t + step

#at this point should have initial guess, so we can start the iteration
#we first need to make the electric and gravitational grids    
GravityGrid=[]
ElectricGrid=[]

for j in range(6001):
    D_G = 0
    D_E=0
    r1=j*0.001
    for i in range(6001):
        sum=0
        r2=i*0.001

        for n in range(25):
            sum+=F(r2, r1, 2*n)
        if (i == 0 or i == 6000):
            D_G += 7 *(me*eList[i]**(3/2)+ma*aList[i]**2)
            D_E += 7 *(2*aList[i]**2-eList[i]**(3/2))
        elif (i%4 == 0) :
            D_G += 14*(me*eList[i]**(3/2)+ma*aList[i]**2)
            D_E += 14*(2*aList[i]**2-eList[i]**(3/2))
        elif (i% 4 == 2): 
            D_G += 12*(me*eList[i]**(3/2)+ma*aList[i]**2)
            D_E += 12*(2*aList[i]**2-eList[i]**(3/2))
        elif (i% 2 == 1 or i% 2 == -1) :
            D_G += 32*(me*eList[i]**(3/2)+ma*aList[i]**2)
            D_E += 32*(2*aList[i]**2-eList[i]**(3/2))
    GravityGrid.append(2*math.pi*0.002/(6)*D_G)
    ElectricGrid.append(2*math.pi*0.002/(6)*D_E)

#now we should have the grids calculated. So we calculate the new densities
neweList=[]

#for the electrons, we first caclulate the lagrange mulitplier
#convert R to the nearest grid point
#there are many errors past this point

eRad=R/0.001
newlambdae=me*G(eRad)+e*e*E(eRad)
for i in range(6001):
    neweList.append(max((-newlambdae+me*G(i)+e*e*E(i))*(2*me/(h*h)*3**(2/3)*math.pi**(4/3)),0))

#now for the alpha particles, we need to use an ODE solver

ua=1
newaList=[1]
i=1
sa=step*(2*ma/(h*h)*(-ma*G(0)*ua+2*e*e*E(0)*ua))
t=0.001

while (i < 6002):
    if (ua + step * sa1 <= 0):
        ua = 0
        sa=0
    nsa = sa + step *((2/t)*sa+2*ma/(h*h)*(-ma*G(i)*ua+2*e*e*E(i)*ua))
    nua = ua + step *sa
    ua=nua
    sa=nsa
    if (ua <= 0):
        ua = 0
        sa=0
    
   
    newaList.append(ua)
    i += 1
    t = t + step

#this should be one iteration.


lambdadiff=newlambdae
adiff=0
ediff=0
amax=0
emax=0
for i in range(6001):
    if abs(aList[i]-newaList[i])>adiff:
        adiff=abs(aList[i]-newaList[i])
    if abs(eList[i]-neweList[i])>ediff:
        ediff=abs(eList[i]-neweList[i])
    if aList[i]>amax:
        amax=newaList[i]
    if eList[i]>emax:
        emax=neweList[i]
iter=1
while (iter<10 and (adiff/amax>Tol or ediff/emax>Tol or abs(lambdadiff)>Tol)):
    iter+=1
    print(iter)
    eList=[]
    aList=[]
    lambdae=newlambdae
    for i in range(6001):
        eList.append(neweList[i])
        aList.append(newaList[i])

    GravityGrid=[]
    ElectricGrid=[]

    for j in range(6001):
        D_G = 0
        D_E=0
        r1=j*0.001
        for i in range(6001):
            sum=0
            r2=i*0.001

            for n in range(25):
                sum+=F(r2, r1, 2*n)
            if (i == 0 or i == 6000):
                D_G += 7 *(me*eList[i]**(3/2)+ma*aList[i]**2)
                D_E += 7 *(2*aList[i]**2-eList[i]**(3/2))
            elif (i%4 == 0) :
                D_G += 14*(me*eList[i]**(3/2)+ma*aList[i]**2)
                D_E += 14*(2*aList[i]**2-eList[i]**(3/2))
            elif (i% 4 == 2): 
                D_G += 12*(me*eList[i]**(3/2)+ma*aList[i]**2)
                D_E += 12*(2*aList[i]**2-eList[i]**(3/2))
            elif (i% 2 == 1 or i% 2 == -1) :
                D_G += 32*(me*eList[i]**(3/2)+ma*aList[i]**2)
                D_E += 32*(2*aList[i]**2-eList[i]**(3/2))
        GravityGrid.append(2*math.pi*0.002/(6)*D_G)
        ElectricGrid.append(2*math.pi*0.002/(6)*D_E)

    neweList=[]

    #there are many errors past this point
    newlambdae=me*G(eRad)+e*e*E(eRad)
    for i in range(6001):
        neweList.append(max((-newlambdae+me*G(i)+e*e*E(i))*(2*me/(h*h)*3**(2/3)*math.pi**(4/3)),0))

    #now for the alpha particles, we need to use an ODE solver

    ua=1
    newaList=[1]
    i=1
    sa=step *(2*ma/(h*h)*(-ma*G(0)*ua+2*e*e*E(0)*ua))
    t=0.001

    while (i < 6002):

        if (ua + step * sa1 <= 0):
            ua = 0
            sa=0

        nsa = sa + step *((2/t)*sa+2*ma/(h*h)*(-ma*G(i)*ua+2*e*e*E(i)*ua))
        nua = ua + step *sa

        ua=nua
        sa=nsa
        if (ua <= 0):
            ua = 0
            sa=0
    
        i += 1
        aList.append(ua)
        t = t + step

    lambdadiff=abs((lambdae-newlambdae)/newlambdae)
    adiff=0
    ediff=0
    amax=0
    emax=0
    for i in range(6001):
        if abs(aList[i]-newaList[i])>adiff:
            adiff=abs(aList[i]-newaList[i])
        if abs(eList[i]-neweList[i])>ediff:
            ediff=abs(eList[i]-neweList[i])
        if aList[i]>amax:
            amax=newaList[i]
        if eList[i]>emax:
            emax=neweList[i]

