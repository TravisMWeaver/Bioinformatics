import sys
import string
import math

#-------------------------------------------
def gammln(xx):
    cof = [76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179E-2,-0.5395239384953E-5]
    x = xx
    y = xx
    tmp = x + 5.5
    tmp = (x+0.5)*math.log(tmp) - tmp
    ser = 1.000000000190015
    for i in range(0,6):
        y = y + 1.0
        ser = ser + cof[i]/y

    gammln = tmp+math.log(2.5066282746310005*ser/x)
    return gammln

#-------------------------------------------
def betai(a,b,x):
    betai =0.0
    if x < 0.0 or x > 1.0:
        sys.stderr.write("bad argument x in betai\n")
        sys.exit(0)
    if x == 0.0 or x == 1.0:
        bt = 0.0
    else:
        bt = math.exp(gammln(a+b) - gammln(a) - gammln(b) + a*math.log(x) + b*math.log(1.0 -x))

    if x < float(a +1.0)/(a+b +2.0):
        betai = bt*betacf(a,b,x)/a
        return betai
    else:
        betai = 1.0 - bt*betacf(b,a,1.0-x)/b
        return betai

#-------------------------------------------
def avevar (data):
    ave = 0.0
    N = len(data)
    for j in range(0, N):
        ave = ave + data[j]

    ave = ave/N

    var = 0.0
    ep = 0.0
    for j in range(0, N):
         s = data[j] -ave
         ep = ep + s
         var = var + s*s

    var = (var -ep*ep/N)/(N-1)

    return (ave, var)

#-------------------------------------------
def betacf(a,b,x):
    betacf = 0.0
    maxit = 100
    eps = 3.0e-7
    fpmin = 1.0e-30
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1.0
    d = 1.0 - float(qab*x)/qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0/float(d)
    h = d
    m = 1
    while m <= maxit:
        m2 = 2*m
        aa = float(m*(b-m)*x)/((qam + m2)*(a+m2))
        d = 1.0 + aa*d
        if abs(d) < fpmin:
            d = fpmin
            
        c = 1.0+float(aa)/c

        if abs(c) < fpmin:
            c = fpmin

        d = 1.0/float(d)

        h = h*d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d = 1.0+aa*d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + float(aa)/c

        if abs(c) < fpmin:
            c = fpmin
        d = 1.0/float(d)
        dEL = d*c

        h=h*dEL

        if abs(dEL - 1.0) < eps:
            betacf = h
            return betacf

    sys.stderr.write("a or b too big, or MAXIT too small in betacf\n")
    sys.exit(0)


#-----------------------------------
def tutest (x1, x2):
    N1 = len(x1)
    N2 = len(x2)
    
    AVEVAR1 = avevar(x1)
    ave1 = AVEVAR1[0]
    var1 = AVEVAR1[1]
    
    AVEVAR2 = avevar(x2)
    ave2 = AVEVAR2[0]
    var2 = AVEVAR2[1]
    
    df = N1 + N2 - 2
    
    var = ((N1 - 1)*var1+(N2 - 1)*var2)/df
    t = (ave1-ave2)/math.sqrt(var*(1.0/N1+1.0/N2))
    prob = betai(0.5*df,0.5, df/(df+t*t))
    
    return t, prob
