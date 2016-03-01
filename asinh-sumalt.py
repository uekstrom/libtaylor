from math import *
import cheb_expand

def compute_fraction(a,b):
    """Compute a (continued) fraction with coefficients a and b. a[0] is not used.
    F = b[0] + a[1]/(b[1] + a[2]/(..))
    """
    AA = [1.0, b[0]]
    BB = [0.0, 1.0 ]
    for j in range(1,len(a)):
        AA.append(b[j]*AA[-1] + a[j]*AA[-2])
        BB.append(b[j]*BB[-1] + a[j]*BB[-2])
    return AA[-1]/BB[-1]



def F(b,c,z,n=10):
    """Return an approximation to 2F1(1,b,c,z) evaluated using
    an n level partial fraction"""
    # In the notation and numbering of the wiki page:
    # b = [0,1,c,c+1,c+2,..] 
    # a = [0,1,-bz, (b-c)*z, -c(b+1)*z, 2*(b-c-1)*z, .. -(c+k)*(b+k+1)*z, (b-c-k)*(k+1)*z ..]
    #                k = 0     k = 0     k = 1
    # i:   0 1   2     3         4         5
    bs = [0.0]*n
    bs[1] = 1.0
    for k in range(2,n):
        bs[k] = c + k - 2
    a = [0.0]*n
    a[1] = 1.0
    a[2] = -b*z
    for k in range(n/2):
        i = 2*k + 3
        if i >= n:
            break
        a[i] = (b-c-k)*(k+1)*z
        i += 1
        if i >= n:
            break
        a[i] = -(c+k)*(b+k+1)*z
    return compute_fraction(a,bs)


def F_tail(b,c,z,m=5,n=10):
    """Return an approximation to 2F1(1,b,c,z) evaluated using
    an n level partial fraction"""
    # In the notation and numbering of the wiki page:
    # b = [0,1,c,c+1,c+2,..] 
    # a = [0,1,-bz, (b-c)*z, -c(b+1)*z, 2*(b-c-1)*z, .. -(c+k)*(b+k+1)*z, (b-c-k)*(k+1)*z ..]
    #                k = 0     k = 0     k = 1
    # i:   0 1   2     3         4         5
    bs = [0.0]*n
    bs[1] = 1.0
    for k in range(2,n):
        bs[k] = c + k - 2
    a = [0.0]*n
    a[1] = 1.0
    a[2] = -b*z
    for k in range(n/2):
        i = 2*k + 3
        if i >= n:
            break
        a[i] = (b-c-k)*(k+1)*z
        i += 1
        if i >= n:
            break
        a[i] = -(c+k)*(b+k+1)*z
    a = a[m:]
    b = bs[m:]
    b[0] = 0
    return compute_fraction(a,b)


def sumalt(x,e=1.0):
    n = len(x)
#    d = (3+sqrt(8))**n
    if e == 'auto':
        if x[-1] != 0:
            e = abs(x[-2]/x[-1])
        else:
            e = 1.0
    d = (2*e+1+2*sqrt(e*(e+1)))**n
    d = (d+1/d)/2
    b = 1.0
    c = d
    s = 0.0
    for k in range(n):
        c -= b
        s += c*x[k]
        b *= (2*(n+k)*(n-k)) / float((2*k+1)*(k+1))
    return s/d

def hyper_taylor(a,b,c,x,n):
    """Return the terms of the taylor expansion of F(a,b;c;x) at x=0"""
    c = float(c)
    terms = [1.0]
    for k in range(1,n):
        terms.append(terms[-1]*(a+k-1)*(b+k-1)*x/((c+k-1)*k))
    return terms


def sumalt_f(z,nderiv,n):
    """sumalt(hyper_taylor(1,1,3/2+n,-z))"""
    d = (3+sqrt(8))**n
    e = 1/(z+1e-8)
    #e = 1
    d = (2*e+1+2*sqrt(e*(e+1)))**n
    d = (d+1/d)/2 # for fixed n and e the constant d can be precomputed
    # can use d = d/2 as a good approximation
    b = 1.0
    c = d
    s = 0.0
    x = 1.0
    for k in range(n):
        c -= b
        s += c*x
        b *= (2*(n+k)*(n-k)) / float((2*k+1)*(k+1)) * e
        x *= -(k+1)/(1.5 + nderiv + k)*z
    return s/d


def asinh_sqrtx_div_sqrtx_smallx(x,nderiv,nterms=20):
    # Return the derivatives (not taylor coeffs) of 
    # the function asinh(sqrt(x))/sqrt(x) at x.
    # With 20 terms the fifth derivative is accurate
    # to 4e-17 at x = 1.
    # Using 11 terms the accuracy matches the largex function at x=1/2

    if nderiv == 0:
        vals = [sqrt(1+x)*sumalt_f(x,0,nterms)]
    else:
        # Compute the h(n)(x) functions
        vals = [1.0/sqrt(1+x)]
        for n in range(1,nderiv):
            vals.append(-0.5*(2*n-1)/(1+x)*vals[-1])
        # And the highest element
        vals.append(-0.5*(2*nderiv-1.0)/(2*nderiv+1.0)*vals[-1]*
                     sumalt_f(x,nderiv,nterms))
        # Recurse downwards
        for n in range(nderiv-1,-1,-1):
            vals[n] = (vals[n]-2*x*vals[n+1])/(2*n+1)
    return vals

def asinh_sqrtx_div_sqrtx_smallx_noaccel(x,nderiv,nterms=20):
    # Return the derivatives (not taylor coeffs) of 
    # the function asinh(sqrt(x))/sqrt(x) at x.
    # With 20 terms the fifth derivative is accurate
    # to 4e-17 at x = 1.
    # Using 11 terms the accuracy matches the largex function at x=1/2

    if nderiv == 0:
        vals = [sqrt(1+x)*sum(hyper_taylor(1.0,1.0,1.5+nderiv,-x,nterms))]
    else:
        # Compute the h(n)(x) functions
        vals = [1.0/sqrt(1+x)]
        for n in range(1,nderiv):
            vals.append(-0.5*(2*n-1)/(1+x)*vals[-1])
        # And the highest element
        vals.append(-0.5*(2*nderiv-1.0)/(2*nderiv+1.0)*vals[-1]*
                     sum(hyper_taylor(1.0,1.0,1.5+nderiv,-x,nterms)))
        # Recurse downwards
        for n in range(nderiv-1,-1,-1):
            vals[n] = (vals[n]-2*x*vals[n+1])/(2*n+1)
    return vals


def hn(x,nderiv):
    # Compute the h(n)(x) functions
    vals = [1.0/sqrt(1+x)]
    for n in range(1,nderiv):
        vals.append(-0.5*(2*n-1)/(1+x)*vals[-1])
    return vals


def asinh_sqrtx_div_sqrtx_smallx_fake(x,nderiv,fn):
    # Return the derivatives (not taylor coeffs) of 
    # the function asinh(sqrt(x))/sqrt(x) at x.
    # fn is the n:th derivative at x (must be given)

    if nderiv == 0:
        vals = [sqrt(1+x)*sum(hyper_taylor(1.0,1.0,1.5+nderiv,-x,nterms))]
    else:
        # Compute the h(n)(x) functions
        vals = [1.0/sqrt(1+x)]
        for n in range(1,nderiv):
            vals.append(-0.5*(2*n-1)/(1+x)*vals[-1])
        # And the highest element
        vals.append(fn)
        # Recurse downwards
        for n in range(nderiv-1,-1,-1):
            vals[n] = (vals[n]-2*x*vals[n+1])/(2*n+1)
    return vals





def asinh_sqrtx_div_sqrtx_largex(x,nderiv):
    """Somewhat stable at x=1/2, but best at x >= 1"""
    s = sqrt(x)
    t = sqrt(1+x)
    vals = [log(s+t)/s]
    hn = 1/t
    for n in range(nderiv):
        vals.append((hn-(2*n+1)*vals[-1])/(2*x))
        hn *= -0.5*(2*n+1)/(1+x)
    return vals


def asinh_sqrtx_div_sqrtx_updown_simulation(x,nderiv):
    fn = asinh_sqrtx_div_sqrtx_largex(x,nderiv)[-1]
    return asinh_sqrtx_div_sqrtx_smallx_fake(x,nderiv,fn)


def hyper_taylor(a,b,c,x,n):
    """Return the terms of the taylor expansion of F(a,b;c;x) at x=0"""
    c = float(c)
    terms = [1.0]
    for k in range(1,n):
        terms.append(terms[-1]*(a+k-1)*(b+k-1)*x/((c+k-1)*k))
    return terms    

def asinh_sqrtx_div_sqrtx_updown(x,nderiv,xc=0.1,nterms=12):
    """
    With parameters selected for the cut from upward to downward recursion at x = xc = 0.01.
    This value of x is only appropriate up to order 4 (abs error 1e-6). For high orders
    the x cut must be larger (and number of Taylor terms a bit larger).
    It is possible to use a smaller number of Taylor terms for higher derivatives,
    but not done to simplify the code.

    Possible reasonable combinations:
    xc = 0.01 nterms = 6 (up to nderiv=4)
    xc = 0.1  nterms = 12 (up to nderiv=7-9)
    xc = 0.5  nterms = 40 (~16 terms with sumalt)
    xc = 1.0  nterms = oo (20 terms with sumalt, not included here.)
    """
    if x < 1e-40:
        x = 1e-40
    s = sqrt(x)
    t = sqrt(1+x)
    fn = log(s+t)/s
    vals = [0]*(nderiv+1)
    vals[0] = fn
    if nderiv < 1:
        return vals
    # Compute hn using upward recursion
    hn = 1/t
    for n in range(nderiv-1):
        fn = (hn-(2*n+1)*fn)/(2*x)
        hn *= -0.5*(2*n+1)/(1+x)
        vals[n+1] = hn
    # Final upward recursion fn
    fn_up = (hn-(2*nderiv-1)*fn)/(2*x)
    # Compute the nterms taylor approximation to fn.     
    # The number of Taylor terms can depend on nderiv, but not on x.
    fn_taylor = 1.0
    t = 1.0 # Reusing the t array for terms in the sum
    for k in range(1,nterms):
        t *= -k/(0.5+nderiv+k)*x
        fn_taylor += t
    if nderiv > 1:
        fn_taylor *= -0.5*(2*nderiv-1.0)/(2*nderiv+1.0)*vals[-2]
    else:
        fn_taylor *= -0.5/3.0*hn
        
    # Determine (using x) if we should use fn (upward value) or Taylor approx for fn.
    if x < xc:
        vals[-1] = fn_taylor
    else:
        vals[-1] = fn_up
    # Recurse downwards but stop before vals[0] is modified.
    for n in range(nderiv-1,0,-1):
        vals[n] = (vals[n]-2*x*vals[n+1])/(2*n+1)
    return vals


def test_cf():
    n = 0
    m = 16
    print "Continued fraction"
    print "%i terms, x, abs err, rel err" % m
    x = 0.0
    while x<=1.0:
        ref = F(1,n+1.5,-x,30)
        val = F(1,n+1.5,-x,m)
        print x,val-ref,(val-ref)/ref
        x += 0.1
    print "Taylor+sumalt"
    x = 0.0
    while x<=1.0:
        ref = F(1,n+1.5,-x,30)
        val = sumalt(hyper_taylor(1,1,n+1.5,-x,m))
        val = sumalt_f(x,n,m)
        print x,val-ref,(val-ref)/ref
        x += 0.1

def test_final():
    """
    Small x cuts may also be acceptable, i.e.

    At x = 0.01
    Errors per order, upwards, downwards (4 term Taylor)
    0 -1.11e-16  0.00e+00
    1  4.19e-15  1.39e-15
    2 -6.27e-13 -2.07e-13
    3  1.57e-10  5.17e-11

    At x = 0.001
    Errors per order, upwards, downwards
    0  2.89e-15 0.00e+00
    1 -1.42e-12 0.00e+00
    2  2.14e-09 0.00e+00
    3 -5.34e-06 5.27e-15

    At x = 0.0001
    Errors per order, upwards, downwards
    0 -9.10e-15 0.00e+00
    1  4.53e-11 0.00e+00
    2 -6.80e-07 0.00e+00
    3  1.70e-02 0.00e+00 [3 or even 2 terms are enough (2e-10 error)]


    Accurate/high order version:
    x < 1  : 20 term taylor with sumalt
    x >= 1 : Upward recursion. (loses 1 digit per order at x = 1)

    Fast version up to order ~ 4-6 ok when starting at high order and recurse down:
    x < 0.01 : 7 term direct Taylor (as accurate as upwards for low orders)
    x >= 0.01: Upwards recursion.

    Adaptive version: Replace less accurate upward derivatives with short Taylor.
    Is max order independent errors important?

    Example errors 3 term Taylor (T3)
    x = 0.0001 - T3 from first derivative
    x = 0.001  - T3 from second
    x = 0.01   - T3 downwards from 5 is almost the same as upward (and cheaper)
    x = 0.1    - T3 downwards from 9 is better than upward (but none very good at order 9..)
    """
    x = 0.005
    print hn(x,5)
    print asinh_sqrtx_div_sqrtx_largex(x,5)
    print asinh_sqrtx_div_sqrtx_updown(x,5)

#     x = 0.5
#     ref_05_5 = -4.0347461078459588001591671458e-1
#     ref_1_5 = -1.0403759000548234568125678246e-1
#     ref_05_3 = -8.8006016019340102483799244871e-2
#     ref_05_2 = 7.2034308599049756012255050634e-2
#     print "At x =",x
#     print "Small x orders 2,3,5"
#     nterms = 11
#     print "Derivatives",asinh_sqrtx_div_sqrtx_smallx(x,5,nterms=nterms)[2]-ref_05_2
#     print "Derivatives",asinh_sqrtx_div_sqrtx_smallx(x,5,nterms=nterms)[3]-ref_05_3
#     print "Derivatives",asinh_sqrtx_div_sqrtx_smallx(x,5,nterms=nterms)[5]-ref_05_5
#     print "Derivatives",asinh_sqrtx_div_sqrtx_smallx(x,10,nterms=20)[10]
#     print "Large x orders 2,3,5"
#     print "Derivatives",asinh_sqrtx_div_sqrtx_largex(x,5)[2]-ref_05_2
#     print "Derivatives",asinh_sqrtx_div_sqrtx_largex(x,5)[3]-ref_05_3
#     print "Derivatives",asinh_sqrtx_div_sqrtx_largex(x,5)[5]-ref_05_5
#     print "Derivatives",asinh_sqrtx_div_sqrtx_largex(x,10)[10]

    x = 0.099
    print "At x =",x
    print "Errors per order, upwards, updown"
    nderiv = 10
    ref = asinh_sqrtx_div_sqrtx_smallx(x,nderiv,nterms=20)
    lx = asinh_sqrtx_div_sqrtx_largex(x,nderiv)
    udx = asinh_sqrtx_div_sqrtx_updown(x,nderiv)
    for k in range(len(lx)):
        print k,"%.2e %.2e" % (lx[k] - ref[k],udx[k] - ref[k])

#asinh_sqrtx_div_sqrtx_smallx_noaccel(x,k,nterms=nterms)[-1] - ref[k])


test_final()

#print sum(c)
#print "%.16e" % sumalt(c)
def test():
    ref_2 = 7.978381361379436e-1
    ref_5 = 8.768718309484423e-1
    ref_2_01 = 9.72633145557962e-1
    ref_5_01 = 9.850117434593606e-1
    print sumalt(hyper_taylor(1,1,2+1.5,-1,5))-ref_2
    print sumalt(hyper_taylor(1,1,2+1.5,-1,10))-ref_2
    print sumalt(hyper_taylor(1,1,2+1.5,-1,15))-ref_2
    print sumalt(hyper_taylor(1,1,2+1.5,-1,18))-ref_2
    print
    print sumalt(hyper_taylor(1,1,5+1.5,-1,5))-ref_5
    print sumalt(hyper_taylor(1,1,5+1.5,-1,10))-ref_5
    print sumalt(hyper_taylor(1,1,5+1.5,-1,15))-ref_5
    print sumalt(hyper_taylor(1,1,5+1.5,-1,18))-ref_5
    print "Fraction"
    print F(1,5+1.5,-1,5)-ref_5
    print F(1,5+1.5,-1,10)-ref_5
    print F(1,5+1.5,-1,15)-ref_5
    print F(1,5+1.5,-1,18)-ref_5

def print_err(n):
    x = -.2
    while x <= 1.1:
#        print "%.15e %.15e" % (x, F_tail(1,n+1.5,-x,5,30))
        print "%.15e %.15e" % (x,F(1,n+1.5,-x,30))
#        c = hyper_taylor(1,1,n+1.5,-x,15)
#        print x,c[-1]/(1e-15+c[-2])
        x += 0.01


def test_chebfit(n):
    def f(x):
        return F(1,n+1.5,-(x+1.0)/2,30)
    c = cheb_expand.cheb_coeffs(f,15)
    print "Chebyshev coeffs at n = %i" % n, c

#test_chebfit(1)
#test_chebfit(5)
#test_chebfit(10)

#F(1,3+1.5,-1.0,10)


def F_aitken(b,c,z,n=10):
    Fm2 = F(b,c,z,n-4)
    Fm1 = F(b,c,z,n-2)
    F0 = F(b,c,z,n)
    return (F0*Fm2 - Fm1**2)/(F0 - 2*Fm1 + Fm2)


def test_aitken():
    ref_2 = 7.978381361379436e-1
    ref_5 = 8.768718309484423e-1
    ref_2_01 = 9.72633145557962e-1
    ref_5_01 = 9.850117434593606e-1
    print "Standard continued Fraction"
    print F(1,5+1.5,-1,10)-ref_5
    print F(1,5+1.5,-1,15)-ref_5
    print F(1,5+1.5,-1,18)-ref_5
    print "Aitken"
    print F_aitken(1,5+1.5,-1,10)-ref_5
    print F_aitken(1,5+1.5,-1,15)-ref_5
    print F_aitken(1,5+1.5,-1,18)-ref_5


def test_rates():
    n = 5
    for x in [0.5,1.0]:
        print "Rate of convergence at x=",x,", error"
        ref = F(1,n+1.5,-x,30)
        for k in range(8,16):
            print k, (F(1,n+1.5,-x,k)-ref)/(F(1,n+1.5,-x,k-2)-ref),(F(1,n+1.5,-x,k)-ref)


