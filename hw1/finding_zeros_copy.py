from pylab import *

def bisection(f, a, b, tol = 1e-10):
    """Finds a zero of f by the bisection method

The sign of f(a) and f(b) must be opposite, and f must be a continuous function.
Input:
   f   function with only one float argument, that returns a float
   a   float, one end of the bracketing interval
   b   float, the other end of the bracketing interval

 tol   stop the iteration when (b-a)/2 < tol

Output:
   x such that f(x) is approximately 0.
"""
    # sanity checks to cover us from user distraction (n.b. usually user=me)
    if a > b:
        c=b
        b=a
        a=c
    if a==b:
        raise ValueError("Bisection called with a==b\n")
    if f(a)*f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))

    while (b-a) > tol:
        print(f"a={repr(a)},   b={repr(b)},   b-a={repr(b-a)}")
        midpoint = (a+b)/2.
        sfm = sign(f(midpoint))
        if sfm == sign(f(a)):
            a = midpoint
        elif sfm == sign(f(b)):
            b = midpoint
        elif sfm == 0:
            #lucky case: we found an exact zero!
            print("Found an exact(?!?) zero!")
            return midpoint
        else:
           raise ValueError("Something went horribly bad: sign(f(midpoint))={}\n".format(sfm))
    return (a+b)/2.

       

def Newton(f, fp, a, tol=1e-10):
    """Finds a zero of f by Newton's method

Input:
   f   function with only one float argument, that returns a float
  fp   the first derivative of f
   a   float

 tol   stop the iteration when |newa-a| < tol, where a, newa are extimates of the zero at consecutive iterations.

Output:
   x such that f(x) is approximately 0.
"""
    while True:
        newa = a - f(a)/fp(a)
        print(f"newa={repr(newa)},  |newa-a|={repr(fabs(a-newa))}" )
        if fabs(a-newa) < tol:
            break
        else:
            a = newa
    return newa



