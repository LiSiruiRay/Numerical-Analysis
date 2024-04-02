from pylab import *
import sys


def bisection_original(f, a, b, tol=1e-10):
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
        c = b
        b = a
        a = c
    if a == b:
        raise ValueError("Bisection called with a==b\n")
    if f(a) * f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))

    counter = 0
    while (b - a) > tol:
        print(f"a={repr(a)},   b={repr(b)},   b-a={repr(b - a)}")
        midpoint = (a + b) / 2.
        sfm = sign(f(midpoint))
        if sfm == sign(f(a)):
            a = midpoint
        elif sfm == sign(f(b)):
            b = midpoint
        elif sfm == 0:
            # lucky case: we found an exact zero!
            print("Found an exact(?!?) zero!")
            return midpoint
        else:
            raise ValueError("Something went horribly bad: sign(f(midpoint))={}\n".format(sfm))
        counter += 1

    print(f"bisection_original, counter={counter}")
    return (a + b) / 2.


def bisection(f, a, b, tol=1e-10, relative_tol=sys.float_info.epsilon):
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
        c = b
        b = a
        a = c
    if a == b:
        raise ValueError("Bisection called with a==b\n")
    if f(a) * f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))
    print(f"check (b - a): {(b - a)}, (b - a) / a: {(b - a) / a}, (b - a) / a: {(b - a) / a}")
    counter = 0
    while (b - a) > tol and abs(b - a) / abs(a) > relative_tol and abs(b - a) / abs(a) > sys.float_info.epsilon:
        print(f"a={repr(a)},   b={repr(b)},   b-a={repr(b - a)}")
        midpoint = (a + b) / 2.
        sfm = sign(f(midpoint))
        if sfm == sign(f(a)):
            a = midpoint
        elif sfm == sign(f(b)):
            b = midpoint
        elif sfm == 0:
            # lucky case: we found an exact zero!
            print("Found an exact(?!?) zero!")
            return midpoint
        else:
            raise ValueError("Something went horribly bad: sign(f(midpoint))={}\n".format(sfm))
        counter += 1
    print(f"bisection, counter={counter}")
    return (a + b) / 2.


def Newton_original(f, fp, a, tol=1e-10):
    """Finds a zero of f by Newton's method

Input:
   f   function with only one float argument, that returns a float
  fp   the first derivative of f
   a   float

 tol   stop the iteration when |newa-a| < tol, where a, newa are extimates of the zero at consecutive iterations.

Output:
   x such that f(x) is approximately 0.
"""
    counter = 0
    # for i in range(10):
    while True:
        # print(f"loop {i}")
        newa = a - f(a) / fp(a)

        print(f"newa={repr(newa)}, a={a}, fp(a)={fp(a)}, f(a)={f(a)},  |newa-a|={repr(fabs(a - newa))}, fabs(a - newa) / fabs(newa) = {fabs(a - newa) / fabs(newa)}")
        if fabs(a-newa) < tol:
            break
        else:
            a = newa
        counter += 1

    print(f"Newton_original, counter={counter}")
    return newa


def Newton(f, fp, a, tol=1e-10, interval=None):
    """Finds a zero of f by Newton's method

Input:
   f   function with only one float argument, that returns a float
  fp   the first derivative of f
   a   float

 tol   stop the iteration when |newa-a| < tol, where a, newa are extimates of the zero at consecutive iterations.

Output:
   x such that f(x) is approximately 0.
"""
    if interval is None:
        interval = [a - 1e2, a + 1e2]
    x_start = interval[0]
    x_end = interval[1]
    f_start = f(x_start)
    f_end = f(x_end)
    print(f"check x_start: {x_start}, x_end: {x_end}")

    if f_start * f_end >= 0.:
        raise ValueError("interval should have different sign")

    # for i in range(500):
    counter = 0

    while True:
        newa = a - f(a) / fp(a)
        print(f"newa={repr(newa)},  |newa-a|={repr(fabs(a - newa))}")

        if fabs(a - newa) / fabs(newa) < tol:
            break

        # bisection step
        if newa < x_start or newa > x_end or abs(newa - a) > 0.5 * abs(x_end - x_start):
            newa = (x_end + x_start) / 2
            if sign(f(newa)) == sign(f(x_start)):
                x_start = newa
            elif sign(f(newa)) == sign(f(x_end)):
                x_end = newa
            elif sign(f(newa)) == 0:
                # exact zero point
                break
            a = newa
        else:
            a = newa
            if f(a) * f(x_start) < 0:
                x_end = newa
            else:
                x_start = newa
        counter += 1
    print(f"Newton, counter={counter}")
    return newa


def x_e_n_x(x: float):
    return x * exp(-x)


def derivative_xenx(x: float):
    return exp(-x) - x * exp(-x)


def regula_falsi(f,
                 a: float, b: float,
                 tol: float = 1e-10,
                 relative_tol=sys.float_info.epsilon) -> float:
    if a > b:
        c = b
        b = a
        a = c
    if a == b:
        raise ValueError("Bisection called with a==b\n")
    if f(a) * f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))

    counter = 0
    while True:
        # for i in range(50):
        c = (b * f(a) - a * f(b)) / (f(a) - f(b))
        print(f"a={repr(a)},   b={repr(b)},   c={repr(c)}")

        min_dis = min(abs(c - a), abs(b - c))
        # if min_dis < tol or (min_dis / abs(c) < relative_tol) or (min_dis / abs(c) < sys.float_info.epsilon):
        if min_dis < tol:
            break

        if sign(f(c)) == sign(f(b)):
            b = c
        elif sign(f(c)) == sign(f(a)):
            a = c
        elif sign(f(c)) == 0:
            break
        counter += 1

    print(f"regula_falsi, counter={counter}")
    return c


## Test code
if __name__ == "__main__":
    bisection_original(x_e_n_x, -0.5, 0.6, tol=1e-10)
    # # Newton(x_e_n_x, derivative_xenx, 500, interval=[-10, 550])
    regula_falsi(x_e_n_x, a=-0.5, b=0.6, tol=1e-10)
    Newton_original(x_e_n_x, derivative_xenx, a=-0.5)
    # Newton_original(sin, cos, 1.5)
    # Newton_original(sin, cos, 3)
    pass
