from finding_zeros import bisection
from pylab import *

if __name__ == '__main__':
    bisection(cos, 1e6 * pi, (1.e6 + 1) * pi, tol=1e-9)
