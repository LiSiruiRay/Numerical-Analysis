from scipy.linalg import lu


def create_matrix(n):
    """Returns an nxn diagonally-dominant matrix to be used for the
exercises.

Input:
    n  a positive integer

Output:
    an nxn matrix

Note: in the PLU decomposition of the output matrix, P is the identity matrix. 
"""
    n = int(n)
    if n<3:
        raise(ValueError("n must be at least 3."))
    A = empty((n,n))#3.*identity(n)
    v = zeros((n,))
    v[0]  = 3.
    v[1]  = 1.
    v[-1] = 1.
    for i in range(n):
        A[i] = roll(v, i)

    return A


#--------------------------------------------------------------------------


def create_matrix_2(n):
    """Returns an nxn matrix to be used for the exercises.

Input:
    n  a positive integer

Output:
    an nxn matrix

Note: in the PLU decomposition of the output matrix, P is not the identity matrix. The matrix returned by this function is not, in general, diagonally dominant.
"""
    A = create_matrix(n)
    primes = [3,5,7,11]
    lp = len(primes)
    #swap the rows in a non-trivial way
    for i in range(n):
        ip = i%lp
        j = (i+primes[ip])%n
        line = copy(A[i])
        A[i] = A[j]
        A[j] = line

    return A
        

#----------------------------------------------------------------        
    

    
