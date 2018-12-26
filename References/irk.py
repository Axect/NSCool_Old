from mpmath import mp, mpf, matrix, lu_solve, factorial, polyroots, eps, sqrt

class irk(object):
    def __init__(self, rang=10, acc=100):
        # Generates the coefficients of Legendre polynomial of n-th order
        # acc is the number of decimal characters of the coefficients
        # self.cf is the list with coefficients
        self.rang = rang
        self.acc = mp.dps = acc
        cn = mpf(0.0)
        k = mpf(0)
        n = mpf(rang)
        m = mpf(n/2)
        cf = []
        for k in range(n+1):
            cn = (-1)**(n+k)*factorial(n+k) / (factorial(n-k)*factorial(k)*factorial(k))
            cf.append(cn)
        cf.reverse()

        # Generates the coefficients of the implicit Runge-Kutta scheme of Gauss-Legendre type
        # Gives back the cortege (r, b, a), the terms of which correspond to Butcher scheme
        #
        # r1 | a11 . . . а1n
        #  . | .          .
        #  . | .          .
        #  . | .          .
        # rn | an1 . . . ann
        # ---+--------------
        #    |  b1 . . .  bn
        
        self.r = polyroots(cf) # Root of Legendre polynomial
        A1 = matrix(rang)
        for j in range(n):
            for k in range(n):
                A1[k, j] = self.r[j] ** k
        
        bn = []
        for j in range(n):
            bn.append(mpf(1.0) / mpf(j+1)) # Uniformly divide
        B = matrix(bn)
        
        self.b = lu_solve(A1, B)
        self.a = matrix(rang)
        for i in range(1, n+1):
            A1 = matrix(rang)
            cil = []
            for l in range(1, n+1):
                cil.append(mpf(self.r[i-1])**l / mpf(l))
                for j in range(n):
                    A1[l-1, j] = self.r[j] ** (l-1)
                Cil = matrix(cil)
                an = lu_solve(A1, cil)
                for k in range(n):
                    self.a[i-1,k] = an[k]

    def init(self, f, t, h, initvalues):
        self.size = len(initvalues)
        self.f = f
        self.t = t
        self.h = h
        self.yb = matrix(initvalues)
        self.ks = matrix(self.size, self.rang)
        for k in range(self.size):
            for i in range(self.rang):
                self.ks[k,i] = self.r[i]
        self.tn = matrix(1, self.rang)
        for i in range(self.rang):
            self.tn[i] = t + hj*self.r[i]
        self.y = matrix(self.size, self.rang)
        for k in range(self.size):
            for i in range(self.rang):
                self.y[k,i] = self.yb[k]
                temp = mpf(0.0)
                for j in range(self.rang):
                    temp += self.a[i, j]*self.ks[k,j]
                self.y[k,i] += temp
        self.yn = matrix(self.yb)

    def iterate(self, tn, y, yn, ks):
        # Generates the coefficients of the implicit Runge-Kutta scheme for the given step
        # with the method of the simple iteration with an initial value, coninciding with the coefficients,
        # calculated at the previous step. At sufficiently small step this must
        # work. There exists such a value of the step, Under which convergence is guaranteed.
        # No automatic re-setup of the step is foreseen in this procedure.
        mp.dps = self.acc
        y0 = matrix(yn)
        norme = mpf(1.0)
        eps0 = sqrt(eps)

