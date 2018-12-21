import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy import optimize
from math import pi

global K, n, Gamma, rho0_c

def RK4Step(u,t,dt,rhs):
    n = len(u)
    up = np.zeros(n)
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)
    k1 = dt*rhs(u,t)
    k2 = dt*rhs(u + 0.5*k1, t + 0.5*dt)
    k3 = dt*rhs(u + 0.5*k2, t + 0.5*dt)
    k4 = dt*rhs(u + k3, t + dt)
    up = u + (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return up # next u

# TOV equation for polytropic EOS p(rho) = K rho0^Gamma
# BC: m(0) = 0, p(0) = p_c
def rho0fp(p):
    if p < 0.0:
        print('Negative p encountered: %.2e')%(p)
        p = abs(p)
    return (p/K)**(1.0/Gamma) + p/(Gamma - 1.0)

def TOVrhs(u,r):
    [m, p] = u
    rho = rho0fp(p)

    if r < 1e-3:
        #use Taylor series
        m_rrr = 4.0/3.0*pi * (rho0_c**2/(Gamma - 1.0)
                              + 4.0**(-1.0/Gamma)*(rho0_c)**(2.0/Gamma))
        m_rr = m_rrr * r
        m_r = m_rr * r
    else:
        m_r = m / r
        m_rr = m_r / r
        m_rrr = m_rr / r

    MRHS = 4.0*pi*r**2 * rho
    PRHS = - rho * m_rr * (1.0 + p/rho) * (1.0 + p/rho) * ((1.0 + 4.0*pi*p / m_rrr)
                                                           *1.0/(1.0 - 2.0*m_r))

    return np.array([MRHS, PRHS])

def tovint(rho0_c, dr):
    nsteps = 100000
    sol = np.zeros((nsteps,2))
    n=0
    r=0.0
    p_c = K*rho0_c**Gamma
    u = np.array([0.0, p_c])
    sol[0] = u
    while sol[n,1] > 1.0e-9 and n < nsteps-1:
        u = RK4Step(u,r,dr,TOVrhs)
        n+=1
        sol[n] = u
        r+=dr
    
    return (r, sol[:n,0], sol[:n,1], rho0fp(p_c)) # (R, M, P, rho_c)

def tovplot(r, m, p):
    ri = np.linspace(0.0,r,len(m))
    plt.figure(1)
    plt.subplot(211)
    plt.plot(ri,m,'b-')
    plt.xlabel('r')
    plt.ylabel('m(r)')
    plt.title("reference model")
    plt.subplot(212)
    plt.plot(ri,p,'r--')
    plt.xlabel('r')
    plt.ylabel('p(r)')
    plt.show()

if __name__ == '__main__':
    
    # reference model
    K = 100
    n = 1
    Gamma = 1.0 + 1.0/float(n)
    rho0_c = 1.28e-3
    M_solar_G_over_c_sq = 1.98892e30 * 6.6384e-11 / 299792458**2
    dr = 1e-3
    (r, m, p, rho_c) = tovint(rho0_c, dr)
    print('stopped at n = %d, r= %.4f km, m = %.4e M_solar and p = %.4e' % (len(m), r*M_solar_G_over_c_sq / 1000.0, m[-1], p[-1]))
    tovplot(r * M_solar_G_over_c_sq / 1000.0, m, p)
