"""iod.py

Author: Michael Roberts-Tsoukkas
Created: 4/14/2025
Description: Contains implementation of Gaussian and Gibbs IOD algorithms."""

# Python Standard Imports

# Third Party Imports
import numpy as np

# Local Imports
from ..constants import MU_SUN, EARTH_SMA

def iod(L:np.ndarray,
        ROG:np.ndarray,
        t1:float,
        t3:float,
        mu:float=MU_SUN,
        ) -> tuple[np.ndarray]:
    """Gaussian and Gibbs IOD.
    
    Args:
        L (np.ndarray): Matrix of unit vectors [L1 L2 L3] directed from observer to target,
        ROG (np.ndarray): Matrix of position vectors of the observer relative to the origin.
        t1 (float): Time of first observation relative to the second observation.
        t3 (float): Time of third observation relative to the second observation.
        mu (float, Optional): Specific Gravitational Parameter for the system. Defaults to MU_Sun

    Returns:
        tuple[np.ndarray]: Position, velocity vectors in the intertial frame centered about the origin.
    """

    # GAUSS METHOD -------------------------------------
    # Computing the a coefficients
    a1 = t3/(t3-t1) 
    a1u = t3*((t3-t1)**2 - t3**2)/(6*(t3-t1))
    a3 = -t1/(t3-t1)
    a3u = -(t1*((t3-t1)**2 - t1**2))/(6*(t3-t1))  

    # Computing the M matrix
    Linv = np.linalg.inv(L)    
    M = Linv*ROG

    # Computing the d coefficients
    d1 = M[1,0]*a1 - M[1,1] + M[1,2]*a3   
    d2 = M[1,0]*a1u + M[1,2]*a3u

    # Extracting information on the second observation
    L2 = L[:,1] 
    ROG2 = ROG[:,2]  

    # C = dot(L2,ROG2) essentially
    # C = (np.transpose(L2)*ROG2)[0,0] # Michael original stuff
    C = np.dot(L2, ROG2)

    # Computing r2 polynomial coefficients
    p6 = -(d1**2 + 2*C*d1 + (np.linalg.norm(ROG2))**2)
    p3 = -2*mu*(C*d2 + d1*d2)
    p0 = -mu**2*d2**2

    # Solving the polynomial and extracting the (assumedly) 1 physical solution
    P = [1,0,p6,0,0,p3,0,0,p0]
    r2 = np.roots(P)
    r2 = r2[np.isreal(r2)]
    r2 = np.real(min(r2))

    # Computing the u parameter
    u = -mu/(r2**3)

    # Computing the vector of c coefficients
    Cvec = np.matrix([a1+a1u*u, -1, a3+a3u*u])
    c1 = Cvec[0,0]
    c2 = Cvec[0,1]
    c3 = Cvec[0,2]

    # Solving for the rho magnitudes
    X = -M*np.transpose(Cvec)
    rho1 = 1/c1*X[0,0]
    rho2 = 1/c2*X[1,0]
    rho3 = 1/c3*X[2,0]

    # Solving for all three position vectors
    r1 = rho1*L[:,0] + ROG[:,0]
    r2 = rho2*L[:,1] + ROG[:,1]
    r3 = rho3*L[:,2] + ROG[:,2]
    print(f"Position vectors solved!")
    print(f"r1 = {r1}")
    print(f"r2 = {r2}")
    print(f"r3 = {r3}")
    print(f"|r1| = {np.linalg.norm(r1)/EARTH_SMA} AU")
    print(f"|r2| = {np.linalg.norm(r2)/EARTH_SMA} AU")
    print(f"|r3| = {np.linalg.norm(r3)/EARTH_SMA} AU")
    # END OF GAUSS METHOD -------------------------------------------------------

    # GIBBS METHOD -----------------------------------------------------------
    # Computing all three position vector magnitudes
    r1m = np.linalg.norm(r1)
    r2m = np.linalg.norm(r2)
    r3m = np.linalg.norm(r3)

    # Computing the Z vectors
    Z12 = np.cross(np.transpose(r1),np.transpose(r2))
    Z23 = np.cross(np.transpose(r2),np.transpose(r3))
    Z31 = np.cross(np.transpose(r3),np.transpose(r1))

    # Computing the N,D,S vectors
    N = r1m*Z23 + r2m*Z31 + r3m*Z12
    D = Z12 + Z23 + Z31
    S = (r2m-r3m)*r1 + (r3m-r1m)*r2 + (r1m-r2m)*r3

    # Computing the B vector
    B = np.cross(D,np.transpose(r2))

    # Computing the N and D magnitudes
    Nm = np.linalg.norm(N)
    Dm = np.linalg.norm(D)

    # Computing the Lg weight
    Lg = np.sqrt(mu/(Nm*Dm))

    # Computing v2
    v2 = Lg/r2m*np.transpose(B) + Lg*S

    # END OF GIBBS METHOD -----------------------------------------------------
    return r2, v2

def iod2(t1:float,
         t3:float,
         obs_pos:np.ndarray,
         rho_vec:np.ndarray,
         mu:float=MU_SUN,
         ) -> np.ndarray:
    """Calculates the position and velocity of an object orbiting the sun using Gaussian and Gibbs Orbit Determination.

    Args:
        t1 (float): Time difference between the second and fist observations, in seconds.
        t3 (float): Time difference between the third and second observations, in seconds.
        obs_pos (np.ndarray): List of positions of the observer, in the solar frame. Each row is a position vector, and each row is in order of observation.
        rho_vec (np.ndarray): List of unit vectors pointing from the observer to the target. Each row is a unit vector, and each row is in order of observation.
        mu (float, Optional): Gravitational parameter of the system in m^3 / s^2

    Returns:
        np.ndarray: State vector containing position (first three elements) and velocity (last three elements.)
    """
    # Step 1: Time scalar factors
    tau1 = t1
    tau3 = t3

    # Observer and LOS unit vectors

    # Calculate a values
    a1 = tau3 / (tau3 - tau1)
    a1u = tau3 * ((tau3 - tau1)**2 - tau3**2) / (6 * (tau3 - tau1))
    a3 = - tau1 / (tau3 - tau1)
    a3u = - tau1 * ((tau3 - tau1)**2 - tau3**2) / (6 * (tau3 - tau1))

    # Invert the LOS, Calculate M Matrix
    L_inv = np.linalg.inv(rho_vec.transpose())
    M = np.matmul(L_inv, obs_pos.transpose())

    # Calculate constants used by the Polynomial
    d1 = M[1][0] * a1 - M[1][1] + M[1][2] * a3
    d2 = M[1][0] * a1u + M[1][2] * a3u
    C = np.dot(rho_vec[1], obs_pos[1])

    # Calculate coefficiants for the polynomial and solve
    p6 = -1 * (d1**2 + 2 * C * d1 + np.dot(obs_pos[1], obs_pos[1]))
    p3 = -2 * mu * (C * d2 + d1*d2)
    p0 = -mu**2 * d2**2
    p = [1,0,p6,0,0,p3,0,0,p0]
    roots:np.ndarray = np.roots(p)
    real_roots = np.array([r2 for r2 in roots if np.isreal(r2)]) # Should be either negative or positive
    u = mu / (np.real(max(real_roots)**3))
    
    # More constants
    c1 = a1 + a1u * u
    c2 = -1
    c3 = a3 + a3u * u

    # Solve linear equation for rho values
    c_rho_vec:np.ndarray = np.matmul(M, np.array([-c1, -c2, -c3]))
    rho1 = c_rho_vec[0] / c1
    rho2 = c_rho_vec[1] / c2
    rho3 = c_rho_vec[2] / c3

    # Calculate Position vectors of the asteroid
    r1 = rho1 * rho_vec[0] + obs_pos[0]
    r2 = rho2 * rho_vec[1] + obs_pos[1]
    r3 = rho3 * rho_vec[2] + obs_pos[2]

    print(f"Position vectors solved!")
    print(f"r1 = {r1}")
    print(f"r2 = {r2}")
    print(f"r3 = {r3}")
    print(f"|r1| = {np.linalg.norm(r1)/EARTH_SMA} AU")
    print(f"|r2| = {np.linalg.norm(r2)/EARTH_SMA} AU")
    print(f"|r3| = {np.linalg.norm(r3)/EARTH_SMA} AU")

