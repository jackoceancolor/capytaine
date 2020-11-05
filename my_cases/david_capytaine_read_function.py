# 2020-05-03: This is a start at a frequency-domain floating support structure model for WEIS-WISDEM.
#             Based heavily on the GA-based optimization framework from Hall 2013

# 2020-05-23: Starting as a script that will evaluate a design based on properties/data provided in separate files.


import numpy as np
import xarray as xr
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import capytaine as cpt

import sys
sys.path.insert(1, '/code/MoorPy')
import MoorPy as mp

# reload the libraries each time in case we make any changes
import importlib
mp = importlib.reload(mp)


## This class represents linear (for now cylinderical) components in the substructure. 
#  It is meant to correspond to Member objects in the WEIS substructure ontology, but containing only 
#  the properties relevant for the Level 1 frequency-domain model, as well as additional data strucutres
#  used during the model's operation.
class Member:

    def __init__(self, strin, nw):
        '''Initialize a Member. For now, this function accepts a space-delimited string with all member properties.
        
        '''
    
        # note: haven't decided yet how to lump masses and handle segments <<<<
    
        entries = strin.split()
        
        self.id = np.int(entries[0])
        self.type = np.int(entries[1])
        self.dA  = np.float(entries[2])                      # diameter of (lower) node
        self.dB  = np.float(entries[3])
        self.rA = np.array(entries[4:7], dtype=np.double)  # x y z of lower node
        self.rB = np.array(entries[7:10], dtype=np.double)
        self.t  = np.float(entries[10])           # shell thickness [m]
        
        self.l_fill = np.float(entries[11])                    # length of member (from end A to B) filled with ballast [m]
        self.rho_fill = np.float(entries[12])               # density of ballast in member [kg/m^3]
        
        rAB = self.rB-self.rA
        self.l = np.linalg.norm(rAB)  # member length
                
        self.q = rAB/self.l           # member axial unit vector
        self.p1 = np.zeros(3)              # member transverse unit vectors (to be filled in later)
        self.p2 = np.zeros(3)              # member transverse unit vectors
        
        
        self.Cd_q  = 0.1  # drag coefficients
        self.Cd_p1 = 0.6
        self.Cd_p2 = 0.6
        self.Ca_q  = 0.0  # added mass coefficients
        self.Ca_p1 = 0.97
        self.Ca_p2 = 0.97
                     
        
        self.n = 10 # number of nodes per member
        self.dl = self.l/self.n   #lumped node lenght (I guess, for now) <<<
        
        self.w = 1    # mass per unit length (kg/m)
        
        self.r  = np.zeros([self.n,3])  # undisplaced node positions along member  [m]
        self.d  = np.zeros( self.n   )  # local diameter along member [m]
        for i in range(self.n):
            self.r[i,:] = self.rA + (i/(self.n-1))*rAB                # spread evenly for now
            self.d[i]   = self.dA + (i/(self.n-1))*(self.dB-self.dA)  # spread evenly since uniform taper
        
        # complex frequency-dependent amplitudes of quantities at each node along member (to be filled in later)
        self.dr = np.zeros([self.n,3,nw], dtype=complex)  # displacement
        self.v  = np.zeros([self.n,3,nw], dtype=complex)  # velocity
        self.a  = np.zeros([self.n,3,nw], dtype=complex)  # acceleration
        self.u  = np.zeros([self.n,3,nw], dtype=complex)  # wave velocity
        self.ud = np.zeros([self.n,3,nw], dtype=complex)  # wave acceleration
        
        
    def getDirections(self):
        '''Returns the direction unit vector for the member based on its end positions'''
        
        q  =  self.q
        p1 = np.cross( np.array([0,1,0]), q) # transverse unit direction vector in X-Z plane
        p2 = np.cross( q, p1 )               # the other transverse unit vector
        
        return q, p1, p2
    
    
    
    def getInertia(self):
        
        # Volume of steel based on the shell thickness [m^3]
        dAi = self.dA - 2*self.t
        dBi = self.dB - 2*self.t
        V_outer = (np.pi/4)*(1/3)*(self.dA**2+self.dB**2+self.dA*self.dB)*self.l
        V_inner = (np.pi/4)*(1/3)*(dAi**2+dBi**2+dAi*dBi)*self.l
        v_steel = V_outer - V_inner         #[m^3] Volume of steel of the member  ((self.t/2)*(self.dA+self.dB)-self.t**2)*np.pi*self.l
        
        # Ballast (future work - this doesn't account for any angle in the member. If the member is horizontal, the ballast doesn't follow gravity)
        dB_fill = (dBi-dAi)*(self.l_fill/self.l) + dAi       # interpolated diameter of member where the ballast is filled to
        v_fill = (np.pi/4)*(1/3)*(dAi**2+dB_fill**2+dAi*dB_fill)*self.l_fill    #[m^3]
        m_fill = self.rho_fill*v_fill                                           #[kg]
        
        # Center of mass
        hco = self.l*((self.dA**2 + 2*self.dA*self.dB + 3*self.dB**2)/(4*(self.dA**2 + self.dA*self.dB + self.dB**2)))  
        hci = self.l*((dAi**2 + 2*dAi*dBi + 3*dBi**2)/(4*(dAi**2 + dAi*dBi + dBi**2)))
        hc_shell = ((hco*V_outer)-(hci*V_inner))/(V_outer-V_inner)  # [m] CoG of member shell in relation to bottom node @ self.rA
        
        hc_fill = self.l_fill*((dAi**2 + 2*dAi*dB_fill + 3*dB_fill**2)/(4*(dAi**2 + dAi*dB_fill + dB_fill**2)))
        
        hc = ((hc_fill*self.rho_fill*v_fill)+(hc_shell*rho_steel*v_steel))/((self.rho_fill*v_fill)+(rho_steel*v_steel))
        
        center = self.rA + (self.q*hc)

        # Moment of Inertia (equations from HydroDyn paper)
        # Calc I@end for outer solid - Calc I@end for "inner" solid - I_outer-I_inner = I_shell @end - PA theorem to calc I@CoG
        r1 = self.dA/2
        r2 = self.dB/2
        m = (r2-r1)/self.l
        r1i = (self.dA/2)-self.t
        r2i = (self.dB/2)-self.t
        mi = (r2i-r1i)/self.l
        if m==0:
            Ir_end_outer = (1/12)*(rho_steel*self.l*np.pi*r1**2)*(3*r1**2 + 4*self.l**2) #[kg-m^2]    about end node
            Ir_end_inner = (1/12)*(rho_steel*self.l*np.pi*r1i**2)*(3*r1i**2 + 4*self.l**2) #[kg-m^2]  about end node
            Ir_end_steel = Ir_end_outer - Ir_end_inner                     # I_outer - I_inner = I_shell -- about end node
            #Ir_end_steel = (1/12)*v_steel*rho_steel*(3*(r1**2 + r1i**2) + 4*self.l**2)
            
            Ir_end_fill = (1/12)*(self.rho_fill*self.l_fill*np.pi*r1i**2)*(2*r1i**2 + 4*self.l_fill**2) #[kg-m^2]  about end node
            
            Ir_end = Ir_end_steel + Ir_end_fill
            
            I_rad = Ir_end - ((rho_steel*v_steel)+m_fill)*hc**2
            
            #I_rad_steel = Ir_end - (rho_steel*v_steel)*hc**2   # about CoG
            
            #I_rad_fill = Ir_end_fill - m_fill*hc**2  # about CoG
            #I_rad = I_rad_steel + I_rad_fill   # sum of all masses about the CoG
            
            I_ax_outer = (1/2)*rho_steel*np.pi*self.l*(r1**4)
            I_ax_inner = (1/2)*rho_steel*np.pi*self.l*(r1i**4)
            I_ax_steel = I_ax_outer - I_ax_inner
            I_ax_fill = (1/2)*self.rho_fill*np.pi*self.l_fill*(r1i**4)
            I_ax = I_ax_steel + I_ax_fill
        else:
            Ir_tip_outer = abs((np.pi/20)*(rho_steel/m)*(1+(4/m**2))*(r2**5-r1**5))                                 # outer, about tip
            Ir_end_outer = abs(Ir_tip_outer - ((rho_steel/(3*m**2))*np.pi*(r2**3-r1**3)*((r1/m)+2*hc)*r1))          # outer, about node
            Ir_tip_inner = abs((np.pi/20)*(rho_steel/mi)*(1+(4/mi**2))*(r2i**5-r1i**5))                             # inner, about tip
            Ir_end_inner = abs(Ir_tip_inner - ((rho_steel/(3*mi**2))*np.pi*(r2i**3-r1i**3)*((r1i/mi)+2*hc)*r1i))    # inner, about node
            Ir_end = Ir_end_outer - Ir_end_inner                                                                    # shell, about node
            I_rad_steel = Ir_end - (rho_steel*v_steel)*hc**2                                                        # shell, about CoG by PAT
            
            I_ax_outer = (rho_steel*np.pi/(10*m))*(r2**5-r1**5)
            I_ax_inner = (rho_steel*np.pi/(10*mi))*(r2i**5-r1i**5)
            I_ax_steel = I_ax_outer - I_ax_inner
            
            if self.l_fill == 0:
                I_rad_fill = 0
                I_ax_fill = 0
            else:
                r2_fill = dB_fill/2
                mi_fill = (r2_fill-r1i)/self.l_fill 
                Ir_tip_fill = abs((np.pi/20)*(rho_steel/mi_fill)*(1+(4/mi_fill**2))*(r2_fill**5-r1i**5))
                Ir_end_fill = abs(Ir_tip_fill - ((self.rho_fill/(3*mi_fill**2))*np.pi*(r2_fill**3-r1i**3)*((r1i/mi_fill)+2*hc)*r1i))    # inner, about node
                I_rad_fill = Ir_end_fill - m_fill*hc**2   # about CoG

                I_ax_fill = (self.rho_fill*np.pi/(10*mi_fill))*(r2_fill**5-r1i**5)
            
            I_rad = I_rad_steel + I_rad_fill # about CoG
            I_ax = I_ax_steel + I_ax_fill 

        return v_steel, center, I_rad, I_ax, m_fill
        
    
    
    def getHydrostatics(self):
        '''Calculates member hydrostatic properties, namely buoyancy and stiffness matrix'''
        
            
        # partially submerged case
        if self.r[0,2]*self.r[-1,2] <= 0:    # if member crosses (or touches) water plane
                
            # angles
            beta = np.arctan2(q[1],q[0])  # member incline heading from x axis
            phi  = np.arctan2(np.sqrt(q[0]**2 + q[1]**2), q[2])  # member incline angle from vertical
            
            # precalculate trig functions
            cosPhi=np.cos(phi)
            sinPhi=np.sin(phi)
            tanPhi=np.tan(phi)
            cosBeta=np.cos(beta)
            sinBeta=np.sin(beta)
            tanBeta=sinBeta/cosBeta
                
            # -------------------- buoyancy and waterplane area properties ------------------------
                
            dWP = np.interp(0, self.r[:,2], self.d)       # diameter of member where its axis crosses the waterplane 
            xWP = np.interp(0, self.r[:,2], self.r[:,0])  # x coordinate where member axis cross the waterplane [m]
            yWP = np.interp(0, self.r[:,2], self.r[:,1])  # y coordinate where member axis cross the waterplane [m]
            AWP = (np.pi/4)*dWP**2                        # waterplane area of member [m^2]
            IWP = (np.pi/64)*dWP**4                       # waterplane moment of inertia [m^4] approximates the waterplane area as the shape of a circle
            
            LWP = abs(self.r[0,2])/cosPhi                 # get length of member that is still underwater. Assumes self.r is about global coords -> z=0 @ SWL
            
            
            # Total enclosed underwater volume
            V_UW = (np.pi/4)*(1/3)*(self.dA**2+dWP**2+self.dA*dWP)*LWP       #[m^3] 
            
            L_center = TaperCV(0.5*self.dA, 0.5*dWP, LWP) # distance from end A to center of buoyancy of member [m]
            
            r_center = self.rA + self.q*L_center          # absolute coordinates of center of volume [m]
        
            
            
        
            # >>>> question: should this function be able to use displaced/rotated values? <<<<
            
            # ------------- get hydrostatic derivatives ---------------- 
            
            # derivatives from global to local 
            dPhi_dThx  = -sinBeta                     # \frac{d\phi}{d\theta_x} = \sin\beta
            dPhi_dThy  =  cosBeta
            #dBeta_dThx = -cosBeta/tanBeta**2
            #dBeta_dThy = -sinBeta/tanBeta**2
            
            # note: below calculations are based on untapered case, but 
            # temporarily approximated for taper by using dWP (diameter at water plane crossing) <<< this is rough
            
            # buoyancy force and moment about end A
            Fz = rho*g* V_UW
            M  = -rho*g*pi*( dWP**2/32*(2.0 + tanPhi**2) + 0.5*(self.rA[2]/cosPhi)**2)*sinPhi  # moment about axis of incline
            Mx = M*dPhi_dThx
            My = M*dPhi_dThy
            
            Fvec = np.zeros(6)                         # buoyancy force (and moment) vector (about PRP)
            Fvec[2] = Fz                               # vertical buoyancy force [N]
            Fvec[3] = Mx + Fz*self.rA[1]                # moment about x axis [N-m]
            Fvec[4] = My - Fz*self.rA[0]                # moment about y axis [N-m]
            
            # derivatives aligned with incline heading
            dFz_dz   = -rho*g*pi*0.25*dWP**2                  /cosPhi
            dFz_dPhi =  rho*g*pi*0.25*dWP**2*self.rA[2]*sinPhi/cosPhi**2
            dM_dz    =  1.0*dFz_dPhi
            dM_dPhi  = -rho*g*pi*0.25*dWP**2 * (dWP**2/32*(2*cosPhi + sinPhi**2/cosPhi + 2*sinPhi/cosPhi**2) + 0.5*self.rA[2]**2*(sinPhi**2+1)/cosPhi**3)
            
            # <<<<<<<<<<<< (end warning) <<<<<<<<<
            
            
            # derivatives in global coordinates about platform reference point
            #dFz_dz is already taken care of
            dFz_dThx = -dFz_dz*self.rA[1] + dFz_dPhi*dPhi_dThx
            dFz_dThy =  dFz_dz*self.rA[0] + dFz_dPhi*dPhi_dThy

            dMx_dz   = -dFz_dz*self.rA[1] + dM_dz*dPhi_dThy  #  = dFz_dThx
            dMy_dz   =  dFz_dz*self.rA[0] + dM_dz*dPhi_dThx  #  = dFz_dThy
            
            dMx_dThx = ( dFz_dz*self.rA[1] + dFz_dPhi*dPhi_dThx)*self.rA[1] + dM_dPhi*dPhi_dThx*dPhi_dThx   
            dMx_dThy = (-dFz_dz*self.rA[0] + dFz_dPhi*dPhi_dThy)*self.rA[1] + dM_dPhi*dPhi_dThx*dPhi_dThy  
            dMy_dThx =-( dFz_dz*self.rA[1] + dFz_dPhi*dPhi_dThy)*self.rA[0] + dM_dPhi*dPhi_dThy*dPhi_dThx  
            dMy_dThy =-(-dFz_dz*self.rA[0] + dFz_dPhi*dPhi_dThy)*self.rA[0] + dM_dPhi*dPhi_dThy*dPhi_dThy  
            
            
            
            # fill in stiffness matrix
            Cmat = np.zeros([6,6]) # hydrostatic stiffness matrix (about PRP)
            '''
            Cmat[2,2] = -dFz_dz
            Cmat[2,3] = -dFz_dThx
            Cmat[2,4] = -dFz_dThy
            Cmat[3,2] = -dMx_dz
            Cmat[4,2] = -dMy_dz   # ignoring symmetries for now, as a way to check equations
            Cmat[3,3] = -dMx_dThx
            Cmat[3,4] = -dMx_dThy
            Cmat[4,3] = -dMy_dThx
            Cmat[4,4] = -dMy_dThy
            '''
            # normal approach to hydrostatic stiffness, using this temporarily until above fancier approach is verified
            Iwp = np.pi*dWP**4/64 # [m^4] Moment of Inertia of the waterplane
            Iwp = np.pi*dWP**4/64 # [m^4] Moment of Inertia of the waterplane
            Cmat[2,2] = -dFz_dz
            Cmat[2,3] = rho*g*(     -AWP*yWP    )
            Cmat[2,4] = rho*g*(      AWP*xWP    )
            Cmat[3,2] = rho*g*(     -AWP*yWP    )
            Cmat[3,3] = rho*g*(IWP + AWP*yWP**2 )
            Cmat[3,4] = rho*g*(      AWP*xWP*yWP)
            Cmat[4,2] = rho*g*(      AWP*xWP    )
            Cmat[4,3] = rho*g*(      AWP*xWP*yWP)
            Cmat[4,4] = rho*g*(IWP + AWP*xWP**2 )
            
            
        
        # fully submerged case 
        elif self.r[0,2] <= 0 and self.r[-1,2] <= 0:
        
            AWP = 0
            IWP = 0
            xWP = 0
            yWP = 0
        
            V_UW  = TaperV( 0.5*self.dA, 0.5*self.dB, self.l)  # displaced volume of member [m^3]
            
            alpha = TaperCV(0.5*self.dA, 0.5*self.dB, 1.0)  # relative location of center of volume from end A (0) to B (1)
            
            r_center = self.rA*(1.0-alpha) + self.rB*alpha  # absolute coordinates of center of volume [m]
        
            # buoyancy force (and moment) vector
            Fvec = translateForce3to6DOF( r_center, np.array([0, 0, rho*g*V_UW]) ) 
  
            # hydrostatic stiffness matrix (about end A)
            Cmat = np.zeros([6,6])  
            Cmat[3,3] = rho*g*V_UW * r_center[2]
            Cmat[4,4] = rho*g*V_UW * r_center[2]
            
        else: # if the members are fully above the surface
            
            AWP = 0
            IWP = 0
            xWP = 0
            yWP = 0
            V_UW = 0
            r_center = np.zeros(3)
            Fvec = np.zeros(6)
            Cmat = np.zeros([6,6])
            
        return Fvec, Cmat, V_UW, r_center, AWP, IWP, xWP, yWP


def TaperV(R1, R2, H):
    '''returns the volume of a cylindrical section, possibly with taper'''
    
    if R1 == R2:             # if just a cylinder
        return np.pi*R1*R1*H
        #taperCV = H/2.0

    elif R1 == 0:             # seperate this case out because it gives a divide by zero in general formula
        return 1./3.*np.pi*R2*R2*H;                                            # cone volume
        #taperCV = 3./4.*H                                                     # from base
    
    else:
        coneH = H/(1.-R2/R1);                                                  # conical height
        coneV = 1./3.*np.pi*R1*R1*coneH;                                       # cone volume
        coneVtip = 1./3.*np.pi*R2*R2*(coneH-H);                                # height from end of taper to cone tip
        return coneV-coneVtip;                                                 # taper volume
        #taperCV = ( coneV*1./4.*coneH - coneVtip*(1./4.*(coneH-H) + H) )/ taperV # from base
    
    return taperV
    
    
def TaperCV(R1, R2, H):
    '''returns the height of the center of buoyancy from the lower node'''
    return H*(R1**2 + 2*R1*R2 + 3*R2**2)*0.25/(R1**2 + R1*R2 + R2**2)
    

def getVelocity(r, Xi, ws):
    '''Get node complex velocity spectrum based on platform motion's and relative position from PRP'''
    
    nw = len(ws)
        
    dr = np.zeros([3,nw], dtype=complex) # node displacement complex amplitudes
    v  = np.zeros([3,nw], dtype=complex) # velocity
    a  = np.zeros([3,nw], dtype=complex) # acceleration
        
    
    for i in range(nw):
        dr[:,i] = Xi[:3,i] + SmallRotate(r, Xi[3:,i])
        v[ :,i] = 1j*ws[i]*dr[:,i]
        a[ :,i] = 1j*ws[i]*v[ :,i]
    

    return dr, v, a # node dispalcement, velocity, acceleration (Each  [3 x nw])
        
    
## Get wave velocity and acceleration complex amplitudes based on wave spectrum at a given location
def getWaveKin(zeta0, w, k, h, r):

    # inputs: wave elevation fft, wave freqs, wave numbers, depth, point position

    beta = 0  # no wave heading for now

    zeta = np.zeros(nw , dtype=complex ) # local wave elevation
    u  = np.zeros([3,nw], dtype=complex) # local wave kinematics velocity
    ud = np.zeros([3,nw], dtype=complex) # local wave kinematics acceleration
    
    for i in range(nw):
    
        # .............................. wave elevation ...............................
        zeta[i] = zeta0[i]* np.exp( -1j*(k[i]*(np.cos(beta)*r[0] + np.sin(beta)*r[1])))  # shift each zetaC to account for location
        
        
        # ...................... wave velocities and accelerations ............................
        
        z = r[2]
        
        # only process wave kinematics if this node is submerged
        if z < 0:
            
            # Calculate SINH( k*( z + h ) )/SINH( k*h ) and COSH( k*( z + h ) )/SINH( k*h )
            # given the wave number, k, water depth, h, and elevation z, as inputs.
            if (    k[i]   == 0.0  ):                   # When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.
                SINHNumOvrSIHNDen = 1.0
                COSHNumOvrSIHNDen = 99999
            elif ( k[i]*h >  89.4 ):                # When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, with h > 14.23*wavelength (since k = 2*Pi/wavelength) we can use the numerically-stable deep water formulation instead.
                SINHNumOvrSIHNDen = np.exp(  k[i]*z );
                COSHNumOvrSIHNDen = np.exp(  k[i]*z );
            else:                                    # 0 < k*h <= 89.4; use the shallow water formulation.
                SINHNumOvrSIHNDen = np.sinh( k[i]*( z + h ) )/np.sinh( k[i]*h );
                COSHNumOvrSIHNDen = np.cosh( k[i]*( z + h ) )/np.sinh( k[i]*h );
            
            # Fourier transform of wave velocities 
            u[0,i] =    w[i]* zeta[i]*COSHNumOvrSIHNDen *np.cos(beta) 
            u[1,i] =    w[i]* zeta[i]*COSHNumOvrSIHNDen *np.sin(beta)
            u[2,i] = 1j*w[i]* zeta[i]*SINHNumOvrSIHNDen

            # Fourier transform of wave accelerations                   
            ud[:,i] = 1j*w[i]*u[:,i]
        
    return u, ud


# calculate wave number based on wave frequency in rad/s and depth
def waveNumber(omega, h, e=0.001):
    
    g = 9.81
    # omega - angular frequency of waves
    
    k = 0                           #initialize  k outside of loop
                                    #error tolerance
    k1 = omega*omega/g                  #deep water approx of wave number
    k2 = omega*omega/(np.tanh(k1*h)*g)      #general formula for wave number
    while np.abs(k2 - k1)/k1 > e:           #repeate until converges
        k1 = k2
        k2 = omega*omega/(np.tanh(k1*h)*g)
    
    k = k2
    
    return k
    

# translate point at location r based on three small angles in th
def SmallRotate(r, th):
    
    rt = np.zeros(3, dtype=complex)    # translated point
    
    rt[0] =              th[2]*r[1] - th[1]*r[2]
    rt[0] = th[2]*r[0]              - th[0]*r[2]
    rt[0] = th[1]*r[0] - th[0]*r[1]
    
    return rt
    
    
# given a size-3 vector, vec, return the matrix from the multiplication vec * vec.transpose
def VecVecTrans(vec):

    vvt = np.zeros([3,3])
    
    for i in range(3):
        for j in range(3):
            vvt[i,j]= vec[i]*vec[j]
            
    return vvt
    

# produce alternator matrix
def getH(r):

    H = np.zeros([3,3])
    H[0,1] =  r[2];
    H[1,0] = -r[2];
    H[0,2] = -r[1];
    H[2,0] =  r[1];
    H[1,2] =  r[0];
    H[2,1] = -r[0];
    
    return H



def translateForce3to6DOF(r, Fin):
    '''Takes in a position vector and a force vector (applied at the positon), and calculates 
    the resulting 6-DOF force and moment vector.    
    
    :param array r: x,y,z coordinates at which force is acting [m]
    :param array Fin: x,y,z components of force [N]
    :return: the resulting force and moment vector
    :rtype: array
    '''
    Fout = np.zeros(6, dtype=Fin.dtype) # initialize output vector as same dtype as input vector (to support both real and complex inputs)
    
    Fout[:3] = Fin
    
    Fout[3:] = np.cross(r, Fin)
    
    return Fout

    
    
# translate mass matrix to make 6DOF mass-inertia matrix
def translateMatrix3to6DOF(r, Min):
    '''Transforms a 3x3 matrix to be about a translated reference point, resulting in a 6x6 matrix.'''
    
    # sub-matrix definitions are accordint to  | m    J |
    #                                          | J^T  I |

    # note that the J term and I terms are zero in this case because the input is just a mass matrix (assumed to be about CG)

    H = getH(r)     # "anti-symmetric tensor components" from Sadeghi and Incecik
    
    Mout = np.zeros([6,6]) #, dtype=complex)
    
    # mass matrix  [m'] = [m]
    Mout[:3,:3] = Min
        
    # product of inertia matrix  [J'] = [m][H] + [J]  
    Mout[:3,3:] = np.matmul(Min, H)
    Mout[3:,:3] = Mout[:3,3:].T
    
    # moment of inertia matrix  [I'] = [H][m][H]^T + [J]^T [H] + [H]^T [J] + [I]
    
    Mout[3:,3:] = np.matmul(np.matmul(H,Min), H.T) 
    
    return Mout


def translateMatrix6to6DOF(r, Min):
    '''Transforms a 6x6 matrix to be about a translated reference point.'''
    
    # sub-matrix definitions are accordint to  | m    J |
    #                                          | J^T  I |

    H = getH(r)     # "anti-symmetric tensor components" from Sadeghi and Incecik
    
    Mout = np.zeros([6,6]) #, dtype=complex)
    
    # mass matrix  [m'] = [m]
    Mout[:3,:3] = Min[:3,:3]
        
    # product of inertia matrix  [J'] = [m][H] + [J]
    Mout[:3,3:] = np.matmul(Min[:3,:3], H) + Min[:3,3:]
    Mout[3:,:3] = Mout[:3,3:].T
    
    # moment of inertia matrix  [I'] = [H][m][H]^T + [J]^T [H] + [H]^T [J] + [I]
    Mout[3:,3:] = np.matmul(np.matmul(H,Min[:3,:3]), H.T) + np.matmul(Min[3:,:3], H) + np.matmul(H.T, Min[:3,3:]) + Min[3:,3:]
        
    return Mout
    

def JONSWAP(ws, Hs, Tp, Gamma=1.0):
    '''Returns the JONSWAP wave spectrum for the given frequencies and parameters.
    
    Parameters
    ----------    
    ws : float | array
        wave frequencies to compute spectrum at (scalar or 1-D list/array) [rad/s]
    Hs : float
        significant wave height of spectrum [m]
    Tp : float
        peak spectral period [s]
    Gamma : float
        wave peak shape parameter []. The default value of 1.0 gives the Pierson-Moskowitz spectrum.
        
    Returns
    -------
    S : array
        wave power spectral density array corresponding to frequencies in ws [m^2/Hz]
        
    
    This function calculates and returns the one-sided power spectral density spectrum 
    at the frequency/frequencies (ws) based on the provided significant wave height,
    peak period, and (optionally) shape parameter gamma.
    
    This formula for the JONSWAP spectrum is adapted from FAST v7 and based
    on what's documented in IEC 61400-3.
    '''
    
    # handle both scalar and array inputs
    if isinstance(ws, (list, tuple, np.ndarray)):
        ws = np.array(ws)
    else:
        ws = np.array([ws])

    # initialize output array
    S = np.zeros(len(ws))    
 
            
    # the calculations
    f        = 0.5/np.pi * ws                         # wave frequencies in Hz
    fpOvrf4  = pow((Tp*f), -4.0)                      # a common term, (fp/f)^4 = (Tp*f)^(-4)
    C        = 1.0 - ( 0.287*np.log(Gamma) )          # normalizing factor
    Sigma = 0.07*(f <= 1.0/Tp) + 0.09*(f > 1.0/Tp)    # scaling factor
    
    Alpha = np.exp( -0.5*((f*Tp - 1.0)/Sigma)**2 )

    return  0.5/np.pi *C* 0.3125*Hs*Hs*fpOvrf4/f *np.exp( -1.25*fpOvrf4 )* Gamma**Alpha
    
    
def printMat(mat):
    '''Print a matrix'''
    for i in range(mat.shape[0]):
        print( "\t".join(["{:+8.3e}"]*mat.shape[1]).format( *mat[i,:] ))
        
def printVec(vec):
    '''Print a vector'''
    print( "\t".join(["{:+8.3e}"]*len(vec)).format( *vec ))
    


# ------------------------------- basic setup -----------------------------------------

nDOF = 6

w = np.arange(.01, 3, 0.01)  # angular frequencies tp analyze (rad/s)
nw = len(w)  # number of frequencies

k= np.zeros(nw)  # wave number


# ----------------------- member-based platform description --------------------------

# (hard-coded for now - set to OC3 Hywind Spar geometry - eventually these will be provided as inputs instead)

# list of member objects
memberList = []

# -------------------- OC3 Hywind Spar ----------------------
'''
# ------------------ turbine Tower description ------------------
# diameters and thicknesses linearly interpolated from dA[0] to dB[-1] and t[0] to t[-1]
#                      number   type    dA      dB      xa      ya     za     xb     yb      zb      t     l_fill  rho_ballast

memberList.append(Member(" 1     1    6.500   6.237    0.0    0.0    10.00   0.0    0.0    17.76   0.0270   0.0    1025.0  ", nw))
memberList.append(Member(" 2     1    6.237   5.974    0.0    0.0    17.76   0.0    0.0    25.52   0.0262   0.0    1025.0  ", nw))
memberList.append(Member(" 3     1    5.974   5.711    0.0    0.0    25.52   0.0    0.0    33.28   0.0254   0.0    1025.0  ", nw))
memberList.append(Member(" 4     1    5.711   5.448    0.0    0.0    33.28   0.0    0.0    41.04   0.0246   0.0    1025.0  ", nw))
memberList.append(Member(" 5     1    5.448   5.185    0.0    0.0    41.04   0.0    0.0    48.80   0.0238   0.0    1025.0  ", nw))

memberList.append(Member(" 6     1    5.185   4.922    0.0    0.0    48.80   0.0    0.0    56.56   0.0230   0.0    1025.0  ", nw))
memberList.append(Member(" 7     1    4.922   4.659    0.0    0.0    56.56   0.0    0.0    64.32   0.0222   0.0    1025.0  ", nw))
memberList.append(Member(" 8     1    4.659   4.396    0.0    0.0    64.32   0.0    0.0    72.08   0.0214   0.0    1025.0  ", nw))
memberList.append(Member(" 9     1    4.396   4.133    0.0    0.0    72.08   0.0    0.0    79.84   0.0206   0.0    1025.0  ", nw))
memberList.append(Member("10     1    4.133   3.870    0.0    0.0    79.84   0.0    0.0    87.60   0.0198   0.0    1025.0  ", nw))

# ---------- spar platform substructure description --------------
memberList.append(Member("11     2    9.400   9.400    0.0    0.0    -120.   0.0    0.0    -12.0   0.0270   52.    1850.0  ", nw))
#memberList.append(Member("11     2    9.400   9.400    0.0    0.0    -120.   0.0    0.0    -12.0   0.066   41.4    2000.0  ", nw))
memberList.append(Member("12     2    9.400   6.500    0.0    0.0    -12.0   0.0    0.0    -4.00   0.0270   0.0    1025.0  ", nw))
memberList.append(Member("13     2    6.500   6.500    0.0    0.0    -4.00   0.0    0.0    10.00   0.0270   0.0    1025.0  ", nw))

#memberList.append(Member("1      2    6.500   6.500    0.0    0.0    -100.00   0.0    0.0    0.00   0.0270   0.0    1025.0  ", nw))

# -------------------------- turbine RNA description ------------------------
# below are rough properties for NREL 5 MW reference turbine
mRNA    = 110000              + 240000  # RNA mass [kg]
IxRNA   = 11776047*(1 + 1 + 1) + 115926   # RNA moment of inertia about local x axis (assumed to be identical to rotor axis for now, as approx) [kg-m^2]
IrRNA   = 11776047*(1 +.5 +.5) + 2607890   # RNA moment of inertia about local y or z axes [kg-m^2]
xCG_RNA = 0                             # x location of RNA center of mass [m] (Close enough to -0.27 m)
hHub    = 90.0                          # hub height above water line [m]

Fthrust = 800e3  # peak thrust force, [N]
'''



# ------------------- DTU 10 MW Spar --------------------------

# ------------------ turbine Tower description ------------------
#                      number   type    dA      dB      xa      ya     za     xb     yb      zb      t     l_fill  rho_ballast

# new version from 11-114.7 (for OpenFAST reasons)
memberList.append(Member(" 1     1    8.00    7.75    0.0    0.0    11.00    0.0    0.0    21.37    0.038   0.0    1025.0  ", nw))
memberList.append(Member(" 2     1    7.75    7.50    0.0    0.0    21.37    0.0    0.0    31.74    0.036   0.0    1025.0  ", nw))
memberList.append(Member(" 3     1    7.50    7.25    0.0    0.0    31.74    0.0    0.0    42.11    0.034   0.0    1025.0  ", nw))
memberList.append(Member(" 4     1    7.25    7.00    0.0    0.0    42.11    0.0    0.0    52.48    0.032   0.0    1025.0  ", nw))
memberList.append(Member(" 5     1    7.00    6.75    0.0    0.0    52.48    0.0    0.0    62.85    0.030   0.0    1025.0  ", nw))

memberList.append(Member(" 6     1    6.75    6.50    0.0    0.0    62.85    0.0    0.0    73.22    0.028   0.0    1025.0  ", nw))
memberList.append(Member(" 7     1    6.50    6.25    0.0    0.0    73.22    0.0    0.0    83.59    0.026   0.0    1025.0  ", nw))
memberList.append(Member(" 8     1    6.25    6.00    0.0    0.0    83.59    0.0    0.0    93.96    0.024   0.0    1025.0  ", nw))
memberList.append(Member(" 9     1    6.00    5.75    0.0    0.0    93.96    0.0    0.0   104.33    0.022   0.0    1025.0  ", nw))
memberList.append(Member("10     1    5.75    5.50    0.0    0.0   104.33    0.0    0.0   114.70    0.020   0.0    1025.0  ", nw))

# =============================================================================
# # old version from 13-116.63
# memberList.append(Member(" 1     1    8.00    7.75    0.0    0.0    13.000    0.0    0.0    23.363   0.038   0.0    1025.0  ", nw))
# memberList.append(Member(" 2     1    7.75    7.50    0.0    0.0    23.363    0.0    0.0    33.726   0.036   0.0    1025.0  ", nw))
# memberList.append(Member(" 3     1    7.50    7.25    0.0    0.0    33.726    0.0    0.0    44.089   0.034   0.0    1025.0  ", nw))
# memberList.append(Member(" 4     1    7.25    7.00    0.0    0.0    44.089    0.0    0.0    54.452   0.032   0.0    1025.0  ", nw))
# memberList.append(Member(" 5     1    7.00    6.75    0.0    0.0    54.452    0.0    0.0    64.815   0.030   0.0    1025.0  ", nw))
# 
# memberList.append(Member(" 6     1    6.75    6.50    0.0    0.0    64.815    0.0    0.0    75.178   0.028   0.0    1025.0  ", nw))
# memberList.append(Member(" 7     1    6.50    6.25    0.0    0.0    75.178    0.0    0.0    85.541   0.026   0.0    1025.0  ", nw))
# memberList.append(Member(" 8     1    6.25    6.00    0.0    0.0    85.541    0.0    0.0    95.904   0.024   0.0    1025.0  ", nw))
# memberList.append(Member(" 9     1    6.00    5.75    0.0    0.0    95.904    0.0    0.0   106.267   0.022   0.0    1025.0  ", nw))
# memberList.append(Member("10     1    5.75    5.50    0.0    0.0   106.267    0.0    0.0   116.630   0.020   0.0    1025.0  ", nw))
# =============================================================================



# ---------- spar platform substructure description --------------

# =============================================================================
# Ballast members from Senu's sizing
# memberList.append(Member("11     2    14.75   14.75    0.0    0.0    -90.   0.0    0.0    -85.2   0.046   4.8    3743.42  ", nw))
# memberList.append(Member("12     2    14.75   14.75    0.0    0.0    -85.2   0.0    0.0    -75.708   0.046   9.492    3792.35  ", nw))
# memberList.append(Member("13     2    14.75   14.75    0.0    0.0    -75.708   0.0    0.0    -72.734   0.046   2.974    1883.78  ", nw))
# =============================================================================

# Ballast members from Stein getting weight = displ
memberList.append(Member("11     2    14.75   14.75    0.0    0.0    -90.000   0.0    0.0    -85.200   0.046    4.800    7850.  ", nw))
memberList.append(Member("12     2    14.75   14.75    0.0    0.0    -85.200   0.0    0.0    -75.708   0.046    9.492    2650.  ", nw))
memberList.append(Member("13     2    14.75   14.75    0.0    0.0    -75.708   0.0    0.0    -72.734   0.046    2.974    1025.  ", nw))

memberList.append(Member("14     2    14.75   14.75    0.0    0.0    -72.734   0.0    0.0    -20.   0.046   0.0    1025.  ", nw))
memberList.append(Member("15     2    14.75    8.00    0.0    0.0    -20.   0.0    0.0    -5.   0.063   0.0    1025.0  ", nw))
memberList.append(Member("16     2     8.00    8.00    0.0    0.0    -5.   0.0    0.0    7.   0.068   0.0    1025.0  ", nw))
memberList.append(Member("17     2     8.00    7.00    0.0    0.0    7.   0.0    0.0    11.   0.055   0.0    1025.0  ", nw))




# -------------------------- turbine RNA description ------------------------
mRotor = 227962 #[kg]
mNacelle = 446036 #[kg]
IxHub = 325671 #[kg-m^2]
IzNacelle = 7326346 #[kg-m^2]
IxBlades = 45671252 #[kg-m^2] MOI value from FAST file, don't know where MOI is about. Assuming about the hub
xCG_Hub = -7.07 #[m] from yaw axis
xCG_Nacelle = 2.687 #[m] from yaw axis


mRNA = mRotor + mNacelle #[kg]
IxRNA = IxBlades*(1 + 1 + 1) + IxHub # RNA moment of inertia about local x axis (assumed to be identical to rotor axis for now, as approx) [kg-m^2]
IrRNA = IxBlades*(1 + .5 + .5) + IzNacelle # RNA moment of inertia about local y or z axes [kg-m^2]


xCG_RNA = ((mRotor*xCG_Hub)+(mNacelle*xCG_Nacelle))/(mRotor+mNacelle)          # x location of RNA center of mass [m]

#hHub    = 119.0                          # hub height above water line [m]
hHub    = 118.0 



# ------- Wind conditions
Fthrust = 800e3  # peak thrust force, [N]
Mthrust = hHub*Fthrust  # overturning moment from turbine thrust force [N-m]





# ---------------- (future work) import hydrodynamic coefficient files ----------------



# ---------------------------- environmental conditions -------------------------------

depth = 200
rho = 1025
g = 9.81

pi = np.pi

# environmental condition(s)
Hs = [8 ];
Tp = [12];
windspeed = [8 ];

S   = np.zeros([len(Hs), nw])   # wave spectrum
S2  = np.zeros([len(Hs), nw])   # wave spectrum
zeta= np.zeros([len(Hs), nw])  # wave elevation

T = 2*np.pi/w # periods



for imeto in range(len(Hs)):       # loop through each environmental condition

    # make wave spectrum (setting to 1 gives approximate RAOs)
    S[imeto,:]  =  JONSWAP(w, Hs[imeto], Tp[imeto])
    
    # wave elevation amplitudes (these are easiest to use) - no need to be complex given frequency domain use
    zeta[imeto,:] = np.sqrt(S[imeto,:])
    

# get wave number
for i in range(nw):
    k[i] = waveNumber(w[i], depth)

'''
plt.plot(w/np.pi/2,  S[2,:], "r")
plt.plot(w/np.pi/2, S2[2,:], "b")
plt.xlabel("Hz")
'''

# ignoring multiple DLCs for now <<<<<
#for imeto = 1:length(windspeeds)
imeto = 0



# ---------------------- set up key arrays -----------------------------

# structure-related arrays
M_struc = np.zeros([6,6])                # structure/static mass/inertia matrix [kg, kg-m, kg-m^2]
B_struc = np.zeros([6,6])                # structure damping matrix [N-s/m, N-s, N-s-m] (may not be used)
C_struc = np.zeros([6,6])                # structure effective stiffness matrix [N/m, N, N-m]
W_struc = np.zeros([6])                  # static weight vector [N, N-m]

# hydrodynamics-related arrays
A_hydro = np.zeros([6,6,nw])             # hydrodynamic added mass matrix [kg, kg-m, kg-m^2]
B_hydro = np.zeros([6,6,nw])             # hydrodynamic added mass matrix [kg, kg-m, kg-m^2]
C_hydro = np.zeros([6,6])                # hydrostatic stiffness matrix [N/m, N, N-m]
W_hydro = np.zeros(6)                    # buoyancy force/moment vector [N, N-m]
F_hydro = np.zeros([6,nw],dtype=complex) # excitation force/moment complex amplitudes vector [N, N-m]

# moorings-related arrays
A_moor = np.zeros([6,6,nw])              # mooring added mass matrix [kg, kg-m, kg-m^2] (may not be used)
B_moor = np.zeros([6,6,nw])              # mooring damping matrix [N-s/m, N-s, N-s-m] (may not be used)
C_moor = np.zeros([6,6])                 # mooring stiffness matrix (linearized about platform offset) [N/m, N, N-m]
W_moor = np.zeros(6)                     # mean net mooring force/moment vector [N, N-m]

# final coefficient arrays
M_tot = np.zeros([6,6,nw])               # total mass and added mass matrix [kg, kg-m, kg-m^2]
B_tot = np.zeros([6,6,nw])               # total damping matrix [N-s/m, N-s, N-s-m]
C_tot = np.zeros([6,6,nw])               # total stiffness matrix [N/m, N, N-m]
F_tot = np.zeros([6,nw], dtype=complex)  # total excitation force/moment complex amplitudes vector [N, N-m]



# --------------- add in linear hydrodynamic coefficients here if applicable --------------------

def interp_hydro_data(wCapy, wDes, addedMass, damping, fEx):
    addedMassInterp = spi.interp1d(wCapy, addedMass, axis=2)(wDes)
    dampingInterp = spi.interp1d(wCapy, damping, axis=2)(wDes)
    fExcitationInterp = spi.interp1d(wCapy, fEx)(wDes)
    return addedMassInterp, dampingInterp, fExcitationInterp


def read_capy_nc(capyFName, wDes=None):
    '''
    Read Capytaine .nc file and return hydrodynamic coefficients in suitable
    format for FrequencyDomain.py

    Parameters
    ----------
    capyFName : str
        String containing path to desired capytaine .nc file
    wDes: array
        Array of desired frequency points - can be different resolution to that
        used by Capytaine (but must be within the range computed by Capytaine)

    Returns
    -------
    frequency range & hydrodynamic coefficients. If wDes passed as argument,
    hydrodynamic coefficients are interpolated for these points (and frequency
    range is not returned)

    Raises
    ------
    ValueError
        when max/min wDes frequency range is outside wCapy range computed by
        Capytaine
    '''

    # read the Capytaine NetCDF data into xarray
    capyData = xr.open_dataset(capyFName)
    wCapy = capyData['omega'].values

    # check desired frequency range is within computed range (for interpolation)
    if wDes is not None:
        if wDes[0]<wCapy[0]:
            raise ValueError(f'\nMinimum desired frequency = {wDes[0]:.2f} rad/s \n'
                             f'Range of computed frequencies = {wCapy[0]:.2f} - {wCapy[-1]:.2f} rad/s \n'
                             f'[out of range]')
        if wDes[-1]>wCapy[-1]:
            raise ValueError(f'\nMaximum desired frequency = {wDes[-1]:.2f} rad/s \n'
                             f'Range of computed frequencies = {wCapy[0]:.2f} - {wCapy[-1]:.2f} rad/s \n'
                             f'[out of range]')

    # convert hydrodynamic coefficients to FrequencyDomain.py format
    addedMass = capyData['added_mass'].values.transpose(1,2,0)
    damping = capyData['radiation_damping'].values.transpose(1,2,0)
    fEx = (np.squeeze(capyData['diffraction_force'].values[0,:].T)
           + 1j*np.squeeze(capyData['diffraction_force'].values[1,:].T))

    # (optional) interpolate coefficients; return appropriate hydrodynamic
    # coefficients 
    if wDes is not None:
        addedMassIntp, dampingIntp, fExIntp = interp_hydro_data(wCapy,
                                                                wDes,
                                                                addedMass,
                                                                damping,
                                                                fEx)
        return wDes, addedMassIntp, dampingIntp, fExIntp
    else:
        return wCapy, addedMass, damping, fEx


def call_capy(meshFName, wCapy, CoG=[0,0,0], headings=[0.0], saveNc=False,
              ncFName=None, wDes=None):
    '''
    call Capytaine for a given mesh, frequency range and wave headings

    Parameters
    ----------
    meshFName : str
        string containing path to hydrodynamic mesh.
        mesh must be cropped at waterline (OXY plane) and have no lid
    wCapy: array
        array of frequency points to be computed by Capytaine
    CoG: list
        3x1 vector of body's CoG
    headings: list
        list of wave headings to compute
    saveNc: Bool
        save results to .nc file
    ncFName: str
        name of .nc file
    wDes: array
        array of desired frequency points
        (for interpolation of wCapy-based Capytaine data)

    Returns
    -------
    hydrodynamic coefficients; as computed or interpolated

    Notes
    -----
    TODO:
    - expand to multibody problems
    '''

    # create capytaine body object
    body = cpt.FloatingBody.from_file(meshFName)
    body.add_all_rigid_body_dofs()
    body.center_of_mass = CoG

    # define the hydrodynamic problems
    problems = [cpt.RadiationProblem(body=body,
                                     radiating_dof=dof,
                                     omega=w,
                                     sea_bottom=-np.infty,
                                     g=9.81,
                                     rho=1025)
                                     for dof in body.dofs for w in wCapy]

    problems += [cpt.DiffractionProblem(body=body,
                                        omega=w,
                                        wave_direction=heading,
                                        sea_bottom=-np.infty,
                                        g=9.81,
                                        rho=1025)
                                        for w in wCapy for heading in headings]

    # call Capytaine solver
    print(f'\n-------------------------------\n'
          f'Calling Capytaine BEM solver...\n'
          f'-------------------------------\n'
          f'mesh = {meshFName}\n'
          f'w range = {wCapy[0]:.3f} - {wCapy[-1]:.3f} rad/s\n'
          f'dw = {(wCapy[1]-wCapy[0]):.3f} rad/s\n'
          f'no of headings = {len(headings)}\n'
          f'no of radiation & diffraction problems = {len(problems)}\n'
          f'-------------------------------\n')

    solver = cpt.BEMSolver()
    results = [solver.solve(problem) for problem in sorted(problems)]
    capyData = cpt.assemble_dataset(results)

    # convert hydrodynamic coefficients to FrequencyDomain.py format
    addedMass = capyData['added_mass'].values.transpose(1,2,0)
    damping = capyData['radiation_damping'].values.transpose(1,2,0)
    fEx = np.squeeze(capyData['diffraction_force']).values.T

    # (optional) save to .nc file and interpolate coefficients; return
    # appropriate hydrodynamic coefficients
    if saveNc == True:
        cpt.io.xarray.separate_complex_values(capyData).to_netcdf(ncFName)

    if wDes is not None:
        addedMassIntp, dampingIntp, fExIntp = interp_hydro_data(wCapy,
                                                                wDes,
                                                                addedMass,
                                                                damping,
                                                                fEx)
        return wDes, addedMassIntp, dampingIntp, fExIntp
    else:
        return wCapy, addedMass, damping, fEx


# --------------- Get general geometry properties including hydrostatics ------------------------

# initialize some variables for running totals
VTOT = 0.                   # Total underwater volume of all members combined
mTOT = 0.                   # Total mass of all members [kg]
AWP_TOT = 0.                # Total waterplane area of all members [m^2]
IWPx_TOT = 0                # Total waterplane moment of inertia of all members about x axis [m^4]  
IWPy_TOT = 0                # Total waterplane moment of inertia of all members about y axis [m^4]  
Sum_V_rCB = np.zeros(3)     # product of each member's buoyancy multiplied by center of buoyancy [m^4]
Sum_AWP_rWP = np.zeros(2)   # product of each member's waterplane area multiplied by the area's center point [m^3]
Sum_M_center = np.zeros(3)  # product of each member's mass multiplied by its center of mass [kg-m] (Only considers the shell mass right now)

# loop through each member
for mem in memberList:
    
    q, p1, p2 = mem.getDirections()                # get unit direction vectors
    
    
    # ---------------------- get member's mass and inertia properties ------------------------------
    rho_steel = 8500 #[kg/m^3]

    v_steel, center, I_rad, I_ax, m_fill = mem.getInertia() # calls the getInertia method to calcaulte values
    
    mass = v_steel*rho_steel + m_fill #[kg]
    Mmat = np.diag([mass, mass, mass, I_rad, I_rad, I_ax]) # MOI matrix = Mmat[3:,3:] is 0 on off diags bc symmetry in cylinders
    # @mhall: Mmat as written above is the mass and inertia matrix about the member CG...@shousner: you betcha
  
    # now convert everything to be about PRP (platform reference point) and add to global vectors/matrices
    W_struc += translateForce3to6DOF( center, np.array([0,0, -g*mass]) )  # weight vector
    M_struc += translateMatrix6to6DOF(center, Mmat)                       # mass/inertia matrix
    # @mhall: Using the diagonal Mmat, and calling the above function with the "center" coordinate, will give the mass/inertia about the PRP!
    # @shousner: center is the position vector of the CG of the member, from the global coordinates aka PRP
    Sum_M_center += center*mass
    
    # -------------------- get each member's buoyancy/hydrostatic properties -----------------------
    
    Fvec, Cmat, V_UW, r_CB, AWP, IWP, xWP, yWP = mem.getHydrostatics()  # call to Member method for hydrostatic calculations
    
    # now convert everything to be about PRP (platform reference point) and add to global vectors/matrices <<<<< needs updating (already about PRP)
    W_hydro += Fvec # translateForce3to6DOF( mem.rA, np.array([0,0, Fz]) )  # weight vector
    C_hydro += Cmat # translateMatrix6to6DOF(mem.rA, Cmat)                       # hydrostatic stiffness matrix
    VTOT    += V_UW    # add to total underwater volume of all members combined
    AWP_TOT += AWP
    IWPx_TOT += IWP + AWP*yWP**2
    IWPy_TOT += IWP + AWP*xWP**2
    Sum_V_rCB   += r_CB*V_UW
    Sum_AWP_rWP += np.array([xWP, yWP])*AWP


# ------------------------- include RNA properties -----------------------------

# for now, turbine RNA is specified by some simple lumped properties
Mmat = np.diag([mRNA, mRNA, mRNA, IxRNA, IrRNA, IrRNA])            # create mass/inertia matrix
center = np.array([xCG_RNA, 0, hHub])                                 # RNA center of mass location

# now convert everything to be about PRP (platform reference point) and add to global vectors/matrices
W_struc += translateForce3to6DOF( center, np.array([0,0, -g*mRNA]) )  # weight vector
M_struc += translateMatrix6to6DOF(center, Mmat)                       # mass/inertia matrix
Sum_M_center += center*mRNA


# ----------- process key hydrostatic-related totals for use in static equilibrium solution ------------------

mTOT = M_struc[0,0]
rCG_TOT = Sum_M_center/mTOT 

rCB_TOT = Sum_V_rCB/VTOT       # location of center of buoyancy on platform

if VTOT==0: # if you're only working with members above the platform, like modeling the wind turbine
    zMeta = 0
else:
    zMeta   = rCB_TOT[2] + IWPx_TOT/VTOT  # add center of buoyancy and BM=I/v to get z elevation of metecenter [m] (have to pick one direction for IWP)


C_struc[3,3] = -mTOT*g*rCG_TOT[2]
C_struc[4,4] = -mTOT*g*rCG_TOT[2]
      


# --------------- set up quasi-static mooring system and solve for mean offsets -------------------


import MoorDesign as md

# =============================================================================
# # Inputs for OC3 Hywind
# depth = 320.
# type_string = 'main'
# anchorR = 853.87
# fairR = 5.2
# fair_depth = 70.
# LineLength = 902.2
# LineD = 0.09
# dryMass_L = 77.7066 #[kg/m]
# EA = 384243000 #[N]
# angle = np.array([0, 2*np.pi/3, -2*np.pi/3])
# =============================================================================

# Inputs for DTU 10 MW setup
depth = 600. #[m]
type_string = 'main'
LineD = 0.15 #[m]
wetMass_L = 4.401 #[kg/m]
dryMass_L = wetMass_L + (np.pi/4)*LineD**2*rho
EA = 384243000 #100e6 #[N]
angle = np.array([0, 2*np.pi/3, -2*np.pi/3])
anchorR = 656.139 #[m]
fair_depth = 21. #[m]
fairR = 7.875 #[m]
LineLength = 868.5 #[m]




MooringSystem = md.make3LineSystem(depth, type_string, LineD, dryMass_L, EA, angle, anchorR, fair_depth, fairR, LineLength)

MooringSystem.BodyList[0].m = mTOT
MooringSystem.BodyList[0].v = VTOT
MooringSystem.BodyList[0].rCG = rCG_TOT
MooringSystem.BodyList[0].AWP = AWP_TOT
MooringSystem.BodyList[0].rM = np.array([0,0,zMeta])
MooringSystem.BodyList[0].f6Ext = np.array([Fthrust,0,0, 0,Mthrust,0])

MooringSystem.initialize()



# --------------- Bridle Confirguration --------------------------
# Inputs for DTU 10 MW setup
depth = 600. #[m]
type_string = ['chain','synth']
LineD = [0.17971, 0.15]  #[m]
#wetMass_L = 4.401 #[kg/m]
#dryMass_L = wetMass_L + (np.pi/4)*LineD**2*rho
dryMass_L = [200, 40.18]
EA = [0, 121415000] #100e6 #[N]
EA[0] = 0.854e11*(LineD[0]**2)
angle = np.array([0, 2*np.pi/3, -2*np.pi/3])
anchorR = 656.139 #[m]
fair_depth = 21. #[m]
fairR = 7.875 #[m]

LineLength = 868.5 #[m]

chainLength = 80
bridleLength = 30
synthLength = LineLength-chainLength-bridleLength

synthR = anchorR-0.9*chainLength
synthZ = 0.95*depth
bridleR = ((synthR-fairR)/(synthLength+bridleLength))*bridleLength + fairR
bridleZ = ((synthZ-fair_depth)/(synthLength+bridleLength))*bridleLength + fair_depth




Bridle = md.makeBridleSystem(depth, type_string, LineD, dryMass_L, EA, angle, anchorR, fair_depth, fairR, LineLength, chainLength, bridleLength, synthLength, synthR, synthZ, bridleR, bridleZ)

Bridle.BodyList[0].m = mTOT
Bridle.BodyList[0].v = VTOT
Bridle.BodyList[0].rCG = rCG_TOT
Bridle.BodyList[0].AWP = AWP_TOT
Bridle.BodyList[0].rM = np.array([0,0,zMeta])
Bridle.BodyList[0].f6Ext = np.array([Fthrust,0,0, 0,Mthrust,0])

Bridle.initialize()






# If using the bridle mooring system rather than the original, do a rename so we can refer to it as MooringSystem going forward (otherwise comment the line out)
MooringSystem = Bridle


# ----------------------------- Calculate mooring system characteristics ---------------------


# First get mooring system characteristics about undisplaced platform position (useful for baseline and verification)


# C_moor0 = MooringSystem.BodyList[0].getStiffness2(np.zeros(6), dx=0.01)  # get mooring stiffness (uses new method that accounts for free Points in mooring system)
# W_moor0 = MooringSystem.BodyList[0].getForces(lines_only=True)           # get net forces and moments from mooring lines on Body


# # Now find static equilibrium offsets of platform and get mooring properties about that point
# MooringSystem.solveEquilibrium()                                        # get the system to its equilibrium
# MooringSystem.plot()
# r6eq = MooringSystem.BodyList[0].r6
# print("Equilibirum platform positions/rotations:")
# printVec(r6eq)
# print("Surge: {:.2f}".format(r6eq[0]))
# print("Pitch: {:.2f}".format(r6eq[4]*180/np.pi))
# C_moor = MooringSystem.BodyList[0].getStiffness2(r6eq, dx=0.01)  # get mooring stiffness (uses new method that accounts for free Points in mooring system)
# W_moor = MooringSystem.BodyList[0].getForces(lines_only=True)           # get net forces and moments from mooring lines on Body

# # manually add yaw spring stiffness as compensation until bridle (crow foot) configuration is added
# #C_moor[5,5] += 98340000.0



# #print(stopt)
# # ------------------------------- sum all static matrices -----------------------------------------
# # this is to check totals from static calculations before hydrodynamic terms are added

# M_tot_stat = M_struc             
# C_tot_stat = C_struc + C_hydro + C_moor
# W_tot_stat = W_struc + W_hydro + W_moor




# print("hydrostatic stiffness matrix")
# printMat(C_hydro)    
    
# print("structural stiffness matrix")
# printMat(C_struc)
    
# print("mooring stiffness matrix")
# printMat(C_moor)
    

# print("total static mass matrix")
# printMat(M_tot_stat)
    
# print("total static stiffness matrix")
# printMat(C_tot_stat)
    
# print("total static forces and moments")
# printVec(W_tot_stat)




# --------------------- get constant hydrodynamic values along each member -----------------------------

A_hydro_morison = np.zeros([6,6])                # hydrodynamic added mass matrix, from only Morison equation [kg, kg-m, kg-m^2]
F_hydro_iner    = np.zeros([6,nw],dtype=complex) # inertia excitation force/moment complex amplitudes vector [N, N-m]


# loop through each member
for mem in memberList:
    
    # loop through each node of the member
    for il in range(mem.n):
    
        # only process hydrodynamics if this node is submerged
        if mem.r[il,2] < 0:
        
            q, p1, p2 = mem.getDirections()                # get unit direction vectors
            
            
            # set dl to half if at the member end (the end result is similar to trapezoid rule integration)
            if il==0 or il==mem.n:
                dl = 0.5*mem.dl
            else:
                dl = mem.dl
            
            
            # get wave kinematics spectra given a certain wave spectrum and location
            mem.u[il,:,:], mem.ud[il,:,:] = getWaveKin(zeta[imeto,:], w, k, depth, mem.r[il,:])    


            # local added mass matrix
            Amat = rho*0.25*np.pi*mem.d[il]**2*dl *( mem.Ca_q*VecVecTrans(q) + mem.Ca_p1*VecVecTrans(p1) + mem.Ca_p2*VecVecTrans(p2) )
            
            
            # add to global added mass matrix for Morison members
            A_hydro_morison += translateMatrix3to6DOF(mem.r[il,:], Amat)
            
            
            # local inertial excitation matrix
            Imat = rho*0.25*np.pi*mem.d[il]**2*dl *( (1+mem.Ca_q)*VecVecTrans(q) + (1+mem.Ca_p1)*VecVecTrans(p1) + (1+mem.Ca_p2)*VecVecTrans(p2) )
            
            for i in range(nw):   # for each wave frequency...
            
                # local inertial excitation force complex amplitude in x,y,z
                F_exc_inert = np.matmul(Imat, mem.ud[il,:,i])  
            
                # add to global excitation vector (frequency dependent)
                F_hydro_iner[:,i] += translateForce3to6DOF( mem.r[il,:], F_exc_inert)




# --------------------------------- get system properties in undisplaced position ----------------------------
# these are useful for verification, etc.


# sum matrices to check totals from static calculations before hydrodynamic terms are added

C_tot0 = C_struc + C_hydro # + C_moor0   # total system stiffness matrix about undisplaced position
W_tot0 = W_struc + W_hydro #+ W_moor0   # system mean forces and moments at undisplaced position

M = M_struc + A_hydro_morison          # total mass plus added mass matrix

# do we want to relinearize structural properties about displaced position/orientation?  (Probably not)

'''
print("hydrostatic stiffness matrix")
printMat(C_hydro)    
    
print("structural stiffness matrix")
printMat(C_struc)
    
print("mooring stiffness matrix about undisplaced position")
printMat(C_moor0)
    
print("total static stiffness matrix about undisplaced position")
printMat(C_tot0)
    
    

print("total static mass matrix")
printMat(M_struc)
    
print("total added mass matrix")
printMat(A_hydro_morison)

print("total mass plus added mass matrix")
printMat(M)
    
    
print("total static forces and moments about undisplaced position")
printVec(W_tot0)
'''



# calculate natural frequencies (using eigen analysis to get proper values for pitch and roll - otherwise would need to base about CG if using diagonal entries only)

eigenvals, eigenvectors = np.linalg.eig(np.matmul(np.linalg.inv(M), C_tot0))   # <<< need to sort this out so it gives desired modes, some are currently a bit messy


# alternative attempt to calculate natural frequencies based on diagonal entries (and taking pitch and roll about CG)

zMoorx = C_tot0[0,4]/C_tot0[0,0]  # effective z elevation of mooring system reaction forces in x and y directions
zMoory = C_tot0[1,3]/C_tot0[1,1]
zCG  = rCG_TOT[2]                    # center of mass in z
zCMx = M[0,4]/M[0,0]              # effective z elevation of center of mass and added mass in x and y directions
zCMy = M[1,3]/M[1,1]

print("natural frequencies without added mass")
fn = np.zeros(6)
fn[0] = np.sqrt( C_tot0[0,0] / M_struc[0,0] )/ 2.0/np.pi
fn[1] = np.sqrt( C_tot0[1,1] / M_struc[1,1] )/ 2.0/np.pi
fn[2] = np.sqrt( C_tot0[2,2] / M_struc[2,2] )/ 2.0/np.pi
fn[5] = np.sqrt( C_tot0[5,5] / M_struc[5,5] )/ 2.0/np.pi
zg = rCG_TOT[2]
fn[3] = np.sqrt( (C_tot0[3,3] + C_tot0[1,1]*((zCG-zMoory)**2 - zMoory**2) ) / (M_struc[3,3] - M_struc[1,1]*zg**2 ))/ 2.0/np.pi     # this contains adjustments to reflect rotation about the CG rather than PRP
fn[4] = np.sqrt( (C_tot0[4,4] + C_tot0[0,0]*((zCG-zMoorx)**2 - zMoorx**2) ) / (M_struc[4,4] - M_struc[0,0]*zg**2 ))/ 2.0/np.pi     # this contains adjustments to reflect rotation about the CG rather than PRP
printVec(fn)


print("natural frequencies with added mass")
fn = np.zeros(6)
fn[0] = np.sqrt( C_tot0[0,0] / M[0,0] )/ 2.0/np.pi
fn[1] = np.sqrt( C_tot0[1,1] / M[1,1] )/ 2.0/np.pi
fn[2] = np.sqrt( C_tot0[2,2] / M[2,2] )/ 2.0/np.pi
fn[5] = np.sqrt( C_tot0[5,5] / M[5,5] )/ 2.0/np.pi
fn[3] = np.sqrt( (C_tot0[3,3] + C_tot0[1,1]*((zCMy-zMoory)**2 - zMoory**2) ) / (M[3,3] - M[1,1]*zCMy**2 ))/ 2.0/np.pi     # this contains adjustments to reflect rotation about the CG rather than PRP
fn[4] = np.sqrt( (C_tot0[4,4] + C_tot0[0,0]*((zCMx-zMoorx)**2 - zMoorx**2) ) / (M[4,4] - M[0,0]*zCMx**2 ))/ 2.0/np.pi     # this contains adjustments to reflect rotation about the CG rather than PRP
# note that the above lines use off-diagonal term rather than parallel axis theorem since rotation will not be exactly at CG due to effect of added mass
printVec(fn)





# ------------------- solve for platform dynamics, iterating until convergence --------------------

Z  = np.zeros([6,6,nw], dtype=complex)               # total system impedance matrix

# system response 
Xi = np.zeros([6,nw], dtype=complex) + 0.01    # displacement and rotation complex amplitudes [m, rad]



nIter = 2  # maximum number of iterations to allow

# start fixed point iteration loop for dynamics
for iiter in range(nIter):


    # ------ calculate linearized coefficients within iteration ------- 
    
    B_hydro_drag = np.zeros([6,6])             # hydrodynamic damping matrix (just linearized viscous drag for now) [N-s/m, N-s, N-s-m]

    F_hydro_drag = np.zeros([6,nw],dtype=complex) # excitation force/moment complex amplitudes vector [N, N-m]


    # loop through each member
    for mem in memberList:
        
        q, p1, p2 = mem.getDirections()                # get unit direction vectors
        
        # loop through each node of the member
        for il in range(mem.n):
            
            # node displacement, velocity, and acceleration (each [3 x nw])
            drnode, vnode, anode = getVelocity(mem.r[il,:], Xi, w)      # get node complex velocity spectrum based on platform motion's and relative position from PRP
            
            
            # only process hydrodynamics if this node is submerged
            if mem.r[il,2] < 0:
            
                # water relative velocity over node (complex amplitude spectrum)  [3 x nw]
                vrel = mem.u[il,:] - vnode
                
                # break out velocity components in each direction relative to member orientation [nw]
                vrel_q  = vrel* q[:,None]
                vrel_p1 = vrel*p1[:,None]
                vrel_p2 = vrel*p2[:,None]
                
                # get RMS of relative velocity component magnitudes (real-valued)
                vRMS_q  = np.linalg.norm( np.abs(vrel_q ) )  # equivalent to np.sqrt( np.sum( np.abs(vrel_q )**2) /nw)
                vRMS_p1 = np.linalg.norm( np.abs(vrel_p1) )
                vRMS_p2 = np.linalg.norm( np.abs(vrel_p2) )
                
                # linearized damping coefficients in each direction relative to member orientation [not explicitly frequency dependent...] (this goes into damping matrix)
                Bprime_q  = np.sqrt(8/np.pi) * vRMS_q  * 0.5*rho * np.pi*mem.d[il]*mem.dl * mem.Cd_q 
                Bprime_p1 = np.sqrt(8/np.pi) * vRMS_p1 * 0.5*rho * mem.d[il]*mem.dl * mem.Cd_p1
                Bprime_p2 = np.sqrt(8/np.pi) * vRMS_p2 * 0.5*rho * mem.d[il]*mem.dl * mem.Cd_p2
                
                # convert to global orientation
                Bmat = Bprime_q*VecVecTrans(q) + Bprime_p1*VecVecTrans(p1) + Bprime_p2*VecVecTrans(p2)
                
                # add to global damping matrix for Morison members
                Btemp = translateMatrix3to6DOF(mem.r[il,:], Bmat)
                
                #breakpoint()
                
                B_hydro_drag += Btemp
                
                
                
                # excitation force based on linearized damping coefficients [3 x nw]
                F_exc_drag = np.zeros([3, nw], dtype=complex)  # <<< should set elsewhere <<<
                for i in range(nw):

                    # get local 3d drag excitation force complex amplitude for each frequency
                    F_exc_drag[:,i] = np.matmul(Bmat, mem.u[il,:,i])  
                
                    # add to global excitation vector (frequency dependent)
                    F_hydro_drag[:,i] += translateForce3to6DOF( mem.r[il,:], F_exc_drag[:,i])




    # ----------------------------- solve matrix equation of motion ------------------------------
    
    for ii in range(nw):          # loop through each frequency
        
        
        # sum contributions for each term        
        M_tot[:,:,ii] = M_struc + A_hydro_morison                         # mass
        B_tot[:,:,ii] = B_struc + B_hydro_drag                            # damping
        C_tot[:,:,ii] = C_struc + C_hydro# + C_moor                        # stiffness
        F_tot[:,  ii] = F_hydro_drag[:,ii] + F_hydro_iner[:,ii]           # excitation force (complex amplitude)
        
        
        # form impedance matrix
        Z[:,:,ii] = -w[ii]**2 * M_tot[:,:,ii] + 1j*w[ii]*B_tot[:,:,ii] + C_tot[:,:,ii];
        
        # solve response (complex amplitude)
        Xi[:,ii] = np.matmul(np.linalg.inv(Z[:,:,ii]),  F_tot[:,ii] )
    
    
    '''

    #Xi{imeto} = rao{imeto}.*repmat(sqrt(S(:,imeto)),1,6); # complex response!
    
    #aNacRAO{imeto} = -(w').^2 .* (rao{imeto}(:,1) + hNac*rao{imeto}(:,5));      # Nacelle Accel RAO
    #aNac2(imeto) = sum( abs(aNacRAO{imeto}).^2.*S(:,imeto) ) *(w(2)-w(1));     # RMS Nacelle Accel

    
    # ----------------- convergence check --------------------
    conv = abs(aNac2(imeto)/aNac2last - 1);
    #disp(['at ' num2str(iiter) ' iterations - convergence is ' num2str(conv)])
    if conv < 0.0001

         
         break
    else
         aNac2last = aNac2(imeto);
    end
    
    '''

# ------------------------------ preliminary plotting of response ---------------------------------

fig, ax = plt.subplots(3,1, sharex=True)

ax[0].plot(w, np.abs(Xi[0,:])          , 'b', label="surge")
ax[0].plot(w, np.abs(Xi[1,:])          , 'g', label="sway")
ax[0].plot(w, np.abs(Xi[2,:])          , 'r', label="heave")
ax[1].plot(w, np.abs(Xi[3,:])*180/np.pi, 'b', label="roll")
ax[1].plot(w, np.abs(Xi[4,:])*180/np.pi, 'g', label="pitch")
ax[1].plot(w, np.abs(Xi[5,:])*180/np.pi, 'r', label="yaw")
ax[2].plot(w, zeta[0,:],                  'k', label="wave amplitude (m)")

ax[0].legend()
ax[1].legend()
ax[2].legend()

#ax[0].set_ylim([0, 1e6])
#ax[1].set_ylim([0, 1e9])

ax[0].set_ylabel("response magnitude (m)")
ax[1].set_ylabel("response magnitude (deg)")
ax[2].set_ylabel("wave amplitude (m)")
ax[2].set_xlabel("frequency (Hz)")

plt.show()

 
 # ---------- mooring line fairlead tension RAOs and constraint implementation ----------
'''
 
 for il=1:Platf.Nlines
      
      #aNacRAO{imeto} = -(w').^2 .* (X{imeto}(:,1) + hNac*X{imeto}(:,5));      # Nacelle Accel RAO
        #aNac2(imeto) = sum( abs(aNacRAO{imeto}).^2.*S(:,imeto) ) *(w(2)-w(1));     # RMS Nacelle Accel

    TfairRAO{imeto}(il,:) = C_lf(il,:,imeto)*rao{imeto}(:,:)';  # get fairlead tension RAO for each line (multiply by dofs)
      #RMSTfair{imeto}(il) = sqrt( sum( (abs(TfairRAO{imeto}(il,:))).^2) / length(w) );
      #figure
    #plot(w,abs(TfairRAO{imeto}(il,:)))
      #d=TfairRAO{imeto}(il,:)
      RMSTfair{imeto}(il) = sqrt( sum( (abs(TfairRAO{imeto}(il,:)).^2).*S(:,imeto)') *(w(2)-w(1)) );
      #RMSTfair
      #sumpart = sum( (abs(TfairRAO{imeto}(il,:)).^2).*S(:,imeto)')
      #dw=(w(2)-w(1))
 end        
 
 [Tfair, il] = min( T_lf(:,imeto) );
 if Tfair - 3*RMSTfair{imeto}(il) < 0 && Xm < 1  # taut lines only
      disp([' REJECTING (mooring line goes slack)'])
      fitness = -1;
      return;  # constraint for slack line!!!
 end
 if grads
      disp(['mooring slackness: ' num2str(Tfair - 3*RMSTfair{imeto}(il))])
 end
 
 # ----------- dynamic pitch constraint ----------------------
 #disp('checking dynamic pitch');
 RMSpitch(imeto) = sqrt( sum( ((abs(rao{imeto}(:,5))).^2).*S(:,imeto) ) *(w(2)-w(1)) ); # fixed April 9th :(
 RMSpitchdeg = RMSpitch(imeto)*60/pi;
 if (Platf.spitch + RMSpitch(imeto))*180/pi > 10
      disp([' REJECTING (static + RMS dynamic pitch > 10)'])
      fitness = -1;
      return;  
 end    
 if grads
      disp(['dynamic pitch: ' num2str((Platf.spitch + RMSpitch(imeto))*180/pi)])
 end
 
 #figure(1)
 #plot(w,S(:,imeto))
 #hold on
 
 #figure()
 #plot(2*pi./w,abs(Xi{imeto}(:,5)))
 #ylabel('pitch response'); xlabel('T (s)')
 
 RMSsurge(imeto) = sqrt( sum( ((abs(rao{imeto}(:,1))).^2).*S(:,imeto) ) *(w(2)-w(1)) ); 
 RMSheave(imeto) = sqrt( sum( ((abs(rao{imeto}(:,3))).^2).*S(:,imeto) ) *(w(2)-w(1)) ); 
 
''' 
    