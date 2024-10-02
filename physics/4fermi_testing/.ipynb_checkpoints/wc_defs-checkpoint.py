# This .py file is needed because multiprocessing suddenly doesn't work within jupyter lab.
# So the functions that we want to run in smelli, for instance when using the built-in plotting functionality, have to be defined here.

import numpy as np 
import cmath
from math import e
import matplotlib.pyplot as plt
import numpy as np
#import seaborn


# Angles for VdL. We are using the ansatz that these are equal to the measured
# central values of the CKM matrix except for theta_23 which we allow to float.
# From PDG 2020 except s23
s12 = 0.22650
s13 = 0.00361
delta = 1.196
t12 = np.arcsin(s12)
t13 = np.arcsin(s13)
c12 = np.cos(t12)
c13 = np.cos(t13)
# Note: this means that nearly all CKM mixing is contained by V_{d_L} --- it is only the V_{d_L}_23 which is allowed to float (it controls contribution to C9)
# V_{u_L} is then tuned to cancel out with the free-floating rotation to give the correct CKM matrix
# see between Eqs 10 and 11 in 2103.12056v3 
m13 = np.array([
    [                      c13,      0,   s13*e**complex(0, -delta)],
    [                        0,      1,                           0],
    [-s13*e**complex(0, delta),      0,                         c13]
])
m12 = np.array([
    [ c12,  s12,   0],
    [-s12,  c12,   0],
    [   0,    0,   1]
])
mmm = m13 @ m12 # Matrix multiply the two
l_xi = []
# This is a convenience function that returns the indices with an off-set
# it helps with checking and makes the writing easier
def L_xi(i,j): 
    return l_xi[i-1][j-1]
    #return "L_xi"+str(i-1)+str(j-1)

# For four-fermion operators, combinations like 2212 and 2221 correspond to the same operator - not every single combination is independent
lq1_nonredundant_indices = [2211,2212,2213,2222,2223,2233,3311,3312,3313,3322,3323,3333]

qq1_nonredundant_indices = [1111,1112,1113,1122,1123,1133,1212,1213,1221,1222,1223,1231,1232,1233,1313,1322,1323,1331,1332,1333,2222,2223,2233,2323,2332,2333,3333]

def fill_qq1(alpha):
   qq1_3333_unmixed = -(1/72) * alpha 
   listqq1 = {}
   for i in range(1,4):
       for j in range(1,4):
           for k in range(1,4):
               for l in range(1,4): # Loops over the four indices
                   # The indices as an int
                   intIndex = 1000*i + 100*j + 10*k + l 
                   # Combination from commuting currents
                   intpermIndex = 1000*k + 100*l + 10*i + j 
                   # h.c
                   inthcIndex1 = 1000*l + 100*k + 10*j + i 
                   # h.c + commuting currents
                   inthcIndex2 = 1000*j + 100*i + 10*l + k 
                   # List of equivalent permutations
                   permlist = [intIndex,intpermIndex, inthcIndex1 , inthcIndex2] 
                   # Here I check which combination is the one in the Smelli basis
                   match = [m for m in qq1_nonredundant_indices if m in permlist][0] 
                   name = 'qq1_' + str(match)
                   factor = 2
                   if i==k and j==l:
                        factor = 1
                   # Here I add the contribution to the non-redundant op.
                   if name in listqq1:
                       listqq1[name] += factor * qq1_3333_unmixed * L_xi(i,j) * L_xi(k,l)
                   else:
                       listqq1[name] = factor * qq1_3333_unmixed * L_xi(i,j) * L_xi(k,l)
   return listqq1

VdL = []   

# This function is for setting the input: alpha=gz'^2/MZ'^2, theta=theta_sb input too
# The down quarks are in their mass basis for these unhatted operators, whereas the up
# quarks will be CKM-mixed
# 19.4.24: It's unclear what the comment above meant. Does it say that we're using smelli in the Warsaw-down basis?

def calcinputlist(alpha, theta_23, sphi=0, m_theta = 3000):
        #flavon parameters
        
        # Here a decision needs to be made: to keep v_theta constant or let it change as we vary gzp?
        # 1st option:
        #v_theta = 3000
        # 2nd option:
        if (alpha > 1.0e-12):
            v_theta = 1 / (np.sqrt(alpha))
        else:
            v_theta = 1.0e6 # This is essentially just to avoid a divide by zero error when alpha=0
        # 1st option corresponds to changing q_th in parallel with gzp, such that the v_theta is constant
    
        # 19.4.24: I now think the second option probably makes more sense - MZp ~ v_theta is very natural
        
        # Fix remaining dimensionful parameters
        #m_theta = 1000
        v_higgs = 246
        m_higgs = 125
        
        #Calculate Lambda_{H \theta} using the input parameters
        lam11 = np.sqrt(1 - (1 - 2 * sphi**2)**2) * (m_higgs**2 - m_theta**2) / (2 * v_theta * v_higgs)
        print('lam11 = ', lam11)
        
        # Left-handed down-type quark mixing matrix
        s23 = np.sin(theta_23)
        c23 = np.cos(theta_23)
        m23 = np.array([
            [1,    0,   0],
            [0,  c23, s23],
            [0, -s23, c23]
        ])
        VdL = m23 @ mmm
        # matrix of couplings to down quarks: see (18) in 1905.10327
        l_xi.clear()
        for i in range(0, 3):
            l_xi.append([np.conj(VdL[2, i]) * VdL[2, k] for k in range(0, 3)])# Make central-valued CKM matrix: from PDG 2020      
        inputlist_qq1 = fill_qq1(alpha)   
        # Checked by BCA 1/3/21
        inputlist_explicit = {            
            # LL operators
            'll_2222' : -(1/8)  * alpha,
            'lq1_2211': (1/12)  * alpha * l_xi[0][0],
            'lq1_2212': (1/12)  * alpha * (l_xi[0][1] + l_xi[1][0]),
            'lq1_2213': (1/12)  * alpha * (l_xi[0][2] + l_xi[2][0]),
            'lq1_2222': (1/12)  * alpha * l_xi[1][1],
            'lq1_2223': (1/12)  * alpha * (l_xi[1][2] + l_xi[2][1]),
            'lq1_2233': (1/12)  * alpha * l_xi[2][2], 
            # RR operators
            'ee_3333': -(1/2)  * alpha,
            'uu_3333': -(2/9)  * alpha,
            'dd_3333': -(1/18) * alpha,
            'eu_3333':  (2/3)  * alpha,
            'ed_3333': -(1/3)  * alpha,
            'ud1_3333': (2/9)  * alpha,
            # LR operators
            'le_2233': -(1/2)  * alpha,
            'lu_2233':  (1/3)  * alpha,
            'ld_2233': -(1/6)  * alpha,
            'qe_1133':  (1/6)   * alpha * l_xi[0][0] ,
            'qe_1233':  (1/6)   * alpha * (l_xi[0][1] + l_xi[1][0]),
            'qe_1333':  (1/6)   * alpha * (l_xi[0][2] + l_xi[2][0]),
            'qe_2233':  (1/6)   * alpha * l_xi[1][1],
            'qe_2333':  (1/6)   * alpha * (l_xi[1][2] + l_xi[2][1]),
            'qe_3333':  (1/6)   * alpha * l_xi[2][2],            
            'qu1_1133':-(1/9)   * alpha * l_xi[0][0] ,
            'qu1_1233':-(1/9)   * alpha * (l_xi[0][1] + l_xi[1][0]),
            'qu1_1333':-(1/9)   * alpha * (l_xi[0][2] + l_xi[2][0]),
            'qu1_2233':-(1/9)   * alpha * l_xi[1][1],
            'qu1_2333':-(1/9)   * alpha * (l_xi[1][2] + l_xi[2][1]),
            'qu1_3333':-(1/9)   * alpha * l_xi[2][2],
            'qd1_1133': (1/18)  * alpha * l_xi[0][0],
            'qd1_1233': (1/18)  * alpha * (l_xi[0][1] + l_xi[1][0]),
            'qd1_1333': (1/18)  * alpha * (l_xi[0][2] + l_xi[2][0]),
            'qd1_2233': (1/18)  * alpha * l_xi[1][1],
            'qd1_2333': (1/18)  * alpha * (l_xi[1][2] + l_xi[2][1]),
            'qd1_3333': (1/18)  * alpha * l_xi[2][2],
            # Higgs + two fermions NB WE MUST CHECK THE SIGNS OF THESE CONTRIBUTIONS!!!
                    # LL operators
            'phil1_22': +(1/4)  * alpha,
            'phiq1_11': -(1/12)  * alpha * l_xi[0][0],
            'phiq1_12': -(1/12)  * alpha * (l_xi[0][1] + l_xi[1][0]),
            'phiq1_13': -(1/12)  * alpha * (l_xi[0][2] + l_xi[2][0]),
            'phiq1_22': -(1/12)  * alpha * l_xi[1][1],
            'phiq1_23': -(1/12)  * alpha * (l_xi[1][2] + l_xi[2][1]),
            'phiq1_33': -(1/12)  * alpha * l_xi[2][2],           
            'phie_33':   (1/2)  * alpha,
            'phiu_33':  -(1/3)  * alpha,
            'phid_33':   (1/6)  * alpha,
            'phiD':     -(1/2)  * alpha,
            'phiBox':   -(1/8)  * alpha - lam11**2 * v_theta**2 / (2 * m_theta**4) # The only operator that the flavon contributes to
        } 
        return {**inputlist_explicit, **inputlist_qq1}
