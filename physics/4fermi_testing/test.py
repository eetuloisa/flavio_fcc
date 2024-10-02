#from flavio.math.integrate import nintegrate
#from flavio.physics.zdecays.smeftew import gV_SM, gA_SM, _QN
#from flavio.physics.common import add_dict
#from flavio.classes import Observable, Prediction
#from flavio.physics import ckm as ckm_flavio
#import wilson
from numpy import pi, sqrt
import numpy as np
#from flavio.config import config


def sigma_llqq_tot(sh,l,q1,q2,wc_eff):
    r"""Total cross section of e^+ e^- \to q1 q2

    Returns $\sigma$ in units of GeV$^{-2}$

    Parameters:
    - `sh`: centre of mass energy in GeV$^2$ of the lepton pair
    - `l`: lepton flavour, should be 'e'
    - `q1`, `q2`: outgoing quark flavours
    - `wc_eff`: SM + SMEFT amplitude for the process as dictionary with entries corresponding to different Wilson coefficients
    """

    # this function returns the piece coming out of the integration of the Dirac traces
    def sh_func(sh, coefficient):
        return coefficient * sh**3

    # Define a dictionary of the non-zero Dirac traces when integrating over the full phase space
    f_fully_integrated = {}
    f_fully_integrated['SSLL'] = f_fully_integrated['SSLR'] = f_fully_integrated['SSRL'] = f_fully_integrated['SSRR'] = lambda sh: sh_func(sh, 1/4)
    f_fully_integrated['VVLL'] = f_fully_integrated['VVLR'] = f_fully_integrated['VVRL'] = f_fully_integrated['VVRR'] = lambda sh: sh_func(sh, 1/3)
    f_fully_integrated['TTLR'] = f_fully_integrated['TTRL'] = lambda sh: sh_func(sh, 4/3)

    # Calculate the cross-section sigma
    sigma = 0
    for gamma in 'SVT':
        for gammap in 'SVT':
            for X in 'LR':
                for Y in 'LR':
                    key = '{}{}{}{}'.format(gamma, gammap, X, Y)
                    if key in f_fully_integrated:
                        sigma += 3 / (16 * pi * sh**2) * f_fully_integrated['{}{}{}{}'.format(gamma, gammap, X, Y)](sh) #* wc_eff['C{}{}{}_e{}'.format(gamma,X,Y,q1)][0,0,q1,q2] * np.conjugate(wc_eff['C{}{}{}_e{}'.format(gammap,X,Y,q1)][0,0,q1,q2])

    # Return total cross-section
    return sigma

