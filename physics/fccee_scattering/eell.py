r"""Functions for the process $e^+e^- \to l^+l^-$ for l != e."""
# Written by Eetu Loisa, 2024
from flavio.physics.zdecays.smeftew import gV_SM, gA_SM
import flavio.physics.zdecays.smeftew as smeftew
from flavio.physics.common import add_dict
from flavio.classes import Observable, Prediction
from flavio.physics import ckm as ckm_flavio
#from . import partondist
import wilson
from numpy import pi, sqrt
import numpy as np
from flavio.config import config


def F_eell_SM(l1, Xl1, l2, Xl2, s, wc_obj, par):
    r"""SM $Z$ and $\gamma$ contribution to the $\bar q q\to \ell^+\ell^-$
    amplitude."""
    # This function is adapted from flavio.physics.dileptons.ppll.py 
    E = np.sqrt(s)
    wcxf = wc_obj.get_wcxf(sector='dB=dL=0', scale=E, par=par, eft='SMEFT', basis='Warsaw up')

    # parameters
    Ql = -1
    mZ = par['m_Z']
    GammaZ = 1 / par['tau_Z']
    s2w = par['s2w']
    aEW = par['alpha_e']
    e2 = 4 * pi * aEW
    g2=e2/s2w 
    gp2=e2/(1-s2w)  # We do not run the couplings as even at 365 GeV they change by less than 1%
    g_cw = sqrt(g2+gp2) # g over cw
    s2w = gp2/(g2+gp2)
    par = par.copy()
    par['s2w'] = s2w
    # SM couplings
    gVl1 = g_cw * gV_SM(l1, par)
    gVl2 = g_cw * gV_SM(l2, par)
    gAl1 = g_cw * gA_SM(l1, par)
    gAl2 = g_cw * gA_SM(l2, par)
    # Add SMEFT corrections to couplings
    gVl1 += g_cw * smeftew.d_gVl(l1, l1, par, wcxf)
    gVl2 += g_cw * smeftew.d_gVl(l2, l2, par, wcxf)
    gAl1 += g_cw * smeftew.d_gAl(l1, l1, par, wcxf)
    gAl2 += g_cw * smeftew.d_gAl(l2, l2, par, wcxf)

    if Xl1 == 'L':
        gZl1 = gVl1 + gAl1
    elif Xl1 == 'R':
        gZl1 = gVl1 - gAl1
    if Xl2 == 'L':
        gZl2 = gVl2 + gAl2
    elif Xl2 == 'R':
        gZl2 = gVl2 - gAl2
    # SM contribution
    return g2*s2w * Ql * Ql / s + gZl1 * gZl2 / (s - mZ**2 + 1j * mZ * GammaZ)

# This is sm-like in the sense that it only uses SM-like vertices; but the Zff couplings are modified by SMEFT corrections
def wceff_eell_sm(wc_obj, par, s):
    wc = {}
    for Xl1 in 'LR':
        for Xl2 in 'LR':
            wc['CV{}{}_mu'.format(Xl1, Xl2)] = F_eell_SM('e', Xl1, 'mu', Xl2, s, wc_obj, par) 
            wc['CV{}{}_tau'.format(Xl1, Xl2)] = F_eell_SM('e', Xl1, 'tau', Xl2, s, wc_obj, par) 
            #print('Added CV{}{}_l to wc'.format(Xl1, Xl2) + ' with value' + str(wc['CV{}{}_mu'.format(Xl1, Xl2)]) + 'at s =' + str(s))
    return wc


def wceff_eell_np(wc_obj, par, scale):
    r"""Returns Wilson coefficients of the effective Lagrangian
    $$\mathcal{L} =  C_{el}^{\Gamma_1\Gamma_2}
    (\bar e_i\Gamma_1 e_j)(\bar e_k\Gamma_2 e_l)$$
    """
    # get the dictionary, in a non-redundant basis
    wcxf_dict = wc_obj.get_wcxf(sector='dB=dL=0', scale=scale, par=par, eft='SMEFT', basis='Warsaw up').dict
    wc = {}
    # match to notation used in the observable
    # Note that the approach here is slightly different from flavio.physics.dileptons.ppll.py in that we add hermitian conjugates explicitly into the wc dictionary
    # CVXX_l (8 entries, 4 for each lepton):
    wc['CVLL_mu'] = wcxf_dict.get('ll_1122',0) +  wcxf_dict.get('ll_1221',0)
    wc['CVRR_mu'] = wcxf_dict.get('ee_1122',0)
    wc['CVLR_mu'] = wcxf_dict.get('le_1122',0)
    wc['CVRL_mu'] = wcxf_dict.get('le_2211',0)

    wc['CVLL_tau'] = wcxf_dict.get('ll_1133',0) +  wcxf_dict.get('ll_1331',0)
    wc['CVRR_tau'] = wcxf_dict.get('ee_1133',0)
    wc['CVLR_tau'] = wcxf_dict.get('le_1133',0)
    wc['CVRL_tau'] = wcxf_dict.get('le_3311',0)

    # These are the SMEFT coefficients which contribute to the scattering cross-section but which do not interfere with the SM due to the chirality structure
    wc['C_non_interfering_mu'] = wcxf_dict.get('le_1221',0)
    wc['C_non_interfering_tau'] = wcxf_dict.get('le_1331',0)

    return wc


# translate quark name to LHAPDF flavour index
fermion_indices = {
    'mu': ('e', 1),
    'tau': ('e', 2),
}

# this function returns the piece coming out of the integration of the Dirac traces
def s_func(s, coefficient):
    return coefficient * s**3


# Define a dictionary of the non-zero Dirac traces when integrating over cos(theta) in the range [-0.9, 0.9] 
f_integrated = {}
f_integrated['VVLL_tau'] =  f_integrated['VVRR_tau'] = lambda s: s_func(s, 1143/4000)
f_integrated['VVLR_tau'] = f_integrated['VVRL_tau'] = lambda s: s_func(s, 1143/4000)
f_integrated['non_interfering_tau'] = lambda s: s_func(s, 9/10) # This arises from the SMEFT coefficient C_le_1331 which does not interfere with the SM


# Define a dictionary of the non-zero Dirac traces when integrating over cos(theta) in the range [0, 0.9] MINUS the range [-0.9, 0], to get A_FB
f_integrated_AFB = {}
f_integrated_AFB['VVLL_tau'] = f_integrated_AFB['VVRR_tau'] = lambda s: s_func(s, 81/400)
f_integrated_AFB['VVRL_tau'] = f_integrated_AFB['VVLR_tau'] = lambda s: s_func(s, -81/400)


# Define a dictionary of the non-zero Dirac traces when integrating over cos(theta) in the range [-0.95, 0.95] 
f_integrated['VVLL_mu'] =  f_integrated['VVRR_mu'] = lambda s: s_func(s, 29659/96000)
f_integrated['VVLR_mu'] = f_integrated['VVRL_mu'] = lambda s: s_func(s, 29659/96000)
f_integrated['non_interfering_mu'] = lambda s: s_func(s, 19/20) # This arises from the SMEFT coefficient C_le_1221 which does not interfere with the SM


# Define a dictionary of the non-zero Dirac traces when integrating over cos(theta) in the range [0, 0.95] MINUS the range [-0.95, 0], to get A_FB
f_integrated_AFB['VVLL_mu'] = f_integrated_AFB['VVRR_mu'] = lambda s: s_func(s, 361/1600)
f_integrated_AFB['VVRL_mu'] = f_integrated_AFB['VVLR_mu'] = lambda s: s_func(s, -361/1600)


def sigma_eell_tot(wc_dict, s, l):
    r"""Total cross section of e^+ e^- \to \ell^+ \ell^-

    Returns $\sigma$ in units of GeV$^{-2}$

    Parameters:
    - `s`: centre of mass energy in GeV$^2$ of the electron pair
    - `l`: lepton flavour, should be 'mu' or 'tau'
    - `wc_dict`: dictionary containing the effective SM + SMEFT Wilson coefficients contributing to lepton scattering
    """

    # Calculate the cross-section, sigma
    sigma = 0
    for X in 'LR':
        for Y in 'LR':
            key = 'VV{}{}_{}'.format(X, Y, l)
            if key in f_integrated:
                sigma += 1 / (16 * pi * s**2) * f_integrated['VV{}{}_{}'.format(X, Y, l)](s) * np.abs(wc_dict['CV{}{}_{}'.format(X,Y,l)])**2 
                #print('Added', key, 'to total cross-section')
                #print('sigma =', sigma, 'at this point')
    # Return total cross-section
    sigma += 1 / (16 * pi * s**2) * f_integrated['non_interfering_{}'.format(l)](s) * np.abs(wc_dict['C_non_interfering_{}'.format(l)])**2
    return sigma

def sigma_eell_forward_minus_backward(wc_dict, s, l):
    r"""Cross section of e^+ e^- \to \ell^+ \ell^- in the forward region minus the backward region

    Returns $\sigma_F - \sigma_B$ in units of GeV$^{-2}$

    Parameters:
    - `s`: centre of mass energy in GeV$^2$ of the lepton pair
    - `l`: lepton flavour, should be 'mu' or 'tau'
    - `wc_dict`: SM + SMEFT amplitude for the process as dictionary with entries corresponding to different Wilson coefficients
    """

    #q, i = fermion_indices[q]

    # Calculate the cross-section in the forward minus backward region, sigma_F - sigma_B
    sigma = 0
    for X in 'LR':
        for Y in 'LR':
            key = 'VV{}{}_{}'.format(X, Y, l)
            if key in f_integrated:
                sigma += 1 / (16 * pi * s**2) * f_integrated_AFB['VV{}{}_{}'.format(X, Y, l)](s) * np.abs(wc_dict['CV{}{}_{}'.format(X,Y,l)])**2 
                #print('Added', key, 'to F minus B')
                #print(' =', sigma, 'at this point')
    # Return total cross-section
    return sigma

def sigma_eell_tot_obs(wc_obj, par, E, l):
    # This function takes the BSM Wilson coefficients as input, runs them down, converts to wc_eff and adds to the SM contribution, and finally returns the cross-section
    scale = E
    s = E**2
    wc_eff_np = wceff_eell_np(wc_obj, par, scale)
    #print('wc_eff_np =', wc_eff_np)
    wc_eff_sm = wceff_eell_sm(wc_obj, par, s)
    wc_eff = add_dict((wc_eff_sm, wc_eff_np))
    #print('wc_eff keys' + str(wc_eff.keys()))

    return np.real(sigma_eell_tot(wc_eff, s, l))

def AFB_eell_obs(wc_obj, par, E, l):
    # This function takes the BSM Wilson coefficients as input, runs them down, converts to wc_eff and adds to the SM contribution, and finally returns the forward-backward asymmetry
    scale = E
    s = E**2
    wc_eff_np = wceff_eell_np(wc_obj, par, scale)
    #print('wc_eff_np =', wc_eff_np)
    wc_eff_sm = wceff_eell_sm(wc_obj, par, s)
    wc_eff = add_dict((wc_eff_sm, wc_eff_np))
    #print('wc_eff keys' + str(wc_eff.keys()))

    A_FB = sigma_eell_forward_minus_backward(wc_eff, s, l) / sigma_eell_tot(wc_eff, s, l)

    return np.real(A_FB)

# Function to help generate the total cross-section observable for a given quark q
def generate_sigma_obs(l):
    _selection_efficiencies = {'mu': 0.98, 'tau': 0.9}

    def f(wc_obj, par, E):
        return _selection_efficiencies[l] * sigma_eell_tot_obs(wc_obj, par, E, l)
    return f

# Function to help generate the total cross-section observable for a given quark q
def generate_afb_obs(l):
    def f(wc_obj, par, E):
        return AFB_eell_obs(wc_obj, par, E, l)
    return f

# Observable and Prediction instances
_tex = ['mu','tau']

# Create the observables and predictions for the total cross-section
for l in _tex:
    _process_tex = r"e^+e^- \to \overline{" + l + r"}" + l 
    _process_taxonomy = r'Process :: $e^+e^-$ scattering :: $' + _process_tex + r"$" #TODO Figure out how this should be formatted
    _obs_name = "sigma(ee->{}{})(high_E)".format(l,l)
    _obs = Observable(_obs_name, arguments=['E'])
    _obs.set_description(r"Cross section of $" + _process_tex + r"$")
    _obs.tex = r"$\sigma(" + _process_tex + r")$"
    _obs.add_taxonomy(_process_taxonomy)
    
    Prediction(_obs_name, generate_sigma_obs(l))
    #print("Prediction for ", _obs_name, "added")

# Create the observables and predictions for the forward-backward asymmetry
for l in _tex:
    _process_tex = r"A_{FB}(e^+e^- \to \overline{" + l + r"}" + l + r")"
    _process_taxonomy = r'Process :: $e^+e^-$ scattering :: $' + _process_tex + r"$" #TODO Figure out how this should be formatted
    _obs_name = "AFB(ee->{}{})(high_E)".format(l,l)
    _obs = Observable(_obs_name, arguments=['E'])
    _obs.set_description(r"Forward-backward asymmetry above Z pole $" + _process_tex + r"$")
    _obs.tex = r"$" + _process_tex + r"$"
    _obs.add_taxonomy(_process_taxonomy)
    
    Prediction(_obs_name, generate_afb_obs(l))
    #print("Prediction for ", _obs_name, "added")
