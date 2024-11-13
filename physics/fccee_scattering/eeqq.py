r"""Functions for the process $e^+e^-\to q\bar q$."""
# Written by Eetu Loisa, 2024
from flavio.physics.zdecays.smeftew import gV_SM, gA_SM, _QN
import flavio.physics.zdecays.smeftew as smeftew
from flavio.physics.common import add_dict
from flavio.classes import Observable, Prediction
from flavio.physics import ckm as ckm_flavio
import wilson
from numpy import pi, sqrt
import numpy as np
from flavio.config import config


def F_qqll_SM(q, Xq, l, Xl, s, wc_obj, par):
    r"""SM $Z$ and $\gamma$ contribution to the $\bar q q\to \ell^+\ell^-$
    amplitude."""

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
    Qq = _QN[q]['Q']
    gVq = g_cw * gV_SM(q, par)
    gVl = g_cw * gV_SM(l, par)
    gAq = g_cw * gA_SM(q, par)
    gAl = g_cw * gA_SM(l, par)
    # Add SMEFT corrections to couplings
    gVq += g_cw * smeftew.d_gV(q, q, par, wcxf)
    gVl += g_cw * smeftew.d_gV(l, l, par, wcxf)
    gAq += g_cw * smeftew.d_gA(q, q, par, wcxf)
    gAl += g_cw * smeftew.d_gA(l, l, par, wcxf)

    if Xq == 'L':
        gZq = gVq + gAq
    elif Xq == 'R':
        gZq = gVq - gAq
    if Xl == 'L':
        gZl = gVl + gAl
    elif Xl == 'R':
        gZl = gVl - gAl
    # SM contribution
    return g2*s2w * Qq * Ql / s + gZq * gZl / (s - mZ**2 + 1j * mZ * GammaZ)


def wceff_qqll_sm(wc_obj, par, s):
    E = np.einsum('ij,kl->ijkl', np.eye(3), np.eye(3))
    wc = {}
    for Xl in 'LR':
        for Xq in 'LR':
            for q in 'ud':
                wc['CV{}{}_e{}'.format(Xl, Xq, q)] = F_qqll_SM(q, Xq, 'e', Xl, s, wc_obj, par) * E
                #print('Added CV{}{}_e{} to wc'.format(Xl, Xq, q) + ' with value' + str(wc['CV{}{}_e{}'.format(Xl, Xq, q)]) + 'at s =' + str(s))
    wc['CSRL_ed'] = wc['CSRR_ed'] = wc['CSLR_ed'] = wc['CSLL_ed'] = np.zeros((3, 3, 3, 3)) 
    wc['CSRL_eu'] = wc['CSRR_eu'] = wc['CSLR_eu'] = wc['CSLL_eu'] =  np.zeros((3, 3, 3, 3))
    wc['CTRL_ed'] = wc['CTRR_ed'] = wc['CTLR_ed'] = wc['CTLL_ed'] = np.zeros((3, 3, 3, 3))
    wc['CTRL_ed'] = wc['CTRR_eu'] = wc['CTLR_eu'] = wc['CTLL_eu'] = np.zeros((3, 3, 3, 3))
    return wc


def wceff_qqll_np(wc_obj, par, scale):
    r"""Returns Wilson coefficients of the effective Lagrangian
    $$\mathcal{L} = \sum_{q=u,d} C_{eq}^{\Gamma_1\Gamma_2}
    (\bar e_i\Gamma_1 e_j)(\bar q_k\Gamma_2 q_l)$$
    as a dictionary of arrays with shape (3, 3, 3, 3) corresponding to ijkl.
    """
    # get the dictionary
    wcxf_dict = wc_obj.get_wcxf(sector='dB=dL=0', scale=scale, par=par,
                                eft='SMEFT', basis='Warsaw up').dict
    # go to redundant basis, C has all the smeft coefficients
    C = wilson.util.smeftutil.wcxf2arrays_symmetrized(wcxf_dict)

    wc = {}
    ckm = ckm_flavio.get_ckm(par)
    # match to notation used in the observable
    # Note that the approach here is slightly different from flavio.physics.dileptons.ppll.py in that we add hermitian conjugates explicitly into the wc dictionary
    # CVXX_ex (8 entries):
    wc['CVLL_eu'] = C['lq1'] - C['lq3']
    wc['CVLL_ed'] = np.einsum('ijmn,mk,nl->ijkl',C['lq1'] + C['lq3'],np.conjugate(ckm),ckm)
    wc['CVRR_eu'] = C['eu']
    wc['CVRR_ed'] = C['ed']
    wc['CVLR_eu'] = C['lu']
    wc['CVLR_ed'] = C['ld']
    wc['CVRL_eu'] = np.einsum('klij->ijkl', C['qe'])
    wc['CVRL_ed'] = np.einsum('mnij,mk,nl->ijkl', C['qe'],np.conjugate(ckm),ckm)
    # CSXX_ex (8 entries):
    wc['CSRL_ed'] = np.einsum('ijkm,ml->ijkl',C['ledq'],ckm)
    wc['CSLR_ed'] = np.conjugate(np.transpose(np.einsum('ijkm,ml->ijkl',C['ledq'],ckm), [1,0,3,2])) # Hermitian conjugate of the line above
    wc['CSRR_ed'] = wc['CSLL_ed'] = np.zeros((3,3,3,3))
    wc['CSRR_eu'] = -C['lequ1']
    wc['CSLL_eu'] = np.conjugate(np.transpose(-C['lequ1'],[1,0,3,2])) # Hermitian conjugate of the line above
    wc['CSLR_eu'] = wc['CSRL_eu'] = np.zeros((3,3,3,3))
    # CTXX_xx (8 entries)
    wc['CTRR_eu'] = -C['lequ3']
    wc['CTLL_eu'] = np.conjugate(np.transpose(-C['lequ3'] ,[1,0,3,2]))
    wc['CTLR_eu'] = wc['CTRL_eu'] = np.zeros((3,3,3,3))
    wc['CTRR_ed'] = wc['CTLL_ed'] = np.zeros((3,3,3,3))
    wc['CTLR_ed'] = wc['CTRL_ed'] = np.zeros((3,3,3,3))

    return wc


# translate quark name to LHAPDF flavour index
fermion_indices = {
    'c': ('u', 1),
    'b': ('d', 2),
}

# this function returns the piece coming out of the integration of the Dirac traces
def s_func(s, coefficient):
    return coefficient * s**3

# Define a dictionary of the non-zero Dirac traces when integrating over cos(theta) in the range [-0.9, 0.9] 
f_integrated_ctheta_09 = {}
f_integrated_ctheta_09['SSLL'] = f_integrated_ctheta_09['SSLR'] = f_integrated_ctheta_09['SSRL'] = f_integrated_ctheta_09['SSRR'] = lambda s: s_func(s, 9/40)
f_integrated_ctheta_09['VVLL'] = f_integrated_ctheta_09['VVLR'] = f_integrated_ctheta_09['VVRL'] = f_integrated_ctheta_09['VVRR'] = lambda s: s_func(s, 1143/4000)
f_integrated_ctheta_09['TTLR'] = f_integrated_ctheta_09['TTRL'] = lambda s: s_func(s, 243/250)

# Define a dictionary of the non-zero Dirac traces when integrating over cos(theta) in the range [0, 0.9] MINUS the range [-0.9, 0], to get A_FB
f_AFB_ctheta_09 = {}
f_AFB_ctheta_09['VVLL'] = f_AFB_ctheta_09['VVRR'] = lambda s: s_func(s, 81/400)
f_AFB_ctheta_09['VVRL'] = f_AFB_ctheta_09['VVLR'] = lambda s: s_func(s, -81/400)
f_AFB_ctheta_09['STLR'] = f_AFB_ctheta_09['STRL'] = lambda s: s_func(s, -81/200)

def sigma_llqq_tot(wc_eff, s, q):
    r"""Total cross section of e^+ e^- \to q1 q2

    Returns $\sigma$ in units of GeV$^{-2}$

    Parameters:
    - `s`: centre of mass energy in GeV$^2$ of the lepton pair
    - `l`: lepton flavour, should be 'e'
    - `q1`, `q2`: outgoing quark flavours
    - `wc_eff`: SM + SMEFT amplitude for the process as dictionary with entries corresponding to different Wilson coefficients
    """
    q, i = fermion_indices[q]

    # Calculate the cross-section, sigma
    sigma = 0
    for gamma in 'SVT':
        for gammap in 'SVT':
            for X in 'LR':
                for Y in 'LR':
                    key = '{}{}{}{}'.format(gamma, gammap, X, Y)
                    if key in f_integrated_ctheta_09:
                        #print('key =', key)
                        #print('Does this return an integer or an array or something else: ', wc_eff['C{}{}{}_e{}'.format(gamma,X,Y,q)])
                        sigma += 3 / (16 * pi * s**2) * f_integrated_ctheta_09['{}{}{}{}'.format(gamma, gammap, X, Y)](s) * wc_eff['C{}{}{}_e{}'.format(gamma,X,Y,q)][0,0,i,i] * np.conjugate(wc_eff['C{}{}{}_e{}'.format(gammap,X,Y,q)][0,0,i,i])
                        #print('Added', key, 'to total cross-section')
                        #print('sigma =', sigma, 'at this point')
    # Return total cross-section
    return sigma

def sigma_llqq_forward_minus_backward(wc_eff, s, q):
    r"""Cross section of e^+ e^- \to q1 q2 in the forward region minus the backward region

    Returns $\sigma_F - \sigma_B$ in units of GeV$^{-2}$

    Parameters:
    - `s`: centre of mass energy in GeV$^2$ of the lepton pair
    - `l`: lepton flavour, should be 'e'
    - `q1`, `q2`: outgoing quark flavours
    - `wc_eff`: SM + SMEFT amplitude for the process as dictionary with entries corresponding to different Wilson coefficients
    """

    q, i = fermion_indices[q]

    # Calculate the cross-section in the forward minus backward region, sigma_F - sigma_B
    sigma = 0
    for gamma in 'SVT':
        for gammap in 'SVT':
            for X in 'LR':
                for Y in 'LR':
                    key = '{}{}{}{}'.format(gamma, gammap, X, Y)
                    if key in f_AFB_ctheta_09:
                        #print('key =', key)
                        #print('Does this return an integer or an array or something else: ', wc_eff['C{}{}{}_e{}'.format(gamma,X,Y,q)])
                        sigma += 3 / (16 * pi * s**2) * f_AFB_ctheta_09['{}{}{}{}'.format(gamma, gammap, X, Y)](s) * wc_eff['C{}{}{}_e{}'.format(gamma,X,Y,q)][0,0,i,i] * np.conjugate(wc_eff['C{}{}{}_e{}'.format(gammap,X,Y,q)][0,0,i,i])
                        #print('Added', key, 'to total cross-section')
                        #print('sigma =', sigma, 'at this point')
    # Return total cross-section
    return sigma


def sigma_llqq_tot_obs(wc_obj, par, E, q):
    # This function takes the BSM Wilson coefficients as input, runs them down, converts to wc_eff and adds to the SM contribution, and finally returns the cross-section
    scale = E
    s = E**2
    wc_eff_np = wceff_qqll_np(wc_obj, par, scale)
    #print('wc_eff_np =', wc_eff_np)
    wc_eff_sm = wceff_qqll_sm(wc_obj, par, s)
    wc_eff = add_dict((wc_eff_sm, wc_eff_np))
    #print('wc_eff keys' + str(wc_eff.keys()))

    return np.real(sigma_llqq_tot(wc_eff, s, q))

def AFB_llqq_obs(wc_obj, par, E, q):
    # This function takes the BSM Wilson coefficients as input, runs them down, converts to wc_eff and adds to the SM contribution, and finally returns the forward-backward asymmetry
    scale = E
    s = E**2
    wc_eff_np = wceff_qqll_np(wc_obj, par, scale)
    #print('wc_eff_np =', wc_eff_np)
    wc_eff_sm = wceff_qqll_sm(wc_obj, par, s)
    wc_eff = add_dict((wc_eff_sm, wc_eff_np))
    #print('wc_eff keys' + str(wc_eff.keys()))
    A_FB = sigma_llqq_forward_minus_backward(wc_eff, s, q) / sigma_llqq_tot(wc_eff, s, q)

    return np.real(A_FB)

# Function to help generate the total cross-section observable for a given quark q
def generate_sigma_obs(q):
    _selection_efficiencies = {'c': 0.03, 'b': 0.15}

    def f(wc_obj, par, E):
        return _selection_efficiencies[q] * sigma_llqq_tot_obs(wc_obj, par, E, q)
    return f

# Function to help generate the total cross-section observable for a given quark q
def generate_afb_obs(q):
    def f(wc_obj, par, E):
        return AFB_llqq_obs(wc_obj, par, E, q)
    return f

# Observable and Prediction instances
_tex = ['c','b']

# Create the observables and predictions for the total cross-section
for q in _tex:
    _process_tex = r"e^+e^- \to \overline{" + q + r"}" + q 
    _process_taxonomy = r'Process :: $e^+e^-$ scattering :: $' + _process_tex + r"$" #TODO Figure out how this should be formatted
    _obs_name = "sigma(ee->{}{})(high_E)".format(q,q)
    _obs = Observable(_obs_name, arguments=['E'])
    _obs.set_description(r"Cross section of $" + _process_tex + r"$")
    _obs.tex = r"$\sigma(" + _process_tex + r")$"
    _obs.add_taxonomy(_process_taxonomy)
    
    Prediction(_obs_name, generate_sigma_obs(q))
    #print("Prediction for ", _obs_name, "added")

# Create the observables and predictions for the forward-backward asymmetry
for q in _tex:
    _process_tex = r"A_{FB}(e^+e^- \to \overline{" + q + r"}" + q + r")"
    _process_taxonomy = r'Process :: $e^+e^-$ scattering :: $' + _process_tex + r"$" #TODO Figure out how this should be formatted
    _obs_name = "AFB(ee->{}{})(high_E)".format(q,q)
    _obs = Observable(_obs_name, arguments=['E'])
    _obs.set_description(r"Forward-backward asymmetry above Z pole $" + _process_tex + r"$")
    _obs.tex = r"$" + _process_tex + r"$"
    _obs.add_taxonomy(_process_taxonomy)
    
    Prediction(_obs_name, generate_afb_obs(q))
    #print("Prediction for ", _obs_name, "added")
