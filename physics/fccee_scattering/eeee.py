r"""Functions for the process $e^+e^- \to e^+e^-$."""
# Written by Eetu Loisa, 2024
# partially based on flavio.physics.dileptons by Greljo, Salko, Smolkovic and Stangl

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

def sigma_eeee_tot(wc_obj, par, E):
    s = E**2

    # get the Wilson coefficient dictionary, in a non-redundant basis
    wcxf = wc_obj.get_wcxf(sector='dB=dL=0', scale=E, par=par, eft='SMEFT', basis='Warsaw up')
    wcxf_dict = wc_obj.get_wcxf(sector='dB=dL=0', scale=E, par=par, eft='SMEFT', basis='Warsaw up').dict

    # Get the Wilson coefficients that contribute to Bhabha scattering, categorised by chirality
    wc_4e = {'RR': 2 * wcxf_dict.get('ee_1111',0), 'LL': 2 * wcxf_dict.get('ll_1111',0), 'LR': wcxf_dict.get('le_1111',0), 'RL': wcxf_dict.get('le_1111',0)}

    #SM parameters
    Ql = -1
    mZ = par['m_Z']
    mZ2 = mZ**2 # mZ squared and GammaZ squared are defined here because they appear so often. mZ3, mZ4, etc. could also be defined if needed for speed.
    GammaZ = 1 / par['tau_Z']
    GammaZ2= GammaZ**2  
    s2w = par['s2w']
    aEW = par['alpha_e']
    e2 = 4 * pi * aEW
    e4 = e2 * e2
    g2=e2/s2w 
    gp2=e2/(1-s2w)  # We do not run the couplings as even at 365 GeV they change by less than 1%
    g_cw = sqrt(g2+gp2) # g over cw
    s2w = gp2/(g2+gp2)
    par = par.copy()
    par['s2w'] = s2w
    # Electron vector and axial couplings in the SM
    gVe =  g_cw * gV_SM('e', par)
    gAe =  g_cw * gA_SM('e', par)
    # Add SMEFT contributions to the electron couplings
    gVe += g_cw * smeftew.d_gVl('e', 'e', par, wcxf)
    gAe += g_cw * smeftew.d_gAl('e', 'e', par, wcxf)


    geX_dict = {'L': gVe + gAe, 'R': gVe - gAe}
    # For faster evaluation, define the squares and product of the Z couplings
    geLR = geX_dict['L'] * geX_dict['R']

    def F_eeee_s_channel(wc_4e, Xl1, Xl2, s):
        r"""$Z$ and $\gamma$ contribution to the s-channel $e^+ e^- \to e^+ e^-$ amplitude, 
        including SMEFT contributions.
        - Xl1 and Xl2: in- and out-going electron chiralities, 'L' or 'R'
        - s: Mandelstam variable s = E^2"""
        # NOTE: Only the s-channel is implemented as its own function as it does not contain any t-dependence, and passes through to the total cross-section untouched.

        return g2*s2w * Ql * Ql / s + geX_dict[Xl1] * geX_dict[Xl2] / (s - mZ2 + 1j * mZ * GammaZ) + wc_4e[f'{Xl1}{Xl2}']

    # This is a dictionary containing the different chirality and s/t-channel contributions that make up the total Bhabha scattering cross-section
    contributions = {}

    # s-channel squared contributions
    for Xl1 in 'LR':
        for Xl2 in 'LR':
            contributions[f"S^2_{Xl1}{Xl2}"] = (1143 * s**3 * np.abs(F_eeee_s_channel(wc_4e, Xl1, Xl2, s))**2 )/4000

    # t-channel squared contributions, LL and RR
    for X in 'LR':
        geX = geX_dict[X]
        geX2 = geX**2 # Define powers of couplings to speed up evaluation
        geX4 = geX**4
        CXY = wc_4e[f'{X}{X}'] / 2
        contributions[f'T^2_{X}{X}'] = (9*(419*e4 + 19*geX4 + 38*e2*geX2)*s)/190. + (geX2*(-2*e2*GammaZ2*(mZ**4 + GammaZ2*mZ2 - s**2) + geX2*(GammaZ2 + mZ2)*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s)))*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))))/(GammaZ*mZ*(GammaZ2 + mZ2)) - (2*geX2*(-9*GammaZ*mZ*s - 10*(-(GammaZ2*mZ2) + (mZ2 + s)**2)*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) + 20*GammaZ*mZ*(mZ2 + s)*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))))*np.imag(CXY))/5. + (1143*s**3*np.imag(CXY)**2)/1000. - (2*(geX2*(geX2*(GammaZ2 + mZ2)*(mZ2 + s) + e2*((mZ2 + s)**2 + GammaZ2*(mZ2 + 2*s)))*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + e2*s*(e2*(GammaZ2 + mZ2) - geX2*s)*np.log(19)))/(GammaZ2 + mZ2) + ((9*s*(2*geX2*mZ2 + 3*(e2 + geX2)*s) - 40*GammaZ*geX2*mZ*(mZ2 + s)*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) - 20*geX2*(mZ*(GammaZ + mZ) + s)*(-(GammaZ*mZ) + mZ2 + s)*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) - 20*e2*s**2*np.log(19))*np.real(CXY))/5. + (1143*s**3*np.real(CXY)**2)/1000.

    # t-channel squared contributions, LR and RL
    CXY = wc_4e['LR']
    contributions['T^2_LR'] = contributions['T^2_RL'] = (360*e4*s)/19. + (geLR*(2*e2*GammaZ2 + geLR*(GammaZ2 + mZ2))*s**2*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))))/(GammaZ*mZ*(GammaZ2 + mZ2)) + 2*geLR*s**2*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s)))*np.imag(CXY) + (9*s**3*np.imag(CXY)**2)/10. + (e2*geLR*s**2*(-2*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + np.log(361)))/(GammaZ2 + mZ2) - 2*s**2*(geLR*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + e2*np.log(19))*np.real(CXY) + (9*s**3*np.real(CXY)**2)/10.

    # s- and t-channel interference contributions
    for X in 'LR':
        geX = geX_dict[X]
        geX2 = geX**2 # Define powers of couplings to speed up evaluation
        geX4 = geX**4
        CXY = wc_4e[f'{X}{X}'] / 2
        contributions[f'ST_{X}{X}'] = 2*np.real((((27*s*(3*e2*s + geX2*(2*mZ*((-1j)*GammaZ + mZ) + 3*s)))/10. - (6j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) + (3429*s**3*np.conjugate(CXY))/1000. - 6*e2*s**2*np.log(19) + 3*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.log((GammaZ2*mZ2 + (mZ2 + s/20.)**2)/(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2)))*F_eeee_s_channel(wc_4e, X, X, s))/6.)

    sigma = 0

    # Add all contributions together to get the total cross-section
    for key, value in contributions.items():
        sigma += 1 / (16 * pi * s**2) * value

    return sigma


def AFB_eeee(wc_obj, par, E):
    s = E**2

    # get the Wilson coefficient dictionary, in a non-redundant basis
    wcxf_dict = wc_obj.get_wcxf(sector='dB=dL=0', scale=E, par=par, eft='SMEFT', basis='Warsaw up').dict

    # Get the Wilson coefficients that contribute to Bhabha scattering, categorised by chirality
    wc_4e = {'RR': 2 * wcxf_dict.get('ee_1111',0), 'LL': 2 * wcxf_dict.get('ll_1111',0), 'LR': wcxf_dict.get('le_1111',0), 'RL': wcxf_dict.get('le_1111',0)}

    #SM parameters
    Ql = -1
    mZ = par['m_Z']
    mZ2 = mZ**2 # mZ squared and GammaZ squared are defined here because they appear so often. mZ3, mZ4, etc. could also be defined if needed for speed.
    GammaZ = 1 / par['tau_Z']
    GammaZ2= GammaZ**2  
    s2w = par['s2w']
    aEW = par['alpha_e']
    e2 = 4 * pi * aEW
    e4 = e2 * e2
    g2=e2/s2w 
    gp2=e2/(1-s2w)  # We do not run the couplings as even at 365 GeV they change by less than 1%
    g_cw = sqrt(g2+gp2) # g over cw
    s2w = gp2/(g2+gp2)
    par = par.copy()
    par['s2w'] = s2w
    # Electron vector and axial couplings
    gVe =  g_cw * gV_SM('e', par)
    gAe =  g_cw * gA_SM('e', par)

    geX_dict = {'L': gVe + gAe, 'R': gVe - gAe}
    # For faster evaluation, define the squares and product of the Z couplings
    geLR = geX_dict['L'] * geX_dict['R']

    def F_eeee_s_channel(wc_4e, Xl1, Xl2, s):
        r"""$Z$ and $\gamma$ contribution to the s-channel $e^+ e^- \to e^+ e^-$ amplitude, 
        including SMEFT contributions.
        - Xl1 and Xl2: in- and out-going electron chiralities, 'L' or 'R'
        - s: Mandelstam variable s = E^2"""
        # NOTE: Only the s-channel is implemented as its own function as it does not contain any t-dependence, and passes through to the total cross-section untouched.

        return g2*s2w * Ql * Ql / s + geX_dict[Xl1] * geX_dict[Xl2] / (s - mZ2 + 1j * mZ * GammaZ) + wc_4e[f'{Xl1}{Xl2}']

    # This is a dictionary containing the different chirality and s/t-channel contributions that make up the total Bhabha scattering cross-section
    contributions = {}

    # s-channel squared contributions
    for X in 'LR':
        contributions[f"S^2_{X}{X}"] = (81 * s**3 * np.abs(F_eeee_s_channel(wc_4e, X, X, s))**2 )/400

    contributions[f"S^2_LR"] = contributions[f"S^2_RL"] = - (81 * s**3 * np.abs(F_eeee_s_channel(wc_4e, 'L', 'R', s))**2 )/400

    # t-channel squared contributions, LL and RR
    for X in 'LR':
        geX = geX_dict[X]
        geX2 = geX**2 # Define powers of couplings to speed up evaluation
        geX4 = geX**4
        CXY = wc_4e[f'{X}{X}'] / 2
        contributions[f'T^2_{X}{X}'] = (81*s**3*np.imag(CXY)**2)/100. + (19*geX2*(-2*e2*GammaZ2*(mZ**4 + GammaZ2*mZ2 - s**2) + geX2*(GammaZ2 + mZ2)*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s)))*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) + GammaZ*mZ*(2*e2*s*(19*geX2*s*np.log(5.2631578947368425) + e2*(GammaZ2 + mZ2)*(162 + 19*np.log(1.9) - 19*np.log(10))) + 19*geX2*(geX2*(GammaZ2 + mZ2)*(mZ2 + s) + e2*(mZ**4 + s*(2*GammaZ2 + s) + mZ2*(GammaZ2 + 2*s)))*np.log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 19*(geX4*(GammaZ2 + mZ2)*(mZ2 + s) + e2*geX2*(mZ**4 + s*(2*GammaZ2 + s) + mZ2*(GammaZ2 + 2*s)))*np.log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))/(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.)))))/(19.*GammaZ*mZ*(GammaZ2 + mZ2)) + 4*geX2*np.imag(CXY)*((mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s))*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - GammaZ*mZ*(mZ2 + s)*np.log((10000*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s))**2)/((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))*(400*mZ**4 + 361*s**2 + 40*mZ2*(10*GammaZ2 + 19*s))))) + ((81*e2*s**2 + 81*geX2*s**2 - 800*GammaZ*geX2*mZ*(mZ2 + s)*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) + 400*e2*s**2*np.log(1.9) - 400*e2*s**2*np.log(10) + 200*geX2*mZ**4*np.log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 200*GammaZ2*geX2*mZ2*np.log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) + 400*geX2*mZ2*s*np.log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) + 200*geX2*s**2*np.log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 200*geX2*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s))*np.log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))/(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.))))*np.real(CXY))/100. + (81*s**3*np.real(CXY)**2)/100.


    # t-channel squared contributions, LR and RL
    CXY = wc_4e['LR']
    contributions['T^2_LR'] = contributions['T^2_RL'] = -0.05263157894736842*(s*(-324*e4*GammaZ**3*mZ - 324*e4*GammaZ*mZ**3 - 38*e2*GammaZ2*geLR*s*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - 19*GammaZ2*geLR**2*s*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - 19*geLR**2*mZ2*s*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - 38*GammaZ*geLR*mZ*(GammaZ2 + mZ2)*s*np.arctan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2)))*np.imag(CXY) + 38*e2*GammaZ*geLR*mZ*s*np.log(1.9) - 19*e2*GammaZ*geLR*mZ*s*np.log(100) + 19*e2*GammaZ*geLR*mZ*s*np.log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))**2/((mZ**4 + mZ2*(GammaZ2 + s/10.) + s**2/400.)*(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.)))) + 19*GammaZ*mZ*(GammaZ2 + mZ2)*s*(e2*np.log(27.700831024930746) - geLR*np.log(((400*GammaZ2*mZ2 + (20*mZ2 + s)**2)*(400*GammaZ2*mZ2 + (20*mZ2 + 19*s)**2))/(10000.*(4*GammaZ2*mZ2 + (2*mZ2 + s)**2)**2)))*np.real(CXY)))/(GammaZ*mZ*(GammaZ2 + mZ2))


    # s- and t-channel interference contributions
    for X in 'LR':
        geX = geX_dict[X]
        geX2 = geX**2 # Define powers of couplings to speed up evaluation
        geX4 = geX**4
        CXY = wc_4e[f'{X}{X}'] / 2
        contributions[f'ST_{X}{X}'] = 2*np.real(((162*(e2 + geX2)*s**2 - (800j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((18*GammaZ*mZ*s)/(40*mZ**4 + s**2 + mZ2*(40*GammaZ2 + 22*s))) + 324*s**3*np.conjugate(CXY) + 400*((2j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((18*GammaZ*mZ*s)/(40*mZ**4 + 19*s**2 + mZ2*(40*GammaZ2 + 58*s))) - 2*e2*s**2*np.log(5.2631578947368425) + geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*(np.log(GammaZ2*mZ2 + (mZ2 + s/20.)**2) - 2*np.log(GammaZ2*mZ2 + (mZ2 + s/2.)**2) + np.log(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2))))*F_eeee_s_channel(wc_4e,X,X,s))/800.)

    sigma_FB = 0 #Forward - backward cross-section

    for key, value in contributions.items():
        sigma_FB += 1 / (16 * pi * s**2) * value  

    sigma_tot = sigma_eeee_tot(wc_obj, par, E) # Get the total cross-section

    AFB = sigma_FB / sigma_tot # Calculate the forward-backward asymmetry

    return AFB


# Function to help generate the total cross-section observable for a given quark q
def generate_sigma_obs():
    _selection_efficiency = 0.98

    def f(wc_obj, par, E):
        return _selection_efficiency * sigma_eeee_tot(wc_obj, par, E)
    return f

# Observable and Prediction instances

# Create the observables and predictions for the total cross-section
_process_tex = r"e^+e^- \to e^+e^-"
_process_taxonomy = r'Process :: $e^+e^-$ scattering :: $' + _process_tex + r"$" #TODO Figure out how this should be formatted
_obs_name = "sigma(ee->ee)(high_E)"
_obs = Observable(_obs_name, arguments=['E'])
_obs.set_description(r"Cross section of $" + _process_tex + r"$")
_obs.tex = r"$\sigma(" + _process_tex + r")$"
_obs.add_taxonomy(_process_taxonomy)

Prediction(_obs_name, generate_sigma_obs())

# Create the observables and predictions for the forward-backward asymmetry
_process_tex = r"A_{FB}(e^+e^- \to e^+e^-)"
_process_taxonomy = r'Process :: $e^+e^-$ scattering :: $' + _process_tex + r"$" #TODO Figure out how this should be formatted
_obs_name = "AFB(ee->ee)(high_E)"
_obs = Observable(_obs_name, arguments=['E'])
_obs.set_description(r"Forward-backward asymmetry above Z pole $" + _process_tex + r"$")
_obs.tex = r"$" + _process_tex + r"$"
_obs.add_taxonomy(_process_taxonomy)

Prediction(_obs_name, AFB_eeee)
