r"""Functions for Higgs production.

Most of the numerical coefficients have been obtained with MadGraph_aMC@NLO v2.6.5
along with SMEFTsim v2 in the alpha scheme.
"""

import flavio

def ggF(C):
    r"""Higgs production from gluon fusion normalized to the SM"""
    # obtained from an analytical one-loop calculation
    flavio.citations.register("Falkowski:2019hvp")
    np = (+35.86 * C['phiG']
          +0.121 * (C['phiBox'] - C['phiD'] / 4.)
          +0.061 * (C['ll_1221'] / 2 - C['phil3_22'] - C['phil3_11'])
          -0.129 * C['uphi_33']
          +0.123 * C['uphi_22']
          +0.239 * C['dphi_33']
          +0.025 * C['dphi_22']
          )
    return 1 + 1e6 * np.real

def hw(C):
    r"""Higgs production associated with a $W$ normalized to the SM"""
    flavio.citations.register("Falkowski:2019hvp")
    np = (+0.891 * C['phiW']
          -0.187 * C['phiWB']
          -0.115 * C['phiD']
          +0.121 * C['phiBox']
          +0.173 * (C['ll_1221'] / 2 - C['phil3_22'] - C['phil3_11'])
          +1.85 * C['phiq3_11']
          +0.126 * C['phiq3_22']
          )
    return 1 + 1e6 * np.real

def hz(C):
    r"""Higgs production associated with a $Z$ normalized to the SM"""
    flavio.citations.register("Falkowski:2019hvp")
    np = (+0.098 * C['phiB']
          +0.721 * C['phiW']
          +0.217 * C['phiWB']
          -0.015 * C['phiD']
          +0.122 * C['phiBox']
          +0.152 * (C['ll_1221'] / 2 - C['phil3_22'] - C['phil3_11'])
          -0.187 * C['phiq1_11']
          +1.699 * C['phiq3_11']
          +0.456 * C['phiu_11']
          -0.148 * C['phid_11']
          +0.044 * C['phiq1_22']
          +0.16 * C['phiq3_22']
          +0.028 * C['phiu_22']
          -0.02 * C['phid_22']
          )
    return 1 + 1e6 * np.real

def hv(C):
    r"""Higgs production associated with a $W$ or $Z$ normalized to the SM"""
    # Wh xsec at 14 TeV in pb, https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#s_13_0_TeV
    xw_sm = 1.380
    # Zh xsec at 14 TeV in pb, https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#s_13_0_TeV
    xz_sm = 0.8696
    d_hw = hw(C) - 1
    d_hz = hz(C) - 1
    return (xw_sm * (1 + d_hw) + xz_sm * (1 + d_hz)) / (xw_sm + xz_sm)

def tth(C):
    r"""Higgs production associated with a top pair normalized to the SM"""
    flavio.citations.register("Falkowski:2019hvp")
    np = (-0.030 * C['phiD']
          +0.118 * C['phiBox']
          -0.853 * C['uG_33']
          +0.146 * C['G']
          +0.557 * C['phiG']
          +0.060 * (C['ll_1221'] / 2 - C['phil3_22'] - C['phil3_11'])
          -0.119 * C['uphi_33'])
    return 1 + 1e6 * np.real

def vv_h(C):
    r"""Higgs production from vector boson fusion normalized to the SM"""
    flavio.citations.register("Falkowski:2019hvp")
    np = (-0.002 * C['phiB']
          -0.088 * C['phiW']
          -0.319 * C['phiWB']
          -0.168 * C['phiD']
          +0.121 * C['phiBox']
          +0.277 * (C['ll_1221'] / 2 - C['phil3_22'] - C['phil3_11'])
          +0.014 * C['phiq1_11']
          -0.384 * C['phiq3_11']
          -0.027 * C['phiu_11']
          +0.008 * C['phid_11']
          -0.004 * C['phiq1_22']
          -0.075 * C['phiq3_22']
          -0.004 * C['phiu_22']
          +0.002 * C['phid_22']
          )
    return 1 + 1e6 * np.real

###############################################
###       FCCee production processes        ###
### Computed with MadGraph_aMC and smeftsim ###
###         by Eetu Loisa, 2024             ###
###############################################

def hz_fccee240(C):
    r"""Higgs production associated with a Z boson normalised to the SM at FCC-ee at 240 GeV"""
    np = (+0.121 * C['phiBox'] 
          -0.006 * C['phiD'] 
          +0.529 * C['phiW'] 
          +0.127 * C['phiB'] 
          +0.236 * C['phiWB'] 
          +0.746 * C['phil3_11'] 
          -0.131 * C['phil3_22'] 
          +0.878 * C['phil1_11']  
          -0.769 * C['phie_11'] 
          +0.132 * C['ll_1221'] 
          ) 
    return 1 + 1e6 * np.real

def hz_fccee365(C):
    r"""Higgs production associated with a Z boson normalised to the SM at FCC-ee at 365 GeV"""
    np = (+0.121 * C['phiBox'] 
          -0.005 * C['phiD'] 
          +0.773 * C['phiW'] 
          +0.182 * C['phiB'] 
          +0.353 * C['phiWB'] 
          +1.782 * C['phil3_11'] 
          -0.131 * C['phil3_22'] 
          +1.912 * C['phil1_11'] 
          -1.677 * C['phie_11'] 
          +0.133 * C['ll_1221'] 
          )
    return 1 + 1e6 * np.real

def hnunu_fccee240(C):
    r"""Higgs production associated with a nu nubar pair, regardless of neutrino flavour, normalised to the SM at FCC-ee at 240 GeV. The cross-section includes the background process e+e- -> Z(->nu nubar)"""
    np = (+0.121 * C['phiBox'] 
          -0.004 * C['phiD'] 
          +0.456 * C['phiW'] 
          +0.110 * C['phiB'] 
          +0.233 * C['phiWB'] 
          +0.001 * C['eW_11'] 
          +0.001 * C['eB_11'] 
          +0.657 * C['phil3_11'] 
          -0.099 * C['phil3_22'] 
          +0.024 * C['phil3_33'] 
          +0.742 * C['phil1_11'] 
          -0.032 * C['phil1_22'] 
          -0.031 * C['phil1_33'] 
          -0.666 * C['phie_11'] 
          +0.005 * C['phie_22'] 
          +0.004 * C['phie_33'] 
          +0.123 * C['ll_1221'] 
          )
    return 1 + 1e6 * np.real

def hnunu_fccee365(C):
    r"""Higgs production associated with a nu nubar pair, regardless of neutrino flavour, normalised to the SM at FCC-ee at 365 GeV. The cross-section includes the background process e+e- -> Z(->nu nubar)"""
    np = (+0.122 * C['phiBox'] 
          -0.047 * C['phiD'] 
          +0.291 * C['phiW'] 
          +0.097 * C['phiB'] 
          +0.114 * C['phiWB'] 
          +0.001 * C['eW_11'] 
          +0.002 * C['eB_11'] 
          +0.531 * C['phil3_11'] 
          -0.154 * C['phil3_22'] 
          +0.013 * C['phil3_33'] 
          +0.805 * C['phil1_11'] 
          -0.017 * C['phil1_22'] 
          -0.016 * C['phil1_33'] 
          -0.762 * C['phie_11'] 
          +0.003 * C['phie_22'] 
          +0.003 * C['phie_33'] 
          +0.168 * C['ll_1221'] 
          )
    return 1 + 1e6 * np.real
