# ttnum11
(9*(419*e4 + 19*gR**4 + 38*e2*gR2)*s)/190. + (gR2*(-2*e2*GammaZ2*(mZ**4 + GammaZ2*mZ2 - s**2) + gR2*(GammaZ2 + mZ2)*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s)))*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))))/(GammaZ*mZ*(GammaZ2 + mZ2)) - (2*gR2*(-9*GammaZ*mZ*s - 10*(-(GammaZ2*mZ2) + (mZ2 + s)**2)*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) + 20*GammaZ*mZ*(mZ2 + s)*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))))*np.imag(Cee))/5. + (1143*s**3*np.imag(Cee)**2)/1000. - (2*(gR2*(gR2*(GammaZ2 + mZ2)*(mZ2 + s) + e2*((mZ2 + s)**2 + GammaZ2*(mZ2 + 2*s)))*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + e2*s*(e2*(GammaZ2 + mZ2) - gR2*s)*np.log(19)))/(GammaZ2 + mZ2) + ((9*s*(2*gR2*mZ2 + 3*(e2 + gR2)*s) - 40*GammaZ*gR2*mZ*(mZ2 + s)*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) - 20*gR2*(mZ*(GammaZ + mZ) + s)*(-(GammaZ*mZ) + mZ2 + s)*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) - 20*e2*s**2*np.log(19))*np.real(Cee))/5. + (1143*s**3*np.real(Cee)**2)/1000.

#ttnum22

(9*(419*e4 + 19*geX4 + 38*e2*geX2)*s)/190. + (geX2*(-2*e2*GammaZ2*(mZ**4 + GammaZ2*mZ2 - s**2) + geX2*(GammaZ2 + mZ2)*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s)))*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))))/(GammaZ*mZ*(GammaZ2 + mZ2)) - (2*geX2*(-9*GammaZ*mZ*s - 10*(-(GammaZ2*mZ2) + (mZ2 + s)**2)*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) + 20*GammaZ*mZ*(mZ2 + s)*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))))*np.imag(CXY))/5. + (1143*s**3*np.imag(CXY)**2)/1000. - (2*(geX2*(geX2*(GammaZ2 + mZ2)*(mZ2 + s) + e2*((mZ2 + s)**2 + GammaZ2*(mZ2 + 2*s)))*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + e2*s*(e2*(GammaZ2 + mZ2) - geX2*s)*np.log(19)))/(GammaZ2 + mZ2) + ((9*s*(2*geX2*mZ2 + 3*(e2 + geX2)*s) - 40*GammaZ*geX2*mZ*(mZ2 + s)*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) - 20*geX2*(mZ*(GammaZ + mZ) + s)*(-(GammaZ*mZ) + mZ2 + s)*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) - 20*e2*s**2*np.log(19))*np.real(CXY))/5. + (1143*s**3*np.real(CXY)**2)/1000.

#ttnum12 = ttnum21
(360*e4*s)/19. + (geLR*(2*e2*GammaZ2 + geLR*(GammaZ2 + mZ2))*s**2*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))))/(GammaZ*mZ*(GammaZ2 + mZ2)) + 2*geLR*s**2*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s)))*np.imag(CXY) + (9*s**3*np.imag(CXY)**2)/10. + (e2*geLR*s**2*(-2*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + np.log(361)))/(GammaZ2 + mZ2) - 2*s**2*(geLR*np.arctanh((180*s*(2*mZ2 + s))/(400*mZ**4 + 181*s**2 + 400*mZ2*(GammaZ2 + s))) + e2*np.log(19))*np.real(CXY) + (9*s**3*np.real(CXY)**2)/10.


# edits:
# - Replace Arctan with np.arctan
# - replace gL*gR with geLR
# - gL with geX, and the same for powers of gL


#stnum11
2*np.real((((27*s*(3*e2*s + geX2*(2*mZ*((-1j)*GammaZ + mZ) + 3*s)))/10. - (6j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) + (3429*s**3*np.conjugate(CXY))/1000. - 6*e2*s**2*np.log(19) + 3*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.log((GammaZ2*mZ2 + (mZ2 + s/20.)**2)/(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2)))*F_eeee_s_channel(wc_4e, X, X, s))/6.)

#stnum22
2*np.real((((27*s*(3*e2*s + gL2*(2*mZ*((-1j)*GammaZ + mZ) + 3*s)))/10. - (0,6)*gL2*((-1j)*GammaZ*mZ + mZ2 + s)**2*ArcTan((360*GammaZ*mZ*s)/(400*mZ**4 + 19*s**2 + 400*mZ2*(GammaZ2 + s))) + (3429*s**3*Conjugate(Cll))/1000. - 6*e2*s**2*Log(19) + 3*gL2*((-1j)*GammaZ*mZ + mZ2 + s)**2*Log((GammaZ2*mZ2 + (mZ2 + s/20.)**2)/(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2)))*Ss(s,2,2))/6.)
