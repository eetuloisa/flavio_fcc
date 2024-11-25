import re

# Script to automatically turn expressions in Fortran form into desirerable python form using regular expressions

fbTT11 = "(81*s**3*Im(Cee)**2)/100. + (19*gR2*(-2*e2*GammaZ2*(mZ**4 + GammaZ2*mZ2 - s**2) + gR2*(GammaZ2 + mZ2)*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s)))*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) + GammaZ*mZ*(2*e2*s*(19*gR2*s*Log(5.2631578947368425) + e2*(GammaZ2 + mZ2)*(162 + 19*Log(1.9) - 19*Log(10))) + 19*gR2*(gR2*(GammaZ2 + mZ2)*(mZ2 + s) + e2*(mZ**4 + s*(2*GammaZ2 + s) + mZ2*(GammaZ2 + 2*s)))*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 19*(gR4*(GammaZ2 + mZ2)*(mZ2 + s) + e2*gR2*(mZ**4 + s*(2*GammaZ2 + s) + mZ2*(GammaZ2 + 2*s)))*Log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))/(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.)))))/(19.*GammaZ*mZ*(GammaZ2 + mZ2)) + 4*gR2*Im(Cee)*((mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s))*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - GammaZ*mZ*(mZ2 + s)*Log((10000*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s))**2)/((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))*(400*mZ**4 + 361*s**2 + 40*mZ2*(10*GammaZ2 + 19*s))))) + ((81*e2*s**2 + 81*gR2*s**2 - 800*GammaZ*gR2*mZ*(mZ2 + s)*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) + 400*e2*s**2*Log(1.9) - 400*e2*s**2*Log(10) + 200*gR2*mZ**4*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 200*GammaZ2*gR2*mZ2*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) + 400*gR2*mZ2*s*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) + 200*gR2*s**2*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 200*gR2*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s))*Log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))/(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.))))*Re(Cee))/100. + (81*s**3*Re(Cee)**2)/100."

fbTT22 = "(81*s**3*Im(Cll)**2)/100. + (19*gL2*(-2*e2*GammaZ2*(mZ**4 + GammaZ2*mZ2 - s**2) + gL2*(GammaZ2 + mZ2)*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s)))*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) + GammaZ*mZ*(2*e2*s*(19*gL2*s*Log(5.2631578947368425) + e2*(GammaZ2 + mZ2)*(162 + 19*Log(1.9) - 19*Log(10))) + 19*gL2*(gL2*(GammaZ2 + mZ2)*(mZ2 + s) + e2*(mZ**4 + s*(2*GammaZ2 + s) + mZ2*(GammaZ2 + 2*s)))*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 19*(gL4*(GammaZ2 + mZ2)*(mZ2 + s) + e2*gL2*(mZ**4 + s*(2*GammaZ2 + s) + mZ2*(GammaZ2 + 2*s)))*Log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))/(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.)))))/(19.*GammaZ*mZ*(GammaZ2 + mZ2)) + 4*gL2*Im(Cll)*((mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s))*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - GammaZ*mZ*(mZ2 + s)*Log((10000*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s))**2)/((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))*(400*mZ**4 + 361*s**2 + 40*mZ2*(10*GammaZ2 + 19*s))))) + ((81*e2*s**2 + 81*gL2*s**2 - 800*GammaZ*gL2*mZ*(mZ2 + s)*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) + 400*e2*s**2*Log(1.9) - 400*e2*s**2*Log(10) + 200*gL2*mZ**4*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 200*GammaZ2*gL2*mZ2*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) + 400*gL2*mZ2*s*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) + 200*gL2*s**2*Log((400*mZ**4 + s**2 + 40*mZ2*(10*GammaZ2 + s))/(100.*(4*mZ**4 + s**2 + 4*mZ2*(GammaZ2 + s)))) - 200*gL2*(mZ**4 + s**2 + mZ2*(-GammaZ2 + 2*s))*Log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))/(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.))))*Re(Cll))/100. + (81*s**3*Re(Cll)**2)/100."

fbTT12 = "-0.05263157894736842*(s*(-324*e4*GammaZ**3*mZ - 324*e4*GammaZ*mZ**3 - 38*e2*GammaZ2*gL*gR*s*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - 19*GammaZ2*gL2*gR2*s*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - 19*gL2*gR2*mZ2*s*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2))) - 38*GammaZ*gL*gR*mZ*(GammaZ2 + mZ2)*s*ArcTan((324*GammaZ*mZ*s**2*(2*mZ2 + s))/(1600*GammaZ**4*mZ**4 + (2*mZ2 + s)**2*(20*mZ2 + s)*(20*mZ2 + 19*s) + 4*GammaZ2*mZ2*(800*mZ**4 + 800*mZ2*s + 281*s**2)))*Im(Cle) + 38*e2*GammaZ*gL*gR*mZ*s*Log(1.9) - 19*e2*GammaZ*gL*gR*mZ*s*Log(100) + 19*e2*GammaZ*gL*gR*mZ*s*Log((mZ**4 + s**2/4. + mZ2*(GammaZ2 + s))**2/((mZ**4 + mZ2*(GammaZ2 + s/10.) + s**2/400.)*(mZ**4 + (361*s**2)/400. + mZ2*(GammaZ2 + (19*s)/10.)))) + 19*GammaZ*mZ*(GammaZ2 + mZ2)*s*(e2*Log(27.700831024930746) - gL*gR*Log(((400*GammaZ2*mZ2 + (20*mZ2 + s)**2)*(400*GammaZ2*mZ2 + (20*mZ2 + 19*s)**2))/(10000.*(4*GammaZ2*mZ2 + (2*mZ2 + s)**2)**2)))*Re(Cle)))/(GammaZ*mZ*(GammaZ2 + mZ2))"

fbST11 = "((162*(e2 + gR2)*s**2 - (800j)*gR2*((-1j)*GammaZ*mZ + mZ2 + s)**2*ArcTan((18*GammaZ*mZ*s)/(40*mZ**4 + s**2 + mZ2*(40*GammaZ2 + 22*s))) + 324*s**3*Conjugate(Cee) + 400*((2j)*gR2*((-1j)*GammaZ*mZ + mZ2 + s)**2*ArcTan((18*GammaZ*mZ*s)/(40*mZ**4 + 19*s**2 + mZ2*(40*GammaZ2 + 58*s))) - 2*e2*s**2*Log(5.2631578947368425) + gR2*((-1j)*GammaZ*mZ + mZ2 + s)**2*(Log(GammaZ2*mZ2 + (mZ2 + s/20.)**2) - 2*Log(GammaZ2*mZ2 + (mZ2 + s/2.)**2) + Log(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2))))*Ss(s,1,1))/800."

fbST22 = "((162*(e2 + gL2)*s**2 - (800j)*gL2*((-1j)*GammaZ*mZ + mZ2 + s)**2*ArcTan((18*GammaZ*mZ*s)/(40*mZ**4 + s**2 + mZ2*(40*GammaZ2 + 22*s))) + 324*s**3*Conjugate(Cll) + 400*((2j)*gL2*((-1j)*GammaZ*mZ + mZ2 + s)**2*ArcTan((18*GammaZ*mZ*s)/(40*mZ**4 + 19*s**2 + mZ2*(40*GammaZ2 + 58*s))) - 2*e2*s**2*Log(5.2631578947368425) + gL2*((-1j)*GammaZ*mZ + mZ2 + s)**2*(Log(GammaZ2*mZ2 + (mZ2 + s/20.)**2) - 2*Log(GammaZ2*mZ2 + (mZ2 + s/2.)**2) + Log(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2))))*Ss(s,2,2))/800."

def apply_regex_rules(input_string, diag):
    # Define the regex rules as tuples of (pattern, replacement)
    if diag == True:
        rules = [
            (r"ArcTan\(", "np.arctan("),
            (r"Log\(", "np.log("),
            (r"ArcTanh\(", "np.arctanh("),
            (r"gL\*gR", "geLR"),
            (r"gR4", "geX4"),     # Order matters here: handle gR4 before gR
            (r"gR2", "geX2"),     # Handle gR2 before gR
            (r"gR", "geX"),
            (r"gL4", "geX4"),     # Order matters here: handle gL4 before gR
            (r"gL2", "geX2"),     # Handle gR2 before gR
            (r"gL", "geX"),
            (r"Re\(", "np.real("),
            (r"Im\(", "np.imag("),
            (r"Cll|Cee|Cle", "CXY")
        ]
    else:
        rules = [
            (r"ArcTan\(", "np.arctan("),
            (r"Log\(", "np.log("),
            (r"ArcTanh\(", "np.arctanh("),
            (r"gL\*gR", "geLR"),
            (r"gL2\*gR2", "geLR**2"),
            (r"Re\(", "np.real("),
            (r"Im\(", "np.imag("),
            (r"Cll|Cee|Cle", "CXY")
        ]
    
    # Apply each regex rule sequentially
    for pattern, replacement in rules:
        input_string = re.sub(pattern, replacement, input_string)
    
    return input_string

def save_to_file(output_string, filename="regex_output.txt"):
    with open(filename, "a") as file:
        file.write(output_string + "\n\n")
    print(f"Output saved to {filename}")

# Example usage
input_string = fbTT11
output_string = apply_regex_rules(input_string, True)
save_to_file(output_string, "afb_regex_output.txt")
input_string = fbTT22
output_string = apply_regex_rules(input_string, True)
save_to_file(output_string, "afb_regex_output.txt")
input_string = fbTT12
output_string = apply_regex_rules(input_string, False)
save_to_file(output_string, "afb_regex_output.txt")

# ST contributions
input_string = fbST11
output_string = apply_regex_rules(input_string, True)
save_to_file(output_string, "afb_regex_output.txt")
input_string = fbST22
output_string = apply_regex_rules(input_string, True)
save_to_file(output_string, "afb_regex_output.txt")

def test_if_equal(string1, string2):
    return string1 == string2

string1 = "((162*(e2 + geX2)*s**2 - (800j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((18*GammaZ*mZ*s)/(40*mZ**4 + s**2 + mZ2*(40*GammaZ2 + 22*s))) + 324*s**3*np.conjugate(CXY) + 400*((2j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((18*GammaZ*mZ*s)/(40*mZ**4 + 19*s**2 + mZ2*(40*GammaZ2 + 58*s))) - 2*e2*s**2*np.log(5.2631578947368425) + geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*(np.log(GammaZ2*mZ2 + (mZ2 + s/20.)**2) - 2*np.log(GammaZ2*mZ2 + (mZ2 + s/2.)**2) + np.log(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2))))*Ss(s,1,1))/800."

string2 = "((162*(e2 + geX2)*s**2 - (800j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((18*GammaZ*mZ*s)/(40*mZ**4 + s**2 + mZ2*(40*GammaZ2 + 22*s))) + 324*s**3*np.conjugate(CXY) + 400*((2j)*geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*np.arctan((18*GammaZ*mZ*s)/(40*mZ**4 + 19*s**2 + mZ2*(40*GammaZ2 + 58*s))) - 2*e2*s**2*np.log(5.2631578947368425) + geX2*((-1j)*GammaZ*mZ + mZ2 + s)**2*(np.log(GammaZ2*mZ2 + (mZ2 + s/20.)**2) - 2*np.log(GammaZ2*mZ2 + (mZ2 + s/2.)**2) + np.log(GammaZ2*mZ2 + (mZ2 + (19*s)/20.)**2))))*Ss(s,1,1))/800."

print("Strings are equal: ", test_if_equal(string1, string2))
