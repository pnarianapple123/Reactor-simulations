from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


def flow_model(v,z,P):
        FA0=z[0]
        FB0=z[1]
        FC0=z[2]
        FI0=z[3]

        r_A=k*PA**(1/3)*PB**(2/3)
        
        dFAdV=rA
        dFBdV=rA*1/2
        dFCdV=-rA
        rI=0
        dFIdV=rI

        return [dFAdV,dFBdV,dFCdV,dFIdV]
#S0=[0.1362,0.068,0,0.256,1.0]

