from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

S0=[0.1362,0.068,0,0.256,1.0]
FA0=S0[0]
FB0=S0[1]
FC0=S0[2]
FI0=S0[3]
Pr=S0[4]
T0=260+275.5 # K (Isothermal)
P0=10e5 # Pa (Isobaric)
Ftot0=FA0+FB0+FC0+FI0
def ethylene_model(W,z,T,P):
        FA=z[0]
        FB=z[1]
        FC=z[2]
        FI=z[3]
        Pr=z[4]

        k=392 # mol(Pa Kg_cat s)^-1
        P=Pr*P0
        R=8.314 # J/mol-K
        Ftot=FA+FB+FC+FI

        Ctot0=P0/(R*T)
        CA=Ctot0*(FA/Ftot)*(P/P0)
        PA=CA*R*T
        CB=Ctot0*(FB/Ftot)*(P/P0)
        PB=CB*R*T
        rA=k*PA**(1/3)*PB**(2/3)
        
        dFAdV=-rA
        dFBdV=-rA*1/2
        dFCdV=rA
        rI=0
        dFIdV=rI
        a=0.0367 # Kg^-1
        dPrdW=-(a/2*Pr)*(Ftot/Ftot0)
        return [dFAdV,dFBdV,dFCdV,dFIdV,dPrdW]

#m_dot=Ftot0
#density=1925 # kg/m^3
#visco=1.825e-5 # Pa-s
#d=6.35e-3 # m
#Re=4*m_dot*d/(np.pi*density*visco)

# Since W is the cumulative weigt of the catalyst down the reator, 
# I will assume a very large mass value, to ensure that the conversion is reached within the range of the integration.
W_max=1 # kg 
W_span=(0,W_max)
Steps=1000
W_array=np.linspace(0,W_max,Steps)
solution = solve_ivp(ethylene_model,t_span=W_span,y0=S0,t_eval=W_array,args=(T0,P0))

W = solution.t

FA = solution.y[0,:]
FB = solution.y[1,:]    
FC = solution.y[2,:]
FI = solution.y[3,:]
Ftot=FA+FB+FC+FI
Pr = solution.y[4,:]
P=Pr*P0

X=1-FA/FA0

yA=FA/Ftot
yB=FB/Ftot
yC=FC/Ftot
yI=FI/Ftot

x_target=0.6
x_diff=np.abs(X-x_target)
target_index=np.argmin(x_diff)
W_actual=W[target_index]
P_final=P[target_index]
F_A_final=FA[target_index]
F_B_final=FB[target_index]
F_C_final=FC[target_index]
F_I_final=FI[target_index]
F_tot_final=F_A_final+F_B_final+F_C_final+F_I_final
yi_A_final=F_A_final/F_tot_final
yi_B_final=F_B_final/F_tot_final
yi_C_final=F_C_final/F_tot_final
yi_I_final=F_I_final/F_tot_final

print('At the conversion  X={0:7.6f}, the properties are:'.format(x_target))
print('      - W = {0:7.6f}'.format(W_actual))
print('      - P = {0:7.6f}'.format(P_final))
print('      - F_A = {0:7.6f}'.format(F_A_final))
print('      - F_B = {0:7.6f}'.format(F_B_final))
print('      - F_C = {0:7.6f}'.format(F_C_final))
print('      - F_I = {0:7.6f}'.format(F_I_final))
print('      - F_tot = {0:7.6f}'.format(F_tot_final))
print('      - y_A = {0:7.6f}'.format(yi_A_final))
print('      - y_B = {0:7.6f}'.format(yi_B_final))
print('      - y_C = {0:7.6f}'.format(yi_C_final))
print('      - y_I = {0:7.6f}'.format(yi_I_final))

# Create a plot object. The size is in inches
fig=plt.figure(figsize=(10,4))

# Create two subplots in the figure. 
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

plt.subplots_adjust(wspace=0.3,hspace=0.7) # This adjust the spacing between plots

# Plot the molar flow rates
ax1.plot(W,Ftot, '-k',label='Ftot')
ax1.plot(W,FA, '--b',label='FA')
ax1.plot(W,FB, '--r',label='FB')
ax1.plot(W,FC, '--g',label='FC')
ax1.plot(W,FI, '--m',label='FI')
ax1.set_xlabel('catalyst weight (kg)')
ax1.set_ylabel('molar flow rate')
ax1.set_xlim([0,W_max])
ax1.set_ylim([0,0.6])
ax1.legend()

# Plot conversion
ax2.plot(W,X, '-k',label='X')
ax2.set_xlabel('catalyst weight (kg)')
ax2.set_ylabel('conversion')
ax2.set_xlim([0,W_max])
ax2.set_ylim([0,1])
ax2.legend()

# Plot pressure
ax3.plot(W,P, '-k', label='P')
ax3.set_xlabel('catalyst weight (kg)')
ax3.set_ylabel('pressure')
ax3.set_xlim([0,W_max])
ax3.set_ylim([0,10e5])
ax3.legend()

# Plot mole fractions
ax4.plot(W,yA, '-b', label='yA')
ax4.plot(W,yB, '-r', label='yB')
ax4.plot(W,yC, '-g', label='yC')
ax4.plot(W,yI, '-m', label='yI')
ax4.set_xlabel('catalyst weight (kg)')
ax4.set_ylabel('mole fraction')
ax4.set_xlim([0,W_max])
ax4.set_ylim([0,1])
ax4.legend()

# Save figure 
plt.savefig("ethylene.pdf")
# Display figure
plt.show()