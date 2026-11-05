#ODEs (numerical)
from matplotlib.pylab import solve
from matplotlib.pylab import solve
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def hiv_model(time,z):  #
    h=z[0]
    i=z[1]
    v=z[2]

    kr1=1e5
    kr2=0.1
    kr3=2e-7
    kr4=0.5
    kr5=5
    kr6=100

    dh_dt=kr1-kr2*h-kr3*h*v
    di_dt=kr3*h*v-kr4*i
    dv_dt=-kr3*h*v-kr5*v+kr6*i

    return [dh_dt,di_dt,dv_dt]

time_span = (0,15)

z_0=[1e6,0,100]

equation=solve(hiv_model,time_span,z_0)



