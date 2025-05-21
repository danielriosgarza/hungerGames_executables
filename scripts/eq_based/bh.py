#!/usr/bin/env python3
"""
bh_system_simulation_updated.py

Simulation of Bh sub-populations and metabolites with pH-dependent kinetics.
Produces four plots: sub-populations, metabolites, active cells, and pH.
Also writes results to CSV.
"""
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import gammaln

# Parameters
PARAMS: Dict[str, float] = {
    "xa_mumax": 0.1927934811,
    "xb_mumax": 0.9781040980,
    "xa_pHopt": 6.8180188470,
    "xb_pHopt": 6.5395960604,
    "xa_pHalpha": 46.2739563004,
    "xb_pHalpha": 62.5189011214,
    "xa_k_s1": 0.2876080490,
    "xa_k_s2": 2.3916509762,
    "xb_k_s3": 0.1039681099,
    "xb_k_s4": 1.2700484896,
    "xb_k_s2": 0.5,
    "xa_g_s1": 1.2417775764,
    "xa_g_s2": 1.0,
    "xa_g_s6_s1": 1.0e-10,
    "xa_g_s5_s1": 2.9999999999,
    "xa_g_s6_s2": 1.8699494242,
    "xb_g_s3": 3.7638796806,
    "xb_g_s4": 0.5,
    "xb_g_s2": 10.0,
    "xb_g_s6_s3_s4": 2.9999999999,
    "xb_g_s6_s2": 4.9999999999,
    "z1_r": 0.0250945181,
    "z2_r": 1.5,
    "z3_r": 1.0e-5,
    "z4_r": 0.2960374542,
    "z5_r": 0.0355569688,
    "z1_l_s1": 4.57e-8,
    "z2_l_s1": 0.21,
    "z4_l_s3_s4": 1.293346e-4,
    "z1_h_s1": 1.0,
    "z2_h_s1": 50.0,
    "z4_h_s3_s4": 2.9055456987,
}

# State indices
IDX = {"xa":0, "xb":1, "xc":2, "xd":3,
       "s1":4, "s2":5, "s3":6, "s4":7, "s5":8, "s6":9,
       "s7":10, "s8":11, "s9":12, "s10":13}

# Helper functions
def _gamma_pdf(x, alpha, beta):
    return 0.0 if x <= 0 else np.exp(alpha*np.log(beta) - gammaln(alpha) + (alpha-1)*np.log(x) - beta*x)

def xi(pH, pHopt, alpha):
    beta = (alpha-1)/pHopt
    return _gamma_pdf(pH, alpha, beta) / _gamma_pdf(pHopt, alpha, beta)

def calc_pH(s6, s8, s10):
    return -0.0491*s6 - 0.0972*s8 - 0.0436*s10 + 6.6369

def Z1(s1, p): return (p["z1_l_s1"]**p["z1_h_s1"]) / (s1**p["z1_h_s1"] + p["z1_l_s1"]**p["z1_h_s1"]) * p["z1_r"]

def Z2(s1, p): return (s1**p["z2_h_s1"]) / (s1**p["z2_h_s1"] + p["z2_l_s1"]**p["z2_h_s1"]) * p["z2_r"]

def Z3(p): return p["z3_r"]

def Z4(s3, s4, p):
    smin = min(s3, s4)
    return (p["z4_l_s3_s4"]**p["z4_h_s3_s4"]) / (smin**p["z4_h_s3_s4"] + p["z4_l_s3_s4"]**p["z4_h_s3_s4"]) * p["z4_r"]

# RHS function
def rhs(t, y, p):
    xa, xb, xc, xd = y[0:4]
    s1, s2, s3, s4, s5, s6, _, s8, _, _ = y[4:14]
    pH = calc_pH(s6, s8, y[13])  # butyrate at index 13
    xi_a = xi(pH, p["xa_pHopt"], p["xa_pHalpha"])
    xi_b = xi(pH, p["xb_pHopt"], p["xb_pHalpha"])
    z1, z2, z3, z4 = Z1(s1, p), Z2(s1, p), Z3(p), Z4(s3, s4, p)
    mu_a = s1/(s1 + p["xa_k_s1"]) + s2/(s2 + p["xa_k_s2"])
    mu_b = (s3/(s3 + p["xb_k_s3"]))*(s4/(s4 + p["xb_k_s4"])) + s2/(s2 + p["xb_k_s2"])
    dy = np.zeros_like(y)

    # Populations
    dy[IDX["xa"]] = xa*(xi_a*p["xa_mumax"]*mu_a - (z1 + z3)) + xb*z2
    dy[IDX["xb"]] = xb*(xi_b*p["xb_mumax"]*mu_b - (z2 + z4)) + xa*z1
    dy[IDX["xc"]] = xa*z3 + xb*z4 - xc*p["z5_r"]
    dy[IDX["xd"]] = xc*p["z5_r"]

    # Metabolites
    dy[IDX["s1"]] = -p["xa_g_s1"]*(s1/(s1 + p["xa_k_s1"]))*xi_a*p["xa_mumax"]*xa
    dy[IDX["s2"]] = -(
        p["xa_g_s2"]*(s2/(s2 + p["xa_k_s2"]))*xi_a*p["xa_mumax"]*xa +
        p["xb_g_s2"]*(s2/(s2 + p["xb_k_s2"]))*xi_b*p["xb_mumax"]*xb
    )
    dy[IDX["s3"]] = -p["xb_g_s3"]*(s3/(s3 + p["xb_k_s3"]))*(s4/(s4 + p["xb_k_s4"]))*xi_b*p["xb_mumax"]*xb
    dy[IDX["s4"]] = -p["xb_g_s4"]*(s3/(s3 + p["xb_k_s3"]))*(s4/(s4 + p["xb_k_s4"]))*xi_b*p["xb_mumax"]*xb
    dy[IDX["s5"]] = p["xa_g_s5_s1"]*(s1/(s1 + p["xa_k_s1"]))*xi_a*p["xa_mumax"]*xa
    dy[IDX["s6"]] = (
        (p["xa_g_s6_s1"]*(s1/(s1 + p["xa_k_s1"])) + p["xa_g_s6_s2"]*(s2/(s2 + p["xa_k_s2"]))) * xi_a*p["xa_mumax"]*xa +
        (p["xb_g_s6_s3_s4"]*(s3/(s3 + p["xb_k_s3"]))*(s4/(s4 + p["xb_k_s4"])) + p["xb_g_s6_s2"]*(s2/(s2 + p["xb_k_s2"]))) * xi_b*p["xb_mumax"]*xb
    )
    return dy

# Simulation function
def simulate(t_end=120.0, num_points=1000):
    y0 = np.array([
        0.0067733333, 0, 0, 0,
        0.684513088333333, 8.1368174275, 8.917955468, 1.0,
        0.272274459166667, 1.99119196216667, 1.0, 0.444971364333333,
        0.699764498666667, 0.003717077
    ])
    t_eval = np.linspace(0, t_end, num_points)
    sol = solve_ivp(
        rhs, (0, t_end), y0, args=(PARAMS,),
        t_eval=t_eval, method='DOP853', atol=1e-12, rtol=1e-12
    )
    df = pd.DataFrame(sol.y.T, columns=list(IDX.keys()), index=sol.t)
    df[df < 1e-8] = 0
    df = df.round(8)
    df['pH'] = calc_pH(df['s6'], df['s8'], df['s10']).clip(3, 10)
    return sol.t, df

# Main script
if __name__ == '__main__':
    t_end = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
    num_pts = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    t, df = simulate(t_end, num_pts)
    plt.figure(); plt.plot(t, df[['xa','xb','xc','xd']]); plt.legend(['x_a','x_b','x_c','x_d']); plt.xlabel('Time [h]'); plt.ylabel('10^5 cells/µL'); plt.title('Bh Sub-populations'); plt.grid(True)
    plt.figure(); plt.plot(t, df[['s1','s2','s3','s4','s5','s6']]); plt.legend(['Trehalose','Pyruvate','Glucose','Glutamate','Lactate','Acetate']); plt.xlabel('Time [h]'); plt.ylabel('mM'); plt.title('Extracellular Metabolites'); plt.grid(True)
    plt.figure(); plt.plot(t, df['xa']+df['xb']); plt.xlabel('Time [h]'); plt.ylabel('10^5 cells/µL'); plt.title('Active Cells'); plt.grid(True)
    plt.figure(); plt.plot(t, df['pH'], color='black'); plt.xlabel('Time [h]'); plt.ylabel('pH'); plt.title('Extracellular pH'); plt.grid(True)
    plt.tight_layout(); plt.show(); df.to_csv(Path('bh_simulation_results.csv'), index_label='time')
