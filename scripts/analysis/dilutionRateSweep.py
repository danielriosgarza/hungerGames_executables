# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:11:52 2023

@author: drgarza
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})

sys.path.append(os.getcwd())
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'core'))
sys.path.append(os.path.join(Path(os.getcwd()).parents[1], 'files'))

from param_change import *
from readModelDB import get_database, createMetabolome, createBacteria
from mainClasses import *
from param_change import *

import cmasher as cmr

from tqdm import tqdm



import plotly.io as pio
pio.renderers.default='browser'

from matplotlib.colors import ListedColormap




def create_heatmap(data, row_labels, x_positions, x_values, xlabel, ylabel, title, fileName = None):
    fig, ax = plt.subplots()
    heatmap = ax.imshow(data, aspect='auto')
    ax.grid(False)

    # Set row labels (y-axis)
    ax.set_yticks(np.arange(data.shape[0]))
    
    yticklabels = ax.set_yticklabels(row_labels)
    yticklabels[-1].set_fontstyle('italic')
    yticklabels[-2].set_fontstyle('italic')
    yticklabels[-3].set_fontstyle('italic')

    # Set specific x-ticks and labels (x-axis)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_values)

    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # Add grid
    ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, data.shape[0], 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    # Add colorbar
    cbar = ax.figure.colorbar(heatmap, ax=ax)
    plt.tight_layout()
    
    if fileName is not None:
        plt.savefig(fileName, transparent=True, dpi=600)

    plt.show()


def mockpHfunc(metObj, pH=7.0):
    def pHfunc(metObj):
        return pH
    return pHfunc

def makeSimulation(pHControl = None, 
                   dilutionFactor = 0.5, 
                   xa = 0.03,
                   xb = 0.001,
                   xe = 0.03,
                   xf = 0.001,
                   xi = 0.03,
                   xj = 0.001):
    
    #pH profile
    ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'pH', 'acids_pH.tsv') 
    databaseName = 'bhbtri_analysis_db.sqlite3'
    databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

    #update database with parameters from a file
    ##########################################

    

       #Load database
    db = get_database(os.path.join(databaseFolder, databaseName))

    #getStarting pH
    wc = createMetabolome(db, 'wc_media')
    
    if pHControl is None:
        predictpH = getpH(wc.metabolites, ipH_path)
    else:
        predictpH = mockpHfunc(wc.metabolites,pH=pHControl)
     
    pH =  predictpH(wc.get_concentration())

    #get the feed media and the reactor media
    wc_feed = createMetabolome(db, 'wc_media', pH, pHFunc=predictpH)
    wc_reactor = createMetabolome(db, 'wc_media', pH, pHFunc=predictpH)

    #get the feed obj. Make it sterile
    feed_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc_media'),
                              'bt':createBacteria(db, 'bt', 'wc_media'),
                              'ri':createBacteria(db, 'ri', 'wc_media')})
    feed_microbiome.subpopD['xa'].count = 0
    feed_microbiome.subpopD['xe'].count = 0
    feed_microbiome.subpopD['xi'].count = 0

    #create the reactor obj, with starting populations
    reactor_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc_media'),
                                 'bt':createBacteria(db, 'bt', 'wc_media'),
                                 'ri':createBacteria(db, 'ri', 'wc_media')})
    reactor_microbiome.subpopD['xa'].count = xa
    reactor_microbiome.subpopD['xb'].count = xb
    reactor_microbiome.subpopD['xe'].count = xe
    reactor_microbiome.subpopD['xf'].count = xf
    reactor_microbiome.subpopD['xi'].count = xi
    reactor_microbiome.subpopD['xj'].count = xj
    
    batchA = Pulse(wc_feed, feed_microbiome, 0, 6000, 100, 0, 0, dilutionFactor,dilutionFactor)
    
        #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[
                                                      batchA
                                                       ], 15)
    
    reactor = update_bh_Params("z1_r", 0.0, reactor)
    reactor = update_bh_Params("z2_r", 0.0, reactor)
    
    reactor = update_bt_Params("z6_r", 0.0, reactor)
    reactor = update_bt_Params("z10_r", 0.0, reactor)
    
    reactor = update_ri_Params("z11_r", 0.0, reactor)
    reactor = update_ri_Params("z15_r", 0.0, reactor)
    reactor.simulate()
    
    
    
    batchA = Pulse(wc_feed, feed_microbiome, 0, 6000, 100, 0, 0, dilutionFactor,dilutionFactor)

    #simulate
    reactor = Reactor(reactor_microbiome, wc_reactor,[
                                                  batchA
                                                   ], 15)
    reactor.simulate()
    b = reactor.cellActive_dyn.T[-1]
    bac_composition = b
    pyru = reactor.met_simul[reactor.metabolome.metabolites.index('pyruvate')]
    gluc = reactor.met_simul[reactor.metabolome.metabolites.index('glucose')]
    treh = reactor.met_simul[reactor.metabolome.metabolites.index('trehalose')]
    mann = reactor.met_simul[reactor.metabolome.metabolites.index('mannose')]
    acet = reactor.met_simul[reactor.metabolome.metabolites.index('acetate')]
    form = reactor.met_simul[reactor.metabolome.metabolites.index('formate')]
    lact = reactor.met_simul[reactor.metabolome.metabolites.index('lactate')]
    succ = reactor.met_simul[reactor.metabolome.metabolites.index('succinate')]
    buty = reactor.met_simul[reactor.metabolome.metabolites.index('butyrate')]

    metsA = np.array([pyru[0],
        gluc[0],
        treh[0],
        mann[0],
        acet[0],
        lact[0],
        form[0],
        succ[0],
        buty[0]
        ])
    metsB = np.array([pyru[-1],
        gluc[-1],
        treh[-1],
        mann[-1],
        acet[-1],
        lact[-1],
        form[-1],
        succ[-1],
        buty[-1]
        ])
    pH = reactor.pH_simul[-1]
    return bac_composition,metsA,metsB, pH



#pH_points = np.linspace(5,6.5,30)
simulation_points = 150 
dilution_rate_points = np.linspace(0,3,simulation_points)


bhbtri2_mA = []
bhbtri2_mB = []
bhbtri2_b = []
bhbtri2_pH = []

for d in tqdm(dilution_rate_points):
    
    
    
    
    bac, metA, metB, pH = makeSimulation(pHControl=None,
                                      dilutionFactor=d)
    bhbtri2_b.append(bac)
    bhbtri2_mA.append(metA)
    bhbtri2_mB.append(metB)
    bhbtri2_pH.append(pH)



bh_bhbtri = np.array([i[0] for i in bhbtri2_b])
bt_bhbtri = np.array([i[1] for i in bhbtri2_b])
ri_bhbtri = np.array([i[2] for i in bhbtri2_b])

pyru_bhbtri = np.array([i[0] for i in bhbtri2_mB])
gluc_bhbtri = np.array([i[1] for i in bhbtri2_mB])
treh_bhbtri = np.array([i[2] for i in bhbtri2_mB])
mann_bhbtri = np.array([i[3] for i in bhbtri2_mB])
acet_bhbtri = np.array([i[4] for i in bhbtri2_mB])
lact_bhbtri = np.array([i[5] for i in bhbtri2_mB])
form_bhbtri = np.array([i[6] for i in bhbtri2_mB])
succ_bhbtri = np.array([i[7] for i in bhbtri2_mB])
buty_bhbtri = np.array([i[8] for i in bhbtri2_mB])


pH_bhbtri = np.array(bhbtri2_pH)


dataM = np.array([pyru_bhbtri, 
                  gluc_bhbtri, 
                  treh_bhbtri, 
                  mann_bhbtri, 
                  acet_bhbtri, 
                  lact_bhbtri, 
                  form_bhbtri, 
                  succ_bhbtri, 
                  buty_bhbtri, 
                  pH_bhbtri, 
                  bh_bhbtri, 
                  bt_bhbtri, 
                  ri_bhbtri])

rows = ['pyruvate', 
        'glucose', 
        'trehalose',
        'mannose',
        'acetate',
        'lactate',
        'succinate',
        'formate',
        'butyrate',
        'pH',
        'Blautia hydrogenotrophica',
        'Bacteroides thetaiotaomicron',
        'Roseburia intestinalis'
        ]

create_heatmap(dataM, 
               rows, 
               np.linspace(0,simulation_points-1, 6), 
               np.round(np.linspace(0, 3/15, 6),3), 
               'dilution rate($h^{-1})$', 
               None, 
               None, 
               fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'dilutionSweep.png'))



