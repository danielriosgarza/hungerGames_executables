# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 10:48:30 2025

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


import cmasher as cmr

from tqdm import tqdm



import plotly.io as pio
pio.renderers.default='browser'

from matplotlib.colors import ListedColormap

#pH profile
ipH_path = os.path.join(Path(os.getcwd()).parents[1], 'files', 'pH', 'acids_pH.tsv') 
databaseName = 'bhbtri_analysis_db.sqlite3'
databaseFolder =  os.path.join(Path(os.getcwd()).parents[1], 'files', 'dbs')

#Load database
db = get_database(os.path.join(databaseFolder, databaseName))

# Create metabolome from media table in your database
wc_media = createMetabolome(db, mediaName='wc_media', pH0=6.5)

# Create function that predicts pH (elastic net)
pHfunction = getpH(wc_media.metabolites, ipH_path)
initial_pH = pHfunction(wc_media.get_concentration())
print(f"initial pH: {initial_pH}")

def get_reactor(species = ['bh'], dilutionFactor = 0, pulses = None, pHFunc = pHfunction, time = None):
  if time is None:
    if dilutionFactor:
      time = 6000
    else:
      time = 120
  else:
    time = time
  #create media for reactor
  wc_media_reactor = createMetabolome(db, mediaName='wc_media', pH0 = initial_pH, pHFunc = pHFunc)

  #Create media for feed
  wc_media_feed = createMetabolome(db, mediaName='wc_media', pH0 = initial_pH, pHFunc = pHFunc)

  #create microbiome for reactor
  microbiome_reactor = Microbiome({'bh': createBacteria(db, speciesID='bh', mediaName='wc_media'),
                               'bt': createBacteria(db, speciesID='bt', mediaName='wc_media'),
                               'ri': createBacteria(db, speciesID='ri', mediaName='wc_media')})

  if 'bh' not in species:
    microbiome_reactor.subpopD['xa'].count = 0

  if 'bt' not in species:
    microbiome_reactor.subpopD['xe'].count = 0

  if 'ri' not in species:
    microbiome_reactor.subpopD['xi'].count = 0



  #create microbiome for feed
  microbiome_feed = Microbiome({'bh': createBacteria(db, speciesID='bh', mediaName='wc_media'),
                               'bt': createBacteria(db, speciesID='bt', mediaName='wc_media'),
                               'ri': createBacteria(db, speciesID='ri', mediaName='wc_media')})

  for subpop in microbiome_feed.subpopD:
      microbiome_feed.subpopD[subpop].count=0

  # Define pulse and reactor
  if pulses is None:
    pulse = Pulse(wc_media_feed, microbiome_feed, 0, time, 1000, 0, 0, dilutionFactor, dilutionFactor)
    pulses = [pulse]
  reactor = Reactor(microbiome_reactor, wc_media_reactor, pulses, 15)

  return reactor


def makeKineticPlot(x,
                    y,
                    color,
                    legend,
                    xlabel,
                    ylabel,
                    title = None,
                    linestyle = 'o-',
                    legendSize = 18):
    '''
    '$10^5$ cells/uL'
    'Time (h)'
    'mM'

    '''

    #fig, ax = plt.subplots()


    plt.plot(x, y, color=color, linestyle = linestyle, lw=3, alpha = 1.0, label=legend)

    if legend is not None:
      legend_properties = {'size': legendSize}
      plt.legend(prop=legend_properties,
                 loc='upper right',        # <-- THIS is the key!
                 bbox_to_anchor=(1.0, 1.0),# anchor at (x=1, y=1) in Axes coords
                 bbox_transform=plt.gca().transAxes,
                 borderaxespad=0.0,        # optional: no extra padding
                 frameon=True)             # optional: keep the black box
    ax = plt.gca()

    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)


    if title is not None:
        plt.title(title, fontsize=16)





    plt.tight_layout()


def plotSpecies(simulated_reactor_obj, withDilution = True):
  if withDilution:
    time = simulated_reactor_obj.time_simul*0.1
  else:
    time = simulated_reactor_obj.time_simul

  makeKineticPlot(x = time,
                  y = simulated_reactor_obj.cellActive_dyn[0],
                  color = '#FF10F0',
                  legend = '$\it{Blautia \ hydrogenotrophica}$',
                  xlabel = 'time (h)',
                  ylabel = '$10^5$ cells/uL',
                  title = None,
                  linestyle = '-',
                  legendSize = 10)

  makeKineticPlot(x = time,
                  y = simulated_reactor_obj.cellActive_dyn[1],
                  color = '#ff8300',
                  legend = '$\it{Bacteroides \ thetaiotaomicron}$',
                  xlabel = 'time (h)',
                  ylabel = '$10^5$ cells/uL',
                  title = None,
                  linestyle = '-',
                  legendSize = 10)
  makeKineticPlot(x = time,
                  y = simulated_reactor_obj.cellActive_dyn[2],
                  color = '#00B8FF',
                  legend = '$\it{Roseburia \ instestinalis}$',
                  xlabel = 'time (h)',
                  ylabel = '$10^5$ cells/uL',
                  title = None,
                  linestyle = '-',
                  legendSize = 10)
  plt.show()
  
  
#################Add collab code here##########
#the dilution factor is the volume removed and added in 1h, since we consider the volume = 15 mL, the rate is [dilution factor]/volume

def pH_f(metabolome):
  return 5.60

reactor = get_reactor(species = ['bh', 'bt', 'ri'], dilutionFactor=1, pHFunc=pH_f)

reactor.microbiome.subpopD['xa'].count = 0.003
reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
reactor.microbiome.subpopD['xj'].count = 0.0001


reactor = update_bh_Params("z1_r", 0.0, reactor)
reactor = update_bh_Params("z2_r", 0.0, reactor)
reactor = update_bt_Params("z6_r", 0.0, reactor)
reactor = update_bt_Params("z10_r", 0.0, reactor)
reactor = update_ri_Params("z11_r", 0.0, reactor)
reactor = update_ri_Params("z15_r", 0.0, reactor)
reactor.simulate()

#to avoid warning from the colab notebook
mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})
plotSpecies(reactor)