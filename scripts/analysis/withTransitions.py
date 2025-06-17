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

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()
import seaborn as sns

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
    


def plotSpecies(simulated_reactor_obj, 
                withDilution = True,
                fileName = None):
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
  if fileName is not None:
      plt.savefig(fileName, dpi = 600, transparent = True)
  plt.show()

  
  
#################Add collab code here##########
#the dilution factor is the volume removed and added in 1h, since we consider the volume = 15 mL, the rate is [dilution factor]/volume

def pH_f(metabolome):
  return 5.60

reactor = get_reactor(species = ['bh', 'bt', 'ri'], dilutionFactor=0.615, time =3000)

reactor.microbiome.subpopD['xa'].count = 0.003
#reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
#reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
#reactor.microbiome.subpopD['xj'].count = 0.0001



reactor.simulate()

#to avoid warning from the colab notebook
#mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', '0.041_with_transition.png')
plotSpecies(reactor, fileName=fileName)
#####################################################################

reactor = get_reactor(species = ['bh', 'bt', 'ri'], dilutionFactor=0.585, time =3000)

reactor.microbiome.subpopD['xa'].count = 0.003
#reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
#reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
#reactor.microbiome.subpopD['xj'].count = 0.0001



reactor.simulate()

#to avoid warning from the colab notebook
#mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', '0.039_with_transition.png')
plotSpecies(reactor, fileName=fileName)

##################################################################################

dilutionFactor = 0.041*15


#getStarting pH
wc = createMetabolome(db, 'wc_media')

predictpH = getpH(wc.metabolites, ipH_path)
pH =  predictpH(wc.get_concentration())


wc_feed = createMetabolome(db, 'wc_media', pH, pHFunc=predictpH)
feed_microbiome = Microbiome({'bh':createBacteria(db, 'bh', 'wc_media'),
                              'bt':createBacteria(db, 'bt', 'wc_media'),
                              'ri':createBacteria(db, 'ri', 'wc_media')})
feed_microbiome.subpopD['xa'].count = 0
feed_microbiome.subpopD['xe'].count = 0
feed_microbiome.subpopD['xi'].count = 0


pulse1 = Pulse(wc_feed, feed_microbiome, 0, 2400, 100, 0, 0, dilutionFactor,dilutionFactor)

pulse2 = Pulse(wc_feed, feed_microbiome, 2400, 2640, 100, 0, 0, 0,0)

pulse3 = Pulse(wc_feed, feed_microbiome, 2640, 5040, 100, 0, 0, dilutionFactor,dilutionFactor)



reactor = get_reactor(species = ['bh', 'bt', 'ri'],
                      dilutionFactor=dilutionFactor,
                      pulses=[pulse1, pulse2, pulse3])

reactor.microbiome.subpopD['xa'].count = 0.003
#reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
#reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
#reactor.microbiome.subpopD['xj'].count = 0.0001


reactor.simulate()
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'dilution_histeresis_with_transition.png')
plotSpecies(reactor)

##########################################

def pH_f(metabolome):
  return 5.47

reactor = get_reactor(species = ['bh', 'bt', 'ri'], dilutionFactor=1.0, time =3000, pHFunc = pH_f)

reactor.microbiome.subpopD['xa'].count = 0.003
#reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
#reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
#reactor.microbiome.subpopD['xj'].count = 0.0001



reactor.simulate()

#to avoid warning from the colab notebook
#mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'pH_5.47_with_transition.png')
plotSpecies(reactor, fileName=fileName)

##########################################

def pH_f(metabolome):
  return 5.50

reactor = get_reactor(species = ['bh', 'bt', 'ri'], dilutionFactor=1.0, time =3000, pHFunc = pH_f)

reactor.microbiome.subpopD['xa'].count = 0.003
#reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
#reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
#reactor.microbiome.subpopD['xj'].count = 0.0001



reactor.simulate()

#to avoid warning from the colab notebook
#mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'pH_5.50_with_transition.png')
plotSpecies(reactor, fileName=fileName)



##########################################

def pH_f(metabolome):
  return 5.60

reactor = get_reactor(species = ['bh', 'bt', 'ri'], dilutionFactor=1.0, time =3000, pHFunc = pH_f)

reactor.microbiome.subpopD['xa'].count = 0.003
#reactor.microbiome.subpopD['xb'].count = 0.0001
reactor.microbiome.subpopD['xe'].count = 0.003
#reactor.microbiome.subpopD['xf'].count = 0.0001
reactor.microbiome.subpopD['xi'].count = 0.003
#reactor.microbiome.subpopD['xj'].count = 0.0001



reactor.simulate()

#to avoid warning from the colab notebook
#mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['DejaVu Sans']})
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'pH_5.60_with_transition.png')
plotSpecies(reactor, fileName=fileName)



################################################################
def mockpHfunc(metObj, pH=6.5):
    def pHfunc(metObj):
        return pH
    return pHfunc

fixedpH1 = 5.60
fixedpH2 = 5.47

#getStarting pH
wc = createMetabolome(db, 'wc_media')
predictpH = getpH(wc.metabolites, ipH_path)
fixedpHA = mockpHfunc(wc.metabolites,pH=fixedpH1)
fixedpHC = mockpHfunc(wc.metabolites,pH=fixedpH2)
pH =  predictpH(wc.get_concentration())

#get the feed media and the reactor media
wc_feedA = createMetabolome(db, 'wc_media', pH, pHFunc=fixedpHA)
wc_feedC = createMetabolome(db, 'wc_media', pH, pHFunc=fixedpHC)

wc_reactorA = createMetabolome(db, 'wc_media', pH, pHFunc=fixedpHA)
wc_reactorC = createMetabolome(db, 'wc_media', pH, pHFunc=fixedpHC)


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
reactor_microbiome.subpopD['xa'].count = 0.003
reactor_microbiome.subpopD['xe'].count = 0.003
reactor_microbiome.subpopD['xi'].count = 0.003


d = 1.0




batchA = Pulse(wc_feedA, feed_microbiome, 0, 600, 100, 0, 0, d,d)

batchB = Pulse(wc_feedC, feed_microbiome, 600, 3000, 100, 0, 0, d,d)


#simulate
reactorA = Reactor(reactor_microbiome, wc_reactorA,[batchA], 15)



reactorA.simulate()
#reactorA.makePlots()



for i in wc_reactorA.metD:
    wc_reactorC.metD[i].update(wc_reactorA.metD[i].concentration)


reactorB = Reactor(reactorA.microbiome, wc_reactorC,[batchB], 15)

reactorB.simulate()


print(f"pH {fixedpH1} and {fixedpH2}")


makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[0],
                color = '#FF10F0',
                legend = 'Blautia hydrogenotrophica',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[0],
                color = '#FF10F0',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)




makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[1],
                color = '#ff8300',
                legend = 'Bacteroides thetaiotaomicron',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[1],
                color = '#ff8300',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)



makeKineticPlot(x = reactorA.time_simul*0.1,
                y = reactorA.cellActive_dyn[2],
                color = '#00B8FF',
                legend = 'Roseburia instestinalis',
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)

makeKineticPlot(x = reactorB.time_simul*0.1,
                y = reactorB.cellActive_dyn[2],
                color = '#00B8FF',
                legend = None,
                xlabel = 'time (h)',
                ylabel = '$10^5$ cells/uL',
                title = None,
                linestyle = '-',
                legendSize = 10)


title = 'stateA'
fileName = os.path.join(Path(os.getcwd()).parents[1], 'files', 'Figures', 'pH_hysteresis_with_transition.png')
plt.savefig(fileName, transparent=True, dpi=600)
plt.show()