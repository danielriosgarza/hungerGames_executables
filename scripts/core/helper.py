# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:09:43 2023

@author: danie
"""

import os
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from aquarel import load_theme

theme = load_theme("boxy_light")
theme.apply()




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
    
        legend_properties = {'size':legendSize}
        plt.legend(fontsize=16, prop=legend_properties, bbox_to_anchor=(1.0, 1.0))
    
    ax = plt.gca()

    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    
    
    if title is not None:
        plt.title(title, fontsize=16)
    
    
    

    
    plt.tight_layout()
    

    
