'''
Tracking and comparing the time-varying vector fields
Compute the bottleneck distance between the persistence diagrams of the vector fields
rig: '2023sl' or '2023kn' means tracking typhoon Kanun or Saola
'''
import numpy as np
import csv
import pandas as pd 
import collections
from numpy import matrix
from numpy.linalg import matrix_rank
import networkx as nx
import os
import heapq
from sympy import *
import gudhi
from gudhi.hera import wasserstein_distance,bottleneck_distance
from scipy.sparse import lil_matrix,save_npz,load_npz
from scipy.linalg import lu
from minBasis import minBasis
from PPH import persHomo
import timeit
from openpyxl import load_workbook
from math import pi,atan
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()

rig='2023kn'
if not os.path.exists(SCRIPT_DIR / rig / 'wfs'):
    print(f"Error: path not found!")


dgm_old=[]
dgm_new=[]
word_label=[]
bot=[]
count=0
old_word='None'
path=SCRIPT_DIR / rig / 'wfs'
files = sorted(os.listdir(path))
for fi in files:
    word=fi[:10]
    word_label.append(word)

    dgm_new=np.loadtxt(SCRIPT_DIR / rig /f"{word}_result"/f"{word}_pair.txt")
    
    if count==0:
        bot_dis=0
    else:
        bot_dis=bottleneck_distance(dgm_old,dgm_new)
        bot.append(bot_dis)

    print(str(old_word)+' ——> '+ word)
    print('Bot_dis = '+str(bot_dis))
    print('**********************************')

    dgm_old=dgm_new
    count+=1
    old_word=word

np.savetxt(SCRIPT_DIR / rig / f'bot_{rig}.txt',bot)
if bot:
    plt.figure(figsize=(10,6))
    plt.plot(bot, label='Bottleneck Distance', marker='o')
    plt.xticks([i for i in range(len(bot))],word_label[1:],rotation=70)
    plt.ylim(0.0,2.0)
    plt.title('Bottleneck Distance for Khanun dataset')
    plt.legend()
    plt.show()
else:
    print("No data to plot.")
