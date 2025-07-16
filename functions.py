import numpy as np
import csv
import pandas as pd 
import collections
from numpy import matrix
from numpy.linalg import matrix_rank
import networkx as nx

import matplotlib.pyplot as plt
import os
import heapq
from sympy import *
import time

from scipy.sparse import lil_matrix,save_npz,load_npz
from scipy.linalg import lu
from minBasis import minBasis
from PPH import persHomo
import timeit
from openpyxl import load_workbook
from math import pi,atan

def cos_theta(x,y):
    l_x=np.sqrt(x.dot(x))
    l_y=np.sqrt(y.dot(y))
    dian=x.dot(y)
    cos_=dian/(l_x*l_y)
    return cos_

def get_sgn(x,y):
    if np.cross(x,y)[2]>0:
        return 1
    elif np.cross(x,y)[2]<0:
        return -1
    else:
        return 0

def get_theta(x,y):
    l_x=np.sqrt(x.dot(x))
    l_y=np.sqrt(y.dot(y))
    dian=x.dot(y)
    cos_=dian/(l_x*l_y)
    #sgn=get_sgn(x,y)
    return np.arccos(cos_)

def get_theta1(x,y):
    l_x=np.sqrt(x.dot(x))
    l_y=np.sqrt(y.dot(y))
    dian=x.dot(y)
    cos_=dian/(l_x*l_y)
    if cos_>1:
        cos_=1
    else:
        cos_=cos_
    sgn=get_sgn(x,y)
    return sgn*np.arccos(cos_)

def find_edge(w,xbegin,ybegin,s,L,W,width,mat):
    i1,j1,i2,j2=xbegin-s,ybegin-s,xbegin-s,ybegin-s
    for i in range(xbegin,xbegin+L,s):
        for j in range(ybegin,ybegin+W,s):
            if i+s<xbegin+L and (mat[i*width+j,(i+s)*width+j]==w or mat[(i+s)*width+j,i*width+j]==w):
                i1,j1,i2,j2=i,j,i+s,j
            if j+s<ybegin+W and (mat[i*width+j,i*width+j+s]==w or mat[i*width+j+s,i*width+j]==w):
                i1,j1,i2,j2=i,j,i,j+s
    if [i1,j1,i2,j2]==[xbegin-s,ybegin-s,xbegin-s,ybegin-s]:
        print(i1,j1,i2,j2)
    return i1,j1,i2,j2
    
def find_non_edge4_rightup(i,j,s,width,mat):
    exist=[]
    non=[]
    v1=i*width+j
    v2=(i+s)*width+j
    v3=(i+s)*width+j+s
    v4=i*width+(j+s)
    if mat[v1,v2]==0 and mat[v2,v1]==0:
        non.append([v1,v2])
    else:
        exist.append([v1,v2])
    if mat[v2,v3]==0 and mat[v3,v2]==0:
        non.append([v2,v3])
    else:
        exist.append([v2,v3])
    if mat[v3,v4]==0 and mat[v4,v3]==0:
        non.append([v4,v3])
    else:
        exist.append([v4,v3])
    if mat[v1,v4]==0 and mat[v4,v1]==0:
        non.append([v1,v4])
    else:
        exist.append([v1,v4])
    return non

def find_non_edge4_rightdown(i,j,s,width,mat):
    exist=[]
    non=[]
    v1=i*width+j-s
    v2=(i+s)*width+j-s
    v3=(i+s)*width+j
    v4=i*width+j
    if mat[v1,v2]==0 and mat[v2,v1]==0:
        non.append([v1,v2])
    else:
        exist.append([v1,v2])
    if mat[v2,v3]==0 and mat[v3,v2]==0:
        non.append([v2,v3])
    else:
        exist.append([v2,v3])
    if mat[v3,v4]==0 and mat[v4,v3]==0:
        non.append([v4,v3])
    else:
        exist.append([v4,v3])
    if mat[v4,v1]==0 and mat[v1,v4]==0:
        non.append([v1,v4])
    else:
        exist.append([v1,v4])
    return non

def find_non_edge4_upright(i,j,s,width,mat):
    exist=[]
    non=[]
    v1=i*width+j
    v2=(i+s)*width+j
    v3=(i+s)*width+j+s
    v4=i*width+j+s
    if mat[v1,v2]==0 and mat[v2,v1]==0:
        non.append([v1,v2])
    else:
        exist.append([v1,v2])
    if mat[v2,v3]==0 and mat[v3,v2]==0:
        non.append([v2,v3])
    else:
        exist.append([v2,v3])
    if mat[v3,v4]==0 and mat[v4,v3]==0:
        non.append([v4,v3])
    else:
        exist.append([v4,v3])
    if mat[v4,v1]==0 and mat[v1,v4]==0:
        non.append([v1,v4])
    else:
        exist.append([v1,v4])
    return non

def find_non_edge4_upleft(i,j,s,width,mat):
    exist=[]
    non=[]
    v1=(i-s)*width+j
    v2=i*width+j
    v3=i*width+j+s
    v4=(i-s)*width+j+s
    if mat[v1,v2]==0 and mat[v2,v1]==0:
        non.append([v1,v2])
    else:
        exist.append([v1,v2])
    if mat[v2,v3]==0 and mat[v3,v2]==0:
        non.append([v2,v3])
    else:
        exist.append([v2,v3])
    if mat[v3,v4]==0 and mat[v4,v3]==0:
        non.append([v4,v3])
    else:
        exist.append([v4,v3])
    if mat[v4,v1]==0 and mat[v1,v4]==0:
        non.append([v1,v4])
    else:
        exist.append([v1,v4])
    return non


def weight_point(listin,VF,mat):
    n=len(listin)
    ecx=[]
    ecy=[]
    W=[]
    Wsum=0
    for i in range(n):
        ecx.append((VF.iloc[listin[i%n],1]+VF.iloc[listin[(i+1)%n],1])/2)
        ecy.append((VF.iloc[listin[i%n],2]+VF.iloc[listin[(i+1)%n],2])/2)
        if mat[listin[i%n],listin[(i+1)%n]]!=0:
            W.append(mat[listin[i%n],listin[(i+1)%n]])
            Wsum=Wsum+mat[listin[i%n],listin[(i+1)%n]]
        else:
            W.append(-mat[listin[(i+1)%n],listin[i%n]])
            Wsum=Wsum-mat[listin[(i+1)%n],listin[i%n]]
    xi,yi=0,0
    for i in range(len(ecx)):
        xi=xi+W[i]*ecx[i]
        yi=yi+W[i]*ecy[i]
    #print(xi/Wsum,yi/Wsum)
    return xi/Wsum,yi/Wsum

def compute_sqar(mat,v1,v2,v3,v4,VF,bps,min_cyc,centers,eps,latbegin,lonbegin):
    flag=0
    V1=np.array([VF.iloc[v1,4],VF.iloc[v1,5],VF.iloc[v1,6]])
    V2=np.array([VF.iloc[v2,4],VF.iloc[v2,5],VF.iloc[v2,6]])
    V3=np.array([VF.iloc[v3,4],VF.iloc[v3,5],VF.iloc[v3,6]])
    V4=np.array([VF.iloc[v4,4],VF.iloc[v4,5],VF.iloc[v4,6]])
    t=(get_theta1(V1,V2)+get_theta1(V2,V3)+get_theta1(V3,V4)+get_theta1(V4,V1))/(2*pi)
    #print(v1,v2,v3,v4,t)
    if t==t:
        index=round(t)
    else:
        index=0
    if (([v1,v2,v3,v4,index] not in min_cyc)) and index<0:
        min_cyc.append([v1,v2,v3,v4,index])
        bps.append((weight_point([v1,v2,v3,v4],VF,mat)))
        rel= [round(eps*weight_point([v1,v2,v3,v4],VF,mat)[1]+lonbegin,2),round(eps*weight_point([v1,v2,v3,v4],VF,mat)[0]+latbegin,2)]
        print(rel)


