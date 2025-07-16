'''
Tracking the typhoon centers from vector field data
rig: '2023sl' or '2023kn' means tracking typhoon Kanun or Saola
'''
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

def find_edge(w):
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
    
def find_non_edge4_rightup(i,j):
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

def find_non_edge4_rightdown(i,j):
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

def find_non_edge4_upright(i,j):
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

def find_non_edge4_upleft(i,j):
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

def weight_point(listin):
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

def compute_sqar(v1,v2,v3,v4,bps,min_cyc,centers):
    flag=0
    V1=np.array([VF.iloc[v1,4],VF.iloc[v1,5],VF.iloc[v1,6]])
    V2=np.array([VF.iloc[v2,4],VF.iloc[v2,5],VF.iloc[v2,6]])
    V3=np.array([VF.iloc[v3,4],VF.iloc[v3,5],VF.iloc[v3,6]])
    V4=np.array([VF.iloc[v4,4],VF.iloc[v4,5],VF.iloc[v4,6]])
    t=(get_theta1(V1,V2)+get_theta1(V2,V3)+get_theta1(V3,V4)+get_theta1(V4,V1))/(2*pi)
    if t==t:
        index=round(t)
    else:
        index=0
    if (([v1,v2,v3,v4,index] not in min_cyc)) and index!=0:
        min_cyc.append([v1,v2,v3,v4,index])
        bps.append((weight_point([v1,v2,v3,v4])))
        
        # the weight_point here has the coordinate of (latitude,longitude) according to the data format from CCMP
        centers.append((round(eps*weight_point([v1,v2,v3,v4])[1]+lonbegin,2),round(eps*weight_point([v1,v2,v3,v4])[0]+latbegin,2)))

def VectorfieldtoDigraph_arc(VF,length,width,high,eps):
    VectorField=pd.read_csv(VF,header=None)
    VectorField=pd.DataFrame(VectorField)

    n=VectorField.shape[0]
    Digraph_mat=lil_matrix((n,n),dtype='float32')

    max_edge=0
    for xi in range(xbegin,xbegin+L,s):
        for yi in range(ybegin,ybegin+W,s):
            i=width*high*xi+high*yi
            vxi=VectorField.iloc[i,4]
            vyi=VectorField.iloc[i,5]

            if xi+s<xbegin+L:
                xj=xi+s
                yj=yi
                j=width*high*xj+high*yj
                vxj=VectorField.iloc[j,4]
                vyj=VectorField.iloc[j,5]
                arr1=np.array([vxi,vyi,0])
                arr2=np.array([vxj,vyj,0])
                if (vxi==0 and vyi==0) or (vxj==0 and vyj==0):
                    Digraph_mat[i,j]=0
                    Digraph_mat[j,i]=0
                else:
                    th=get_theta(arr1,arr2)
                    w_x=get_sgn(arr1,arr2)*th
                    if abs(abs(th)-pi)==0 or th==0 or abs(w_x)>eps:
                        Digraph_mat[i,j]=0
                        Digraph_mat[j,i]=0
                    elif w_x>0:
                        Digraph_mat[i,j]=w_x
                    else:
                        Digraph_mat[j,i]=-w_x
                    max_edge=max(max_edge,w_x)
            if yi+s<ybegin+W:
                yj=yi+s
                xj=xi
                j=width*high*xj+high*yj
                vxj=VectorField.iloc[j,4]
                vyj=VectorField.iloc[j,5]
                arr1=np.array([vxi,vyi,0])
                arr2=np.array([vxj,vyj,0])
                if (vxi==0 and vyi==0) or (vxj==0 and vyj==0):
                    Digraph_mat[i,j]=0
                    Digraph_mat[j,i]=0
                else:
                    th=get_theta(arr1,arr2)
                    w_y=get_sgn(arr1,arr2)*th
                    if abs(abs(th)-pi)==0 or th==0 or abs(w_y)>eps:
                        Digraph_mat[i,j]=0
                        Digraph_mat[j,i]=0
                    elif w_y>0:
                        Digraph_mat[i,j]=w_y
                    else:
                        Digraph_mat[j,i]=-w_y
                    max_edge=max(max_edge,w_y)

    Digraph_mat=Digraph_mat.tocoo()
    save_npz(f'{rig}\\{word}_result\\{word}_Digraph_mat.npz',Digraph_mat)
    #print('get matrix!')
    return max_edge

rig='2023sl'
path=f'{rig}\\wfs'
files = sorted(os.listdir(path))
for fi in files:
    word=fi[:10]
    print(word)
    
    if not os.path.exists(f'{rig}\\{word}_result'):
        os.makedirs(f'{rig}\\{word}_result')
    
    
    VF= path+'\\'+fi
    ls=pd.read_csv(VF,header=None).iloc[:,1]
    ws=pd.read_csv(VF,header=None).iloc[:,2]
    all_len=pd.read_csv(VF,header=None).shape[0]
    #print(max(ls),max(ws))
    s=1
    eps=0.25 
    if rig=='2023sl':
        lonbegin,latbegin=113.125,13.125  #sl
    if rig=='2023kn':
        lonbegin,latbegin=123.625,21.125  #kn

    length,width,high=max(ls)+1,max(ws)+1,1
    L,W=length,width
    xbegin,ybegin=0,0
    max_edge=VectorfieldtoDigraph_arc(VF,length,width,high,10)
    mat=load_npz(f'{rig}\\{word}_result\\{word}_Digraph_mat.npz')
    mat=mat.tolil()

    G=nx.from_scipy_sparse_matrix(mat, create_using=nx.DiGraph)
    g=persHomo(G)
    max_edge=3.15
    g.perHom(max_edge)
    
    e=[]
    min_cyc=[]
    bps=[]
    centers=[]
    pairs=np.array(g.pair)
    for i in range(len(pairs)):
        if pairs[i][1]==max_edge:
            e.append(pairs[i][0])

    VF=pd.read_csv(path+'\\'+fi,header=None)
    while len(e)>0:
        i1,j1,i2,j2=find_edge(e[0])
        w=max(mat[i1*width*high+j1*high,i2*width*high+j2*high],mat[i2*width*high+j2*high,i1*width*high+j1*high])
        flag=0
        if i2-i1==s:
            i=i1
            j=j1
            v1,v2,v3,v4=i*width*high+j*high,(i+s)*width*high+j*high,(i+s)*width*high+(j+s)*high,i*width*high+(j+s)*high
            if i+s<xbegin+L and j+s<ybegin+W:
                non=find_non_edge4_rightup(i,j)
                if len(non)==0:
                    compute_sqar(v1,v2,v3,v4,bps,min_cyc,centers)
            
            v1,v2,v3,v4=i*width*high+(j-s)*high,(i+s)*width*high+(j-s)*high,(i+s)*width*high+j*high,i*width*high+j*high
            if i+s<xbegin+L and j-s>=ybegin:
                non=find_non_edge4_rightdown(i,j)
                if len(non)==0:
                    compute_sqar(v1,v2,v3,v4,bps,min_cyc,centers)

        if j2-j1==s:
            i=i1
            j=j1
            v1,v2,v3,v4=i*width*high+j*high,(i+s)*width*high+j*high,(i+s)*width*high+(j+s)*high,i*width*high+(j+s)*high
            if i+s<xbegin+L and j+s<ybegin+W:
                non=find_non_edge4_upright(i,j)
                if len(non)==0:
                    compute_sqar(v1,v2,v3,v4,bps,min_cyc,centers)

            v1,v2,v3,v4=(i-s)*width*high+j*high,i*width*high+j*high,i*width*high+(j+s)*high,(i-s)*width*high+(j+s)*high
            if i-s>=xbegin and j+s<ybegin+W:
                non=find_non_edge4_upleft(i,j)
                if len(non)==0:
                    compute_sqar(v1,v2,v3,v4,bps,min_cyc,centers)

        del(e[0])
            

    print("singularities:",centers)
    
 