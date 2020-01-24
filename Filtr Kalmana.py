#!/usr/bin/env python
# coding: utf-8

# In[331]:


import matplotlib.pyplot as plt
import numpy as np
import random 
import pandas as pd 
from io import StringIO
import math


# In[344]:



dt=0.1
#t=np.linspace(0,10,100) # Parametr symulacji
t=np.linspace(0,200,2000)
std=10 # odchylenie standardowe

V = np.array([[std*std*dt, 0 ], [0 ,std*std*dt]])
W = std*std


N=2000

ekg = np.loadtxt("ekg.txt", unpack=True)     

#parametry sygnału 
f0=50
fs=250
A1=10
fi=5
wo=(2*3.14*f0)/fs;

x=np.zeros(N)
ekg_sz=np.zeros(N)

#Generacja sygnału zaszumionego
for i in range(N):
    x[i]=np.array(A1*math.cos(wo*i+fi))
    ekg_sz=x+ekg #[:100]


#Warunki początkowe wraz z warunkami określającymi filtr
A = np.array([[2*math.cos(wo), -1],[1 , 0]])
B = np.array([1,0])
C = np.array([1,0])


# In[345]:


# Wartosci poczatkowe
x0 = np.array([0, 0])
P0 = np.array([[1, 0], [0, 1]])
xpri = x0
Ppri = P0
xpost = x0
Ppost = P0
xpost2=x0


# In[346]:


Y = np.zeros(N)
Yf =np.zeros(N)


# In[347]:



for i in range(N):

   Y[i] = ekg_sz[i]
   
   if i > 1:
        
       xpri = A*xpost
       Ppri = A*Ppost*np.transpose(A)+V
       eps = Y[i]-C*xpri
       S = C*Ppri*np.transpose(C) + W
       K = Ppri*C*S**(-1)
       xpost = xpri + K*eps
       Ppost = Ppri - K*S*np.transpose(K)   
       xpost2=xpost.flatten()
       Yf[i] = xpost2[0]

                   


# In[348]:



plt.plot(t, Y, 'blue')
plt.plot(t, Yf, 'red')
plt.legend(('Sygnał Zaszumiony', 'Sygnał Odfiltrowany'),loc='upper right')


# In[ ]:




