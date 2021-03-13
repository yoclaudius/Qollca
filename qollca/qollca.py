import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def sigma_5():
    pass

class Muestra_control():

    def m_c():
        pass

    def identificacion_cruzada():
        pass


def monte_carlo(datos, bines, rango, grafico=False, test=False):
    #Histograma para los datos
    pesos=np.ones_like(datos)/float(len(datos))
    histograma=np.histogram(datos,bins=bines,range=rango,weights=pesos)        
    #Distribucion acumulada
    acumulada=np.cumsum(histograma[0])
    #Se asigna a cada valor aleatorio de la tabla 'Simulacion' su bin 
    #correspondiente
    valaleatorios = np.random.random(300000)
    valencadabin = np.zeros(len(histograma[0]))
    for i in range(len(valencadabin)-1):
        valencadabin[i] = np.sum(np.array([valaleatorios >= acumulada[i]])
                             * np.array([valaleatorios < acumulada[i+1]]))
    valencadabin = np.array(valencadabin, dtype=int)
#Valores aleatorios que forman la distribucion simulada 
    simulacion =( np.repeat(histograma[1][:-1], valencadabin) + 
                 np.random.random() * (histograma[1][1] - histograma[1][0]) )
    
    if grafico == True:
        plt.figure(figsize=(1.4*6.4,1.4*4.8),constrained_layout=True)
        bins=np.linspace(rango[0],rango[1],bines)
        pesos1=(np.ones_like(datos)/float(len(datos)))
        plt.hist(datos,bins,weights=pesos1,color='blue',histtype='step',
                 lw=3,alpha=1,label='O3HB')
        pesos2=(np.ones_like(simulacion)/
                        float(len(simulacion)))
        plt.hist(simulacion,bins,weights=pesos2,
                 color='red',histtype='step',ls='--',lw=3,alpha=1,label='Simulacion')
        
    if test == False:
        return simulacion
    else:
        return histograma


def catalog_corr():
    pass

def bootstrap(datos):
    pass


def jackknife():
    pass

def Luminosity_function():
    pass




