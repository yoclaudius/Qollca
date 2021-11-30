import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext

def sigma_5():
    pass

class Muestra_control():

    def m_c():
        pass

    def r_c():
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

    
def bootstrap_error(columna):
    with NumpyRNGContext(1):
        bootresult =bootstrap(columna,100,bootfunc=np.mean)    
    error=bootresult.std()        
    return error    


def jackknife(columna):
    numdatos = len(columna)
    medias = np.zeros(100)
    for i in range(100):
        medias[i] = (np.random.choice(columna, numdatos-1)).mean()
    error = medias.std(ddof=1)
    return error


def Luminosity_function():
    pass


def hist2dmatriz(columna_x, columna_y, bines_x, bines_y):
    """
    

    Parameters
    ----------
    columna_x : array_like
        DESCRIPTION.
    columna_y : array_like
        DESCRIPTION.
    bines_x : int
        DESCRIPTION.
    bines_y : int
        DESCRIPTION.

    Returns
    -------
    matriz : ndarray
        DESCRIPTION.

    Example
    -------
    >>> import qollca
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> x = np.random.normal(size=1000)
    >>> y = np.random.normal(size=1000)
    >>> matriz = qollca.hist2dmatriz(x,y,10,10)
    >>> plt.contourf(matriz)

    """

    
    limbines_x = np.linspace(columna_x.min(), columna_x.max(), bines_x + 1)
    limbines_y = np.linspace(columna_y.min(), columna_y.max(), bines_y + 1)
    matriz = np.zeros([bines_y,bines_x])
    for liminfcol_i, limsupcol_i, j in zip(limbines_x[:-1],limbines_x[1:],range(bines_x + 1)):
        for liminffil_i, limsupfil_i, i in zip(limbines_y[:-1],limbines_y[1:],range(bines_y + 1)):
            a = (columna_x > liminfcol_i) * (columna_x < limsupcol_i)
            b = (columna_y > liminffil_i) * (columna_y < limsupfil_i)
            mascara = a * b
            matriz[i,j] = np.sum(mascara)            
    return matriz        


