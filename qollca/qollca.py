import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext

def sigma_5():
    pass

def Muestra_control():

    def m_c():
        pass

    def r_c():
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
        #Gráficos 
        #masa estelar
        Columna1=Herg['SMASS_MEDIAN']
        Columna2=noAGNs['SMASS_MEDIAN']
        Columna1=Columna1[Columna1<Columna1.quantile(q=0.99)]
        Columna1=Columna1[Columna1>Columna1.quantile(q=0.01)]
        Columna2=Columna2[Columna2<Columna2.quantile(q=0.99)]
        Columna2=Columna2[Columna2>Columna2.quantile(q=0.01)]
        CantidadBines=int((Columna1.max()-Columna1.min())*((len(Columna1))**(1/3))/
                    (3.5*Columna1.std())+0.5)
        pesos1=np.ones_like(Columna1)/float(len(Columna1))
        pesos2=np.ones_like(Columna2)/float(len(Columna2))
        plt.figure(figsize=(1.4*6.4,1.4*4.8),constrained_layout=True)
        bines=np.linspace(Columna1.min(),Columna1.max(),CantidadBines)
        plt.hist(Columna1,bins=bines,weights=pesos1, color='b',histtype='step',
                lw=3,alpha=0.75,label='HERGs')
        plt.hist(Columna2,bins=bines,weights=pesos2, color='indigo',histtype='step',
                linestyle=':',lw=4,alpha=1,label='ControlHERGs')
        plt.legend(loc='upper right',fontsize=16,frameon=False)
        plt.xlabel(r'$log_{10}(M_{\ast} / M_{\odot})$',fontsize=16)
        plt.tick_params(labelsize=16)
        plt.show()

        """        
        
        #Se establecen rangos de masa y luminosidad
        deltaZ=0.02
        deltaMag=0.2
        deltaMass=0.15
        deltaD4000=0.05

        #Se agragan, a las AGNs, columnas con los limites en redshift, magnitud, masa y
        #d4000
        Herg['Zmin']=Herg['Z']-deltaZ
        Herg['Zmax']=Herg['Z']+deltaZ
        Herg['Magmin']=Herg['ABSR']-deltaMag
        Herg['Magmax']=Herg['ABSR']+deltaMag
        Herg['Massmin']=Herg['SMASS_MEDIAN']-deltaMass
        Herg['Massmax']=Herg['SMASS_MEDIAN']+deltaMass
        Herg['D4000min']=Herg['D4000']-deltaD4000
        Herg['D4000max']=Herg['D4000']+deltaD4000
        #Se agrega una columna para indicar la galaxia AGN de la cual es galaxia
        #de control
        noAGNs['EsControl_oNo']=np.zeros(len(noAGNs))
        Herg['HERGSeleccionada_oNo']=np.zeros(len(Herg))
        noAGNs['HERGdeLaQueEsControl']=np.zeros(len(noAGNs))

        np.random.seed(seed=1)
        #for i in range(int(len(Herg))):
        def fun1(i):   
            #SELECCIONO LOS INDICES DEL GRUPO DE HERGs SIN SELECCIONAR
            HergsNoSelec=(Herg.groupby(['HERGSeleccionada_oNo']).get_group(0)).index
            #ELIJO ALEATORIAMENTE UNA DE LAS HERG SIN CONTROL
            indiceHERG=np.random.choice(HergsNoSelec)
            #Marco a la HERG seleccionada
            Herg['HERGSeleccionada_oNo'].iloc[indiceHERG]=1    
            #SELECCIONON EL GRUPO DE GALAXIAS NO ACTIVAS SIN SELECCIONAR
            noAGNsSinSelec=noAGNs.groupby(['EsControl_oNo']).get_group(0)
            #Se indican las galaxias que sirven como control
            Var=noAGNsSinSelec[((noAGNsSinSelec['SMASS_MEDIAN']>=
                                Herg['Massmin'].iloc[indiceHERG])&
                                (noAGNsSinSelec['SMASS_MEDIAN']<=
                                Herg['Massmax'].iloc[indiceHERG])&
                                (noAGNsSinSelec['ABSR']>=
                                Herg['Magmin'].iloc[indiceHERG])&
                                (noAGNsSinSelec['ABSR']<=
                                Herg['Magmax'].iloc[indiceHERG])&
                                (noAGNsSinSelec['Z']>=
                                Herg['Zmin'].iloc[indiceHERG])&
                                (noAGNsSinSelec['Z']<=
                                Herg['Zmax'].iloc[indiceHERG])&
                                (noAGNsSinSelec['D4000']>=
                                Herg['D4000min'].iloc[indiceHERG])&
                                (noAGNsSinSelec['D4000']<=
                                Herg['D4000max'].iloc[indiceHERG]))]
            #SELECCIONO LOS INDICES DE LA MASCARA ANTERIOR
            Var=Var.index
            #ELIJO ALEATORIAMENTE UNA DE LAS GALAXIAS ACTIVAS SIN SELECCIONAR    
            indiceNoAGN=np.random.choice(Var)
            #SE MARACA A LA GALAXIA NO ACTIVA SELECCIONADA
            noAGNs['EsControl_oNo'].iloc[indiceNoAGN]=1
            #SE INDICA CUAL ES LA GALAXIA HERG DE LA QUE ES CONTROL LA GALAXIA NO ACTIVA
            noAGNs['HERGdeLaQueEsControl'].iloc[indiceNoAGN]=indiceHERG
    
        fun1(Herg['MPA_IDX'])
        
        #Se filtra y eliminan columnas             
        noAGNs=noAGNs[noAGNs['EsControl_oNo']==1]


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


