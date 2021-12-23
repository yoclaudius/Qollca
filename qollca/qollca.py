import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from astropy.cosmology import LambdaCDM

def sigma_5():
    pass

class Muestra_control():

    def m_c():
        pass

    def rand_selec(mp, ms, nc1, nc2, nc3, delta1, delta2, delta3):
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
        from astropy.table import Table
        import qollca
        dir = '/usr/local/datos/Tesis de grado/Datos_iniciales/'
        nom = 'HERG.fits'
        mainsample = Table.read(dir + nom)
        nom = 'noAGNs.fits'
        secondarysample = Table.read(dir + nom)
        mainsample = mainsample.to_pandas()
        secondarysample = secondarysample.to_pandas()
        controls, nocontrols = qollca.Muestra_control.rand_selec(mainsample, 
                                               secondarysample, 'SMASS_MEDIAN',
                                              'ABSR', 'Z', 'D4000', 0.15, 0.2,
                                               0.02, 0.05)
        import matplotlib.pyplot as plt
        import numpy as np
        Columna1 = mainsample['SMASS_MEDIAN']
        Columna2 = controls['SMASS_MEDIAN']
        Columna1 = Columna1[Columna1<Columna1.quantile(q=0.99)]
        Columna1 = Columna1[Columna1>Columna1.quantile(q=0.01)]
        Columna2 = Columna2[Columna2<Columna2.quantile(q=0.99)]
        Columna2 = Columna2[Columna2>Columna2.quantile(q=0.01)]
        pesos1 = np.ones_like(Columna1)/float(len(Columna1))
        pesos2 = np.ones_like(Columna2)/float(len(Columna2))
        plt.figure(figsize=(1.4*6.4,1.4*4.8),constrained_layout=True)
        CantidadBines = 15
        bines = np.linspace(Columna1.min(),Columna1.max(),CantidadBines)
        plt.hist(Columna1,bins=bines,weights=pesos1, color='b',histtype='step',
                lw=3,alpha=0.75,label='HERGs')
        plt.hist(Columna2,bins=bines,weights=pesos2, color='indigo',histtype='step',
                linestyle=':',lw=4,alpha=1,label='ControlHERGs')
        plt.legend(loc='upper right',fontsize=16,frameon=False)
        plt.xlabel(r'$log_{10}(M_{\ast} / M_{\odot})$',fontsize=16)
        plt.tick_params(labelsize=16)
        plt.show()
        
        import pandas as pd
        controlstot = controls.copy()
        nnewcon = 100
        nnew = 101
        while nnew >= nnewcon:
            cond1 = len(nocontrols)
            controls, nocontrols = qollca.Muestra_control.rand_selec(mainsample, 
                                    nocontrols, 
                                    'SMASS_MEDIAN',
                                    'ABSR', 'Z', 'D4000', 0.15, 0.2,
                                                   0.02, 0.05)
            controlstot = pd.concat([controlstot, controls])
            cond2 = len(nocontrols)
            nnew = cond1 - cond2
            print('nnew = ', nnew)
        
        import matplotlib.pyplot as plt
        import numpy as np
        Columna1 = mainsample['SMASS_MEDIAN']
        Columna2 = controlstot['SMASS_MEDIAN']
        Columna1 = Columna1[Columna1<Columna1.quantile(q=0.99)]
        Columna1 = Columna1[Columna1>Columna1.quantile(q=0.01)]
        Columna2 = Columna2[Columna2<Columna2.quantile(q=0.99)]
        Columna2 = Columna2[Columna2>Columna2.quantile(q=0.01)]
        pesos1 = np.ones_like(Columna1)/float(len(Columna1))
        pesos2 = np.ones_like(Columna2)/float(len(Columna2))
        plt.figure(figsize=(1.4*6.4,1.4*4.8),constrained_layout=True)
        CantidadBines = 15
        bines = np.linspace(Columna1.min(),Columna1.max(),CantidadBines)
        plt.hist(Columna1,bins=bines,weights=pesos1, color='b',histtype='step',
                lw=3,alpha=0.75,label='HERGs')
        plt.hist(Columna2,bins=bines,weights=pesos2, color='indigo',histtype='step',
                linestyle=':',lw=4,alpha=1,label='ControlHERGs')
        plt.legend(loc='upper right',fontsize=16,frameon=False)
        plt.xlabel(r'$log_{10}(M_{\ast} / M_{\odot})$',fontsize=16)
        plt.tick_params(labelsize=16)
        plt.show()
        
        """        
        
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


def luminosity_function(ms, n_reddening, n_redshift):
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
    from astropy.table import Table
    import qollca
    dir = '/usr/local/datos/Tesis de grado/llllllllllllllll/'
    nom = 'ms.fits'
    mainsample = Table.read(dir + nom)
    mainsample = mainsample[mainsample['dered_g'] > 0]

    qollca.luminosity_function(mainsample, 'dered_r', 'z')

    import matplotlib.pyplot as plt
    import numpy as np
    frec, lim = np.histogram(mainsample['Mr'], weights=mainsample['1/Vmax']
                             ,bins=30)
    bin_centers = 0.5 * (lim[1:] + lim[:-1])
    ancho_bin = lim[2] - lim[1]
    frec = frec / ancho_bin
    plt.plot(bin_centers, frec, lw=1, ls='--', marker='o', alpha=0.5)
    plt.gca().invert_xaxis()
    plt.yscale('logit')

    """
    H0 = 100.    # Hubble constant
    Omm = 0.3     # Omega matter
    Oml = 0.7     # Omega lambda

    cosmo = LambdaCDM(H0=H0, Om0=Omm, Ode0=Oml)

    m_lim = 17.7
    m_limmin = 14.5
    z_min = 0.03

    # Se calcula la distancia de luminosidad (Mpc)
    ms['DL'] = cosmo.luminosity_distance(ms[n_redshift]).value
    # Se calcula el modulo de distancia
    ms['DM'] = 5 * np.log10(ms['DL'])
    # Se calcula la magnitud absoluta
    ms['Mr'] = ms[n_reddening] - ms['DM'] - 25
    # Se obtiene la distancia maxima a la cual es posible observar cada galaxia    
    ms['Dmax'] = np.zeros(len(ms))
    for i in range(len(ms)):
        print('Objeto trabajado en dist max: ', i)
        m = 0
        #z = ms[n_redshift].iloc[i]
        z = ms[n_redshift][i]
        while m_lim > m:
            # Calculo modulo de distancia
            DL = cosmo.luminosity_distance(z).value
            DM = 5 * np.log10(DL)
            # Se calcula la magnitud aparente
            #m = ms['Mr'].iloc[i] + 25 + DM
            m = ms['Mr'][i] + 25 + DM
            z = z + 0.01
        #ms.loc[i,'Dmin'] = DL
        ms['Dmax'][i] = DL
    # Se obtiene la distancia minima a la cual es posible observar cada galaxia
    ms['Dmin'] = np.zeros(len(ms))
    for i in range(len(ms)):
        print('Objeto trabajado en dist min: ', i)
        m = 99
        #z = ms[n_redshift].iloc[i]
        z = ms[n_redshift][i]
        while m_limmin < m and z > z_min:
            # Calculo modulo de distancia
            DL = cosmo.luminosity_distance(z).value
            DM = 5 * np.log10(DL)
            # Se calcula la magnitud aparente
            #m = ms['Mr'].iloc[i] + 25 + DM
            m = ms['Mr'][i] + 25 + DM
            z = z - 0.01
        #ms.loc[i,'Dmin'] = DL
        ms['Dmin'][i] = DL
    # Se obtiene el Vmax
    ms['Vmax'] = ((ms['Dmax']**3 - ms['Dmin']**3) * 2.23226371519234) / 3
    # ms['Vmax'] = ((ms['Dmax']**3) * 2.23226371519234) / 3
    # Se obtienen los pesos
    ms['1/Vmax'] = 1 / ms['Vmax']


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


def muestra_limpia():
    pass
