import qollca
import numpy as np
import numpy.testing as npt

# ============================================================================
# CONSTANTS
# ============================================================================





# ============================================================================
# FIXTURES
# ============================================================================





# ============================================================================
# TESTS
# ============================================================================


def test_monte_carlo():
    # Se crea una distribución normal
    al = np.random.RandomState(30)
    datos = al.normal(2,3,30000)
    
    histograma = qollca.monte_carlo(datos,bines=20, rango=(-1,1), test=True)    
    npt.assert_almost_equal(histograma, datos)
    
