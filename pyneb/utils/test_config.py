import sys

print('Python version:')
print(sys.version)

try:
    import numpy
    print('Numpy: {0}'.format(numpy.__version__))
except:
    print('NO Numpy')

try:
    import matplotlib
    print('Matplotlib: {0}'.format(matplotlib.__version__))
except:
    print('NO matplotlib')

try:
    import scipy
    print('scipy: {0}'.format(scipy.__version__))
except:
    print('NO scipy')

try:
    import pyCloudy
    print('pyCloudy: {0}'.format(pyCloudy.__version__))
except:
    print('NO pyCloudy')

try:
    import pyneb
    print('pyneb: {0}'.format(pyneb.__version__))
except:
    print('NO pyneb')
    
try:
    import atpy
    print('atpy: {0}'.format(atpy.__version__))
except:
    print('NO atpy')

try:
    import asciitable
    print('asciitable: {0}'.format(asciitable.__version__))
except:
    print('NO asciitable')

try:
    import pyfits
    print('pyfits: {0}'.format(pyfits.__version__))
except:
    try:
        import astropy.io.fits as pyfits
        print('pyfits: from astropy')
    except:
        pn.log_.error('pyfits not installed')