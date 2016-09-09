import os, sys
import pyneb as pn
import numpy as np
from pyneb.utils.misc import int_to_roman, roman_to_int
try:
    import pyfits
except:
    try:
        import astropy.io.fits as pyfits
    except:
        pn.log_.error('pyfits not installed')
from physics import sym2name, gsFromAtom

    
def _updateHeader(hdu, input_file, content, elem, spec):
    """
    Update the header of the output fits file

    Parameters: 
        - hdu          pyfits file
        - input_file   name of the ascii files containing the original data
        - content      type of physical quantities contained in the file
        - elem         symbol of the selected element
        - spec         ionization stage in spectroscopic notation (I = 1, II = 2, etc.)
        
    """
    hdu.header.update('ORIGIN', 'File {0}, PyNeb v{1}'.format(input_file, pn.__version__))
    hdu.header.update('CONTENT', content)
    hdu.header.update('ATOM', sym2name[elem.capitalize()])
    hdu.header.update('SPECTRUM', roman_to_int(spec))
    atom = elem.capitalize() + str(roman_to_int(spec))
    gs = gsFromAtom(atom)
    hdu.header.update('GSCONFIG', gs)
    

def writeAtom(input_file):
    """ 
    Write the Es and As read from an ascii file into a fits file in PyNeb format.
    The bibliographic references are read from comment lines in the file.
        
    Parameters:
        - input_file    name of the ascii input file containing the energy levels (in 1/A or eV) and the Aij,
                        in the format elem_spec_atom_REF-STRING.dat

    Example: 
    writeAtom('o_iii_atom_AB01-CD02-EFG03.dat')

    """
    if not os.path.exists(input_file):
        pn.log_.error('{0} not found.'.format(input_file), calling='writeAtom')
        sys.exit(0)
    data_type = input_file.split(".")[0].split("_")[2]
    if (data_type != 'atom'):
        pn.log_.error('The name of {0} does not conform to the prescribed format.'.format(input_file), calling='writeAtom')
        sys.exit(0)
    elem = input_file.split(".")[0].split("_")[0]
    spec = input_file.split(".")[0].split("_")[1]
    ref = input_file.split(".")[0].split("_")[-1]
    out_fits = '{0}_{1}_{2}_{3}.fits'.format(elem.lower(), spec.lower(), data_type, ref)
    atom_dat = []

    # Read data from ascii file
    #
    # Read energies and stat weights 
    f = open(input_file)
    units = f.readlines()[1].split()[0]
    all_data = np.genfromtxt(input_file, comments='***', skip_header=2, unpack=True)
    # Read Es
    energy = all_data[0, :]
    NLevels = len(energy)
    if units[0] is 'eV': 
        energy = energy / (pn.CST.RYD_EV * pn.CST.RYD_ANG)
    # Read statistical weights
    stat_weight = all_data[1, :]
    atom_dat.append(pyfits.Column(name='Energy', format='1D', unit='1/Ang', array=energy))
    atom_dat.append(pyfits.Column(name='Stat_Weight', format='1J', unit='none', array=stat_weight))    
    # Read As
    for j in np.arange(NLevels - 1):
        lev_j = j + 1
        fieldName = 'A(*->{0})'.format(lev_j)
        atom_dat.append(pyfits.Column(name=fieldName, format='1D', unit='1/sec', array=all_data[1 + lev_j, :]))

    # Create fits file and start adding info to the header
    atom_hdu = pyfits.new_table(atom_dat)   
            
    # Primary header (extension 0)
    hdu = pyfits.PrimaryHDU()
    hdu.header.update('FILENAME', out_fits)
    hdu.header.update('NEXTEND', 1)
 
    # Write references in the header
    refs = []
    for line in open(input_file):
        if '***' in line:
            text = line.strip("***").split("'")
            if ('SOURCE' in text[0]) or ('NOTE' in text[0]):
                refs.append((text[0], text[1]))

    _updateHeader(atom_hdu, input_file, 'Energy levels and transition probabilities', elem, spec)
    _updateHeader(hdu, input_file, 'Energy levels and transition probabilities', elem, spec)
    atom_hdu.header.update('N_LEVELS', NLevels)
    for item in refs:
        atom_hdu.header.update(item[0], item[1])

    hdulist = pyfits.HDUList([hdu])
    hdulist.append(atom_hdu)
    hdulist.writeto(out_fits, clobber=True)


def writeColl(input_file, temp_in_cols = True, chebOrder = None):
    """ 
    Write the Omegas read from an ascii file into a fits file in PyNeb format.
    The bibliographic references are read from comment lines in the file.
        
    Parameters:
        - input_file    name of the ascii input file containing the Omegas,
                        in the format elem_spec_coll_REF-STRING.dat
        - temp_in_cols  If True, each column of the input file corresponds to a given Te, 
                         otherwise each rows corresponds to a given Te (a transposition is 
                         thus applied to the data)
    Example: 
    writeColl('o_iii_coll_AB01-CD02-EFG03.dat')

    """
    def _omegaArray(lev_j, lev_i):
        """
        Return the array of collision strength values for given levels
 
        """       
        ij = (cs[:, 0] == lev_j) & (cs[:, 1] == lev_i)
        if ij.sum() == 0.:
            return np.zeros(n_temp)
        else:
            return cs[ij, 2::].reshape(n_temp)

    if not os.path.exists(input_file):
        pn.log_.error('{0} not found.'.format(input_file), calling='writeColl')
        sys.exit(0)
    data_type = input_file.split(".")[0].split("_")[2]
    if (data_type != 'coll'):
        pn.log_.error('The name of {0} does not conform to the prescribed format.'.format(input_file), calling='writeColl')
        sys.exit(0)
    elem = input_file.split(".")[0].split("_")[0]
    spec = input_file.split(".")[0].split("_")[1]
    ref = input_file.split(".")[0].split("_")[-1]
    out_fits = '{0}_{1}_{2}_{3}.fits'.format(elem.lower(), spec.lower(), data_type, ref)
    coll_dat = []

    # Read Omegas
    cs = np.genfromtxt(input_file, comments='***')
    if not temp_in_cols:
        cs = cs.transpose()
    n_temp = np.size(cs[0]) - 2
    temArray = cs[0, 2::]
    NLevels = int(np.max(cs[:, 1]))
    if chebOrder is None:
        chebOrder = n_temp
    coll_dat.append(pyfits.Column(name='Te', format='1E', unit='log(K)', array=temArray)) # Change if linear T used (see below)
    for j in np.arange(NLevels - 1):
        lev_j = j + 1
        lev_i = lev_j + 1
        while (lev_i < NLevels + 1):
            fieldName = 'Omega({0}->{1})'.format(int(lev_i), int(lev_j))
            omegaArray = _omegaArray(lev_j, lev_i)
            coll_dat.append(pyfits.Column(name=fieldName, format='1E', unit=chebOrder, array=omegaArray))
            lev_i += 1
    coll_hdu = pyfits.new_table(coll_dat)        
            
    # Primary header (extension 0)
    hdu = pyfits.PrimaryHDU()
    hdu.header.update('FILENAME', out_fits)
    hdu.header.update('NEXTEND', 1)
 
    # Write references in the header
    refs = []
    for line in open(input_file):
        if '***' in line:
            text = line.strip("***").split("'")
            if ('SOURCE' in text[0]) or ('NOTE' in text[0]):
                refs.append((text[0], text[1]))

    _updateHeader(coll_hdu, input_file, 'Collision strengths', elem, spec)
    _updateHeader(hdu, input_file, 'Collision strengths', elem, spec)
    coll_hdu.header.update('N_LEVELS', NLevels)
    for item in refs:
        coll_hdu.header.update(item[0], item[1])

    hdulist = pyfits.HDUList([hdu])
    hdulist.append(coll_hdu)
    hdulist.writeto(out_fits, clobber=True)


class Hdr(object):
    """
    Explore and update headers of fits files
    
    Parameters:
        - fitsFile  name of the fits file
        - ext       file extension
    
    Examples:
        f=hdr('o_iii.fits', 2)
        f.show()
        f.up('ATOM', 'kriptonyte')
        f.up('NATURE', 'invented')
        f.save
        
    """    
    def __init__(self, fitsFile, ext):
        self.fitsFile = fitsFile
        self.hdu = pyfits.open(fitsFile, ignore_missing_end=True)
        self.ext = ext
    

    def show(self):
        print self.hdu[self.ext].header
    

    def up(self, key, value):
        """ This will be changed to header.set at some point"""
        self.hdu[self.ext].header.update(key, value)
        

    def delete(self, key):
        del self.hdu[self.ext].header[key]


    def cheb(self, orderFile):    
        """ 
        Update the Chebyshev order of a coll fits file.
        
        Parameters:
            - orderFile    ascii file with upper level, lower level, order
        
        """    
        j_order, i_order, order = np.genfromtxt(orderFile, comments='***', usecols=(0, 1, 2), unpack=True)    
        NLevels = np.max(j_order)
        order_array = np.zeros([NLevels, NLevels])
        for index in range(len(j_order)):
            order_array[int(j_order[index] - 1), int(i_order[index] - 1)] = int(order[index])
            
        # Iterate over the keys to associate each unit with its Omega(*->*) field
        for key in self.hdu[self.ext].header.keys():
            if 'TTYPE' in key:
                num = key.strip('TTYPE')
                unit = 'TUNIT' + num
                field_name = self.hdu[self.ext].header[key]
                # The try is to skip the Te field
                try:
                    lev_j = field_name.split('->')[0].split('Omega(')[1]
                    lev_i = field_name.split('->')[1].split(')')[0]
                    j = int(lev_j) - 1
                    i = int(lev_i) - 1
                    self.up(unit, order_array[j, i])
                except:
                    pass


    def save(self):
        self.hdu.writeto(self.fitsFile, clobber=True)

# Test
# writeAtom('o_ii_atom_TEST.dat')
# writeColl('o_ii_coll_TEST.dat')
