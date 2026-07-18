from numba import njit
import tqdm
import gsd.hoomd
import numpy as np
import memspectrum
import h5py
from scipy import fftpack
from scipy.fftpack import fft, ifft, ifftshift
from itertools import islice
import matplotlib.pyplot as plt
import sys
import os
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Compute IR spectrum from trajectory data')
parser.add_argument('--idx', type=int, required=True, help='Index number for the trajectory')
parser.add_argument('--dir', type=str, required=True, help='Directory containing the trajectory files')

args = parser.parse_args()

# Change to the specified directory
os.chdir(args.dir)
print(f"Working in directory: {os.getcwd()}")

def read_lammps_data(filename):
    """Extracts box dimensions, atom positions, charges, and image flags from a LAMMPS data file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    atoms_section = False
    box_bounds = {}
    charges = []
    positions = []
    images = []
    seeblank = False
    for i, line in enumerate(lines):
        words = line.split()
        if len(words) == 4 and words[2] in ["xlo", "ylo", "zlo"]:
            # Box boundary lines format: xlo xhi, ylo yhi, zlo zhi
            axis = words[2][0]  # 'x', 'y', or 'z'
            box_bounds[axis] = float(words[1]) - float(words[0])  # Compute box length

        if "Atoms" in line:
            atoms_section = True
            continue
        if atoms_section and line.strip() == "" and not seeblank:
            seeblank = True
            continue
        if atoms_section and seeblank:
            data = line.split()
            if len(data) < 10:
                continue  # Ensure the line has enough columns

            atom_id = int(data[0])   # Atom ID
            charge = float(data[3])  # Assuming charge is in 4th column
            x, y, z = map(float, data[4:7])  # Positions
            ix, iy, iz = map(int, data[-3:])  # Image flags
            charges.append(charge)
            positions.append([x, y, z])
            images.append([ix, iy, iz])

    box_lengths = np.array([box_bounds.get(axis, 0) for axis in "xyz"])  # Ensure x, y, z order
    return np.array(charges), np.array(positions), np.array(images), box_lengths

def compute_dipole_moment(charges, positions):
    """Computes dipole moment from charges and positions."""
    return np.sum(charges[:, None] * positions, axis=0)

def unwrap_positions(positions, images, box_lengths):
    """
    Unwrap particle positions across periodic boundaries.
    
    Args:
        positions: Array of wrapped positions (N x 3)
        images: Array of image flags (N x 3)
        box_lengths: Array of box dimensions (3,)
        
    Returns:
        Array of unwrapped positions (N x 3)
    """
    # Convert inputs to cupy arrays if they aren't already
    pos = np.asarray(positions)
    img = np.asarray(images)
    box = np.asarray(box_lengths)
    
    # Unwrap by adding box lengths multiplied by image flags
    return pos + img * box[None, :]

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

def acf(x):
    finalval = 0
    xp = ifftshift(x)#(x - np.average(x))/np.std(x))
    n, = xp.shape
    xp = np.r_[xp[:n//2], np.zeros_like(xp), xp[n//2:]]
    f = fft(xp)
    p = np.absolute(f)**2
    pi = ifft(p)
    finalval += np.real(pi)[:n//2]/(np.arange(n//2)[::-1]+n//2)
    return np.real(finalval)

def fft3(x, timestep):
    # Adding zeros to the end of x
    lineshape = fftpack.dct(x, type=1)
    freq_au = np.linspace(0, 0.5/timestep * 1e15, len(x))
    # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    # Calculate spectra
    field_description =  freq_au**2
    #field_description =  freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra
    #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

def fft3_maxent(x, timestep):
    """# Adding zeros to the end of x
    lineshape = fftpack.dct(x, type=1)
    freq_au = np.linspace(0, 0.5/timestep * 1e15, len(x))
    # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    # Calculate spectra
    field_description =  freq_au**2
    #field_description =  freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra
    """
    m = memspectrum.MESA()
    m.solve(x, method = 'Standard', optimisation_method = 'FPE')#,regularisation=4.0)
    
    freq_au, lineshape = m.spectrum(dt =timestep,onesided=True)
    freq_au = freq_au*1e15
    
    autocorr = m.compute_autocorrelation(dt = timestep,normalize=True)
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    # Calculate spectra
    field_description =  freq_au**2
    #field_description =  freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra #, autocorr
    #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

@njit
def compute_acf_numba(time_series):
    n = len(time_series)
    mean = np.mean(time_series)
    autocorr = np.zeros(n)

    for tau in range(n):
        sum_product = 0.0
        for t in range(n - tau):
            sum_product += (time_series[t] ) * (time_series[t + tau] )
        autocorr[tau] = sum_product / (n - tau)

    # Normalize by R(0)
    autocorr /= autocorr[0]
    return autocorr

freq_hmd = []
spectra_hmd = []

with h5py.File(f'ir-{args.idx}.h5', 'r') as f:
    # Get dipole moments directly from the file
    dgx = f['/hoomd-data/cavitymd/DipoleMoment/dgx'][:]
    dgy = f['/hoomd-data/cavitymd/DipoleMoment/dgy'][:]
    
    # Get timestep information for later use if needed
    timesteps = f['/hoomd-data/Simulation/timestep'][:]
    dt = f['/hoomd-data/Status/Dt'][0]  # assuming dt is constant
    
    print(f"Processing {len(dgx)} frames with timestep {dt} ps")
    
    # Compute autocorrelation for x and y components
    dacf = compute_acf_numba(dgx)
    dacf += compute_acf_numba(dgy)
    
    # Save the normalized autocorrelation function
    np.savetxt(f'dacf_hoomd_{args.idx}.txt', np.vstack((dacf/dacf[0],)).T)
    
    # Optionally compute and save the spectrum
    freq_cminverse, spectrax = fft3_maxent(dgx, 4)
    freq_cminverse, spectray = fft3_maxent(dgy, 4)
    np.savetxt(f'spectrum_hoomd_{args.idx}.txt', 
               np.vstack((freq_cminverse, spectrax+spectray)).T,
               header='Frequency(cm^-1) Intensity')
