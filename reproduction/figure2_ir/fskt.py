import glob
import tqdm
import gsd.hoomd
import numpy as np
import memspectrum
import h5py
from scipy import fftpack
from scipy.fftpack import fft, ifft, ifftshift
from itertools import islice
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set up matplotlib to use LaTeX and classic style
plt.style.use('classic')
mpl.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsmath}',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 12
})

def smooth(x,window_len=11,window='hamming'):
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat': 
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

def acf(x):
    finalval = 0
    xp = ifftshift(x - np.average(x))
    n, = xp.shape
    xp = np.r_[xp[:n//2], np.zeros_like(xp), xp[n//2:]]
    f = fft(xp)
    p = np.absolute(f)**2
    pi = ifft(p)
    finalval += np.real(pi)[:n//2]/(np.arange(n//2)[::-1]+n//2)
    return np.real(finalval)

def fft3(x, timestep):
    lineshape = fftpack.dct(x, type=1)
    freq_au = np.linspace(0, 0.5/timestep * 1e15, len(x))
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    field_description = freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra

def fft3_maxent(x, timestep):
    m = memspectrum.MESA()
    m.solve(x, method = 'Standard', optimisation_method = 'VM')
    
    freq_au, lineshape = m.spectrum(dt=timestep, onesided=True)
    freq_au = freq_au*1e15
    
    autocorr = m.compute_autocorrelation(dt=timestep, normalize=True)
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    field_description = freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra, autocorr

# Create figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
fig.subplots_adjust(wspace=0.3)

# Define line styles and colors for different files
styles = {
        './**/simu_1.xc.xyz.fkt.txt': {'color': 'blue', 'linestyle': '-', 'label': 'LAMMPS+i-PI', 'lw': 2},
        'isf_hoomd_*.txt': {'color': 'green', 'linestyle': '-', 'label': 'HOOMD (GPU)', 'lw': 2},
        'newisf_hoomd_*.txt': {'color': 'red', 'linestyle': '-', 'label': 'HOOMD (CPU)', 'lw': 2}
}

#for filenamee in ['./**/dipole_ipi.txt', 'dipoles_hoomd_*.txt']:
for filenamee in ['./**/simu_1.xc.xyz.fkt.txt','isf_hoomd_*.txt','newisf_hoomd_*.txt']:
    filenames = glob.glob(filenamee,recursive=True)[:1]
    dacf = 0
    dacfs = 0
    time = 0
    count = 0
    for filename in filenames:
        try:
            data = np.loadtxt(filename)
            #if 'hoomd' in filename or 'ipi' in filename:
            if 'hoomd' in filenamee:
                dacf += data[:,0]
                dacfs += data[:,1]
            else:
                dacf += data[:,1]
                dacfs += data[:,2]
        except Exception as e:
            print(e)
        time = np.arange(0, len(dacf)) * 4/1000.0 *250
    # Plot DACF
    ax1.plot(time, dacf[:len(dacf)]/dacf[0], **styles[filenamee])
    ax2.plot(time, dacfs[:len(dacfs)]/dacfs[0], **styles[filenamee])

# Configure IR spectrum plot
ax1.set_xlabel(r'Time (ps)')
ax1.set_ylabel(r'$F(k,t)$')
ax1.set_xlim([0, 100])
ax1.set_ylim(top=1,bottom=-0.01)
ax1.legend()

# Configure DACF plot
ax2.set_xlabel(r'Time (ps)')
ax2.set_ylabel(r'$F_s(k,t)$')
ax2.set_xlim([0, 100])
ax2.set_ylim(top=1,bottom=-0.01)
ax2.legend()

# Adjust layout to prevent label clipping
plt.tight_layout()
plt.show()
