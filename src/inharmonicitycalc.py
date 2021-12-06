"""here the inharmonicity related function lye"""

from typing import List
import numpy as np
from scipy.optimize import least_squares


class fft:
    def __init__(self, fftAmps, fftFreqs, sr):
        self.amps = fftAmps
        self.freqs = fftFreqs
        self.sr = sr
        self.sz = len(fftAmps)

def inharmonicomp(i: int) -> int:
    """function returning computation of inharmonicity"""
    i = i + 2
    return i

def partialTrack(fft: tuple, f0: float) -> list:
    """Function that tracks the partial in the given note instace.
    Works by invoking partialDetect and computeInharm for adjusting 
    partialDetect window according to inharmonicity estimation so far"""
    pass

def castToBin(sr, sz, freq):
    return sr*freq/sz

def partialDetect(fft, rStart: float, rEnd: float) -> float:
    """function that detects the highest peak in this case,
    used for detection of the partial in that range.
    fft: The fft of the note instance under inspection
    rStart: Starting frequency range to detect partial
    rEnd: Ending frequency range to detect partial"""
    def highPeak(fft, rStart: float, rEnd: float) -> float:
        rBStart = castToBin(fft.sr,fft.sz, rStart)
        rBEnd = castToBin(fft.sr, fft.sz, rEnd)
        peakAmp = 0
        for t in zip(fft.amps[rBStart:rBEnd], fft.freqs[rBStart:rBEnd]):
            if t[0] > peakAmp:
                peakAmp, peakFreq = t
        return peakFreq
    return highPeak(fft, rStart, rEnd)

def computeInharm(differences: list, f0: float) -> float:
    """Function that computes the inharmonicity as proposed
    in Barbancho et al with least squares.
    differences: The difference between measured partial and expected
    f0: The fundamental frequency of the note under inspection."""
    def compute_least(u,y):
        def model(x, u):
            return x[0] * u**3 + x[1]*u + x[2]
        def fun(x, u, y):
            return model(x, u)-y
        def jac(x, u, y):
            J = np.empty((u.size, x.size))
            J[:, 0] = u**3
            J[:, 1] = u
            J[:, 2] = 1
            return J
        x0=[0.00001,0.00001,0.000001]
        res = least_squares(fun, x0, jac=jac,bounds=(0,np.inf), 
                            args=(u, y),loss = 'soft_l1', verbose=0)
        return res.x

    differences, orders = zip(*differences)
    u=np.array(orders)+2
    [a,b,_] = compute_least(u,differences) # least_squares
    beta=2*a/(f0+b) # Barbancho et al. (17)
    return beta
