"""here the inharmonicity related function lye"""

from typing import List
import numpy as np
from scipy.optimize import least_squares
from math import sqrt

class Fft:
    def __init__(self, fftAmps, fftFreqs, sr):
        self.amps = fftAmps
        self.freqs = fftFreqs
        self.sr = sr
        self.sz = len(fftAmps)

def inharmonicomp(i: int) -> int:
    """function returning computation of inharmonicity"""
    i = i + 2
    return i
def innerComputeBarbanchoInharm(differences: list, f0: float) -> float:
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

def partialTrack(fft, f0: float, N: int, winSize: float,
                initB = 1) -> list:
    """Function that tracks the partial in the given note instace.
    Works by invoking partialDetect and innerComputeBarbanchoInharm for adjusting 
    partialDetect window according to inharmonicity estimation so far"""
    def computeDiff(pFreq, f0, pOrder):
        """compute and return a tuple with the difference between
        measured and theoretical partial and the associated order"""
        return (pFreq - f0*pOrder, pOrder)
    def dLimits(pFreq: float, winSize: float) -> tuple:
        """function computing the limits used for partial detection."""
        rStart = pFreq - winSize/2
        rEnd = pFreq + winSize/2
        return rStart, rEnd

    partials, differences = [], []
    beta = initB
    for pOrder in range(1, N + 1):
        pFreq = sqrt(1 + beta*pOrder**2)*pOrder*f0
        rStart, rEnd = dLimits(pFreq, winSize)
        partials.append(partialDetect(fft, rStart, rEnd))
        differences.append(computeDiff(pFreq, f0, pOrder))
        beta = innerComputeBarbanchoInharm(differences, f0)
    return partials, differences


def partialDetect(fft, rStart: float, rEnd: float) -> float:
    """function that detects the highest peak in this case,
    used for detection of the partial in that range.
    fft: The fft of the note instance under inspection
    rStart: Starting frequency range to detect partial
    rEnd: Ending frequency range to detect partial"""
    def highPeak(fft, rStart: float, rEnd: float) -> float:
        def castToBin(sr, sz, freq):
            return sr*freq/sz
        rBStart = castToBin(fft.sr,fft.sz, rStart)
        rBEnd = castToBin(fft.sr, fft.sz, rEnd)
        peakAmp = 0
        for t in zip(fft.amps[rBStart:rBEnd], fft.freqs[rBStart:rBEnd]):
            if t[0] > peakAmp:
                peakAmp, peakFreq = t
        return peakFreq
    return highPeak(fft, rStart, rEnd)


def computeBarbanchoInharm(fftAmps, fftFreqs, sr, f0, winSize, N, *kwargs):
    """function that performs partial tracking and returns associated
    inharmonicity coefficient as described in paper.
    fftAmps: is an array of amplitudes of the fft
    fftFreqs: an array of frequencies associated to the fft
    sr: the sampling rate
    f0: the estimated fundamental frequency
    winSize: the partial tracking window size
    N: the number of partials to account for
        additional possible arguements are:
        initB: an estimation of the inharmonicity coefficient in order to 
        refine expected results. If not familar with concepts, don't use.
    """
    fft = Fft(fftAmps=fftAmps, fftFreqs=fftFreqs, sr=sr)
    _, differences = partialTrack(fft, f0, N, winSize, *kwargs)
    return innerComputeBarbanchoInharm(differences, f0)