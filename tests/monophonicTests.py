'''tests for the monophonic implementation. Checking for valid betas etc.'''
import pytest
from .src.barbanchoInharmonicity import computeBarbancho
import numpy as np
def checkBetas(name):
    fft = np.load(name)
    assert 10**(-7) < computeBarbancho(fftAmps, fftFreqs, sr) < 10**(-2)