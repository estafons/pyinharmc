from barbanchoInharmonicity import *
from itertools import combinations
import copy

def elegibilityCheck(f0l, pOrder, winSize, betas) -> list:
    f0bs = zip(f0l, betas)
    fks = [pOrder*f0*sqrt(1 + beta*pOrder**2) for (f0, beta) in f0bs]
    eligibles = zip(copy.deepcopy(fks), f0l)
    lims = [dLimits(fk, winSize) for fk in fks]
    combs = combinations(lims, 2)
    #pFreq = sqrt(1 + beta*pOrder**2)*pOrder*f0
    for comb in combs:
        x, y = comb
        if max(x[0], y[0]) < min(x[1], y[1]) + 1: # intersection
            eligibles.remove(x[2])
            eligibles.remove(y[2])
    return eligibles

def multiPartialTrack(fft, f0l: list, N: int, winSize: float) -> list:
    """needs work. Not robust at all, not implemented correctly 
    for the beta and how it passes through to every iteration"""
    differences = {}
    partials = {}
    betas = {}
    for f0 in f0l:
        differences[f0] = []
        partials[f0] = []
        betas[f0] = []
    # if initB is not None:
    #     betas = copy.deepcopy(initB)
    # else:
    #     betas = []
    for pOrder in range(1, N + 1):
        eligibles = elegibilityCheck(f0l, pOrder, winSize, betas)
        for eligible in eligibles:
            pFreq, f0 = eligible
            rStart, rEnd, _ = dLimits(pFreq, winSize)
            partialF = partialDetect(fft, rStart, rEnd)
            partials[f0].append(partialF)
            differences[f0].append(computeDiff(partialF, f0, pOrder))
            beta = innerComputeBarbanchoInharm(differences, f0)
            betas[f0].append(beta)
    return partials, differences, betas