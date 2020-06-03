import numpy as np
from scipy import linalg
from math import exp, sqrt


def SimilarProd(peaksA, peaksB, sgm):
    # Similarity function for two sets of peaks A and B.
    if len(peaksA) == 0 and len(peaksB) == 0:
        return 1
    elif len(peaksA) == 0 or len(peaksB) == 0:
        return 0
    Norm = 2.*len(peaksA)*len(peaksB)/(len(peaksA)**2 + len(peaksB)**2)
    # Norm = Norm**5
    # Distances from all points in peaksA to closest in peaksB.
    AtoB = []
    dist = 0
    for ptA in peaksA:
        dist = abs(ptA - peaksB[0])
        for ptB in peaksB:
            tmp = abs(ptB - ptA)
            if tmp < dist:
                dist = tmp
        AtoB.append(dist)
    # Distances from all points in peaksA to closest in peaksB.
    BtoA = []
    dist = 0
    for ptB in peaksB:
        dist = abs(ptB - peaksA[0])
        for ptA in peaksA:
            tmp = abs(ptB - ptA)
            if tmp < dist:
                dist = tmp
        BtoA.append(dist)
    # Get Gaussian weights.
    # Form of the weight function is
    # Result = Norm*Prod_i exp(-(AtoB_i)^2/2sgm^2) *
    #         Prod_i exp(-(BtoA_i)^2/2sgm^2)
    #    So for two identical sets A and B Result = 1, but will
    #    rapidly decay for different sets.

    Result = 1
    sgm2 = sgm**2
    for elem in AtoB:
        Result *= exp(-0.5*elem**2/sgm2)
    for elem in BtoA:
        Result *= exp(-0.5*elem**2/sgm2)
    Result *= Norm
    return Result



def SimilarSum(peaksA, peaksB, sgm):
    # Similarity function for two sets of peaks A and B.
    if len(peaksA) == 0 and len(peaksB) == 0:
        return 1
    elif len(peaksA) == 0 or len(peaksB) == 0:
        return 0
    elif abs(len(peaksA)-len(peaksB))>2:
        return 0

    Norm = 2.*len(peaksA)*len(peaksB)/(len(peaksA)**2 + len(peaksB)**2)
    Norm = Norm**4
    # Distances from all points in peaksA to closest in peaksB.
    AtoB = []
    dist = 0
    J = 0
    for i, ptA in enumerate(peaksA):
        J = 0
        dist = abs(ptA - peaksB[0])
        for j, ptB in enumerate(peaksB):
            tmp = abs(ptB - ptA)
            if tmp < dist:
                dist = tmp
                J = j
        AtoB.append(dist*sqrt(float(i+1)*(J+1)))
    # Distances from all points in peaksA to closest in peaksB.
    # Additional weighting factor sqrt(float(i+1)*(J+1) intriduce to prefer low-resolution similarity.
    BtoA = []
    dist = 0
    for i, ptB in enumerate(peaksB):
        J = 0
        dist = abs(ptB - peaksA[0])
        for j, ptA in enumerate(peaksA):
            tmp = abs(ptB - ptA)
            if tmp < dist:
                dist = tmp
                J = j
        BtoA.append(dist*sqrt(float(i+1)*(J+1)))
    # Get Gaussian weights.
    # Form of the weight function is
    # Result = 0.5/len(AtoB)Sum_i exp(-(AtoB_i)^2/2sgm^2) +
    #         0.5/len(BtoA)Sum_i exp(-(BtoA_i)^2/2sgm^2)
    
    Result = 0
    sgm2 = sgm**2
    for elem in AtoB:
        Result += exp(-0.5*elem**2/sgm2)
    Result /= len(AtoB)
    Result = Result**5
    tmpRes = 0
    for elem in BtoA:
        tmpRes += exp(-0.5*elem**2/sgm2)
    tmpRes /= len(BtoA)
    tmpRes = tmpRes**5
    Result += tmpRes
    Result *= 0.5*Norm

    return Result


def GetDW(PeaksData, sgm, metric="sum"):
    # Function returns list of k nearest to i neighbours.
    # The neighbourhood is mutual, so i is neighrour for each
    # of those k points.
    N = len(PeaksData)
    W = np.zeros((N,N))
    for i, elemA in enumerate(PeaksData):
        for j, elemB in enumerate(PeaksData):
            if j > i:
                if metric == "sum":
                    tmp = SimilarSum(elemA, elemB, sgm)
                elif metric == "prod":
                    tmp = SimilarProd(elemA, elemB, sgm)
                    
                W[i][j] = tmp
                W[j][i] = tmp
            elif i == j:
                W[i][i] = 1

    return W




def mknn(W, k):
    '''
    W - similarity matrix.
    k - number of mutual nearest neighbours.
    '''
    Rslt = []
    '''
    k nearest neighbours matrix.
    # Has dimensions of k by len(W), and coordinates [i,j] 
    # of nearest neighbours to the p-th element in row "p".
    '''    
    # Temporary array of k closest neightbours in a row.
    # Smart sorting algo. O(k*n^2) steps instead of O(n^3)
    for i, row in enumerate(W):
        W[i][i] = 0
        tmp_k = [[elem, j] for j, elem in enumerate(row)] 
      
        tmp_k.sort(key = lambda x: x[0], reverse = True)
        Rslt.append(tmp_k[:k]) 
    '''
    At this point Rslt contains k nearest neighbours
    for each data point. Find the mutual ones and 
    modify W accordingly.
    '''
    # knn
    for i, row in enumerate(W):
        for j, elem in enumerate(row):
            if CmprScnd(j, Rslt[i]) and CmprScnd(i, Rslt[j]):
                continue
            else:
                W[i][j] = 0

def CmprScnd(j, lst):
    found = False
    for it in lst:
        if it[1] == j:
            found = True
            break

    return found 
                

   
         
