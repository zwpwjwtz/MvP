#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 09:37:25 2022

Python script for find difference in m/z of two spectra
"""


# Utility functions
CONST_INFINITY = float('Inf')

def compareScalarVectorLessOrEqual(scalar: float, vector: list) -> list:
    return list(map(lambda X, Y: X <= Y, [scalar] * len(vector), vector))
def compareScalarVectorGreaterOrEqual(scalar: float, vector: list) -> list:
    return list(map(lambda X, Y: X >= Y, [scalar] * len(vector), vector))

def vectorAND(vector1: list, vector2: list) -> list:
    return list(map(lambda X, Y: X & Y, vector1, vector2))

def matchOrderedList(probe: float, data: list, tolerance = 0) -> list:
    # Required lower & upper bound
    lowerBound = probe - tolerance
    upperBound = probe + tolerance
    if data[0] > upperBound or data[len(data) - 1] < lowerBound:
        return []
    
    # Bisect search for the actual lower bound
    leftLowerPos = 0
    leftUpperPos = len(data) - 1
    while leftLowerPos < leftUpperPos:
        midPos = int((leftLowerPos + leftUpperPos) * 0.5)
        if data[midPos] > lowerBound:
            leftUpperPos = midPos - 1
        elif data[midPos] < lowerBound:
            leftLowerPos = midPos + 1
        else:
            leftUpperPos = leftLowerPos = midPos
    
    # Bisect search for the actual upper bound
    rightLowerPos = max(leftLowerPos, leftUpperPos)
    rightUpperPos = len(data) - 1
    while rightLowerPos < rightUpperPos:
        midPos = int((rightLowerPos + rightUpperPos) * 0.5)
        if data[midPos] > upperBound:
            rightUpperPos = midPos - 1
        elif data[midPos] < upperBound:
            rightLowerPos = midPos + 1
        else:
            rightUpperPos = rightLowerPos = midPos
    
    # Find the lower and the upper bound
    return list(range(max(leftLowerPos, leftUpperPos), 
                      min(rightLowerPos, rightUpperPos) + 1))

def matchWithAbsoluteTolerance(source: float, 
                               target: list, 
                               tolerance: float,
                               targetIsOrdered = False) -> int:
    if targetIsOrdered:
        matched = matchOrderedList(source, target, tolerance)
        if len(matched) == 0:
            return -1
        delta = list(abs(target[i] - source) for i in matched)
        return matched[delta.index(min(delta))]
    
    # For unordered list: using element-wise comparison
    matched = vectorAND(
                 compareScalarVectorLessOrEqual(source - tolerance, 
                                                target),
                 compareScalarVectorGreaterOrEqual(source + tolerance, 
                                                   target))
    delta = list(map(lambda X, index: abs(target[index] - source) if X 
                                      else CONST_INFINITY, 
                     matched, range(0, len(target))))
    minDelta = min(delta)
    if (minDelta < CONST_INFINITY):
        return delta.index(minDelta)
    else:
        return -1

def matchWithRelativeTolerance(source: float, 
                               target: list, 
                               tolerance: float,
                               targetIsOrdered = False) -> int:
    return matchWithAbsoluteTolerance(source, target, 
                                      source * tolerance, targetIsOrdered)

# The target list must be of ascending order
def mapUniqueMZ(source: list, 
                target: list, 
                tolerance = 1e-3, 
                relativeTolerance = True,
                noMatchAsNone = False) -> list:
    
    # Find all matched target m/z for each source m/z
    matched = list(map(matchWithRelativeTolerance if relativeTolerance
                       else matchWithAbsoluteTolerance, 
                       source, 
                       [target] * len(source), 
                       [tolerance] * len(source),
                       [True] * len(source)))
    
    if noMatchAsNone:
        matched = list(map(lambda X: X if X >= 0 else None, matched))
        
    return list(filter(lambda X: X != -1, matched))

# mzList must be a list of (ascendingly) ordered list of m/z
def clusterMZ(mzList: list, 
              tolerance = 1e-3, 
              relativeTolerance = True) -> list:
    if (len(mzList) < 1):
        return list()
    if (len(mzList) == 1):
        return list(range(0, len(mzList[0])))
    
    mzCluster = mzList[0].copy()
    mzIndexList = list()
    mzIndexList.append(list(range(0, len(mzList[0]))))
    for i in range(1, len(mzList)):
        mzIndex = mapUniqueMZ(mzList[i], mzCluster,
                              tolerance = tolerance, 
                              relativeTolerance = relativeTolerance,
                              noMatchAsNone = True)
        unmappedIndex = list(filter(lambda j: mzIndex[j] == None, 
                                    range(0, len(mzIndex))))
        if len(unmappedIndex) > 0:
            mzCluster.extend(mzList[i][j] for j in unmappedIndex)
            for i, X in \
                zip(unmappedIndex, range(len(mzCluster) - len(unmappedIndex), 
                                         len(mzCluster))):
                mzIndex[i] = X
        mzIndexList.append(mzIndex)
    
    return mzIndexList

def averageMZCluster(index, data, indexList, threshold, algorithm):
     cluster = list()
     for i in range(0, len(indexList)):
         indexMask = list(map(lambda X: X == index, indexList[i]))
         for j in range(0, len(data[i])):
             if indexMask[j]:
                 cluster.append(data[i][j])
     if len(cluster)/len(data) >= threshold:
         return algorithm(cluster)
     else:
         return None
 
def averageMZ(mzList: list, 
              tolerance = 1e-3, 
              relativeTolerance = True,
              minFrequency = 1, 
              algorithm = lambda X: sum(X) / len(X)) -> list:
    if len(mzList) < 1:
        return list()
    
    # Cluster m/z and get index maps for all m/z vectors
    clusterIndexList = clusterMZ(mzList, tolerance, relativeTolerance)
    
    # Get all available cluster indexes
    allIndexes = set(X for Y in clusterIndexList for X in Y)
    
    # Apply the statistical algorithm on each cluster
    averagedMZ = list(map(averageMZCluster,
                          allIndexes,
                          [mzList] * len(allIndexes),
                          [clusterIndexList] * len(allIndexes),
                          [minFrequency] * len(allIndexes),
                          [algorithm] * len(allIndexes)))
    return list(X for X in filter(lambda Y: Y != None, averagedMZ))
