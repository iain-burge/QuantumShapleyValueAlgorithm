
from math import comb
import random as rand
import numpy as np

from functools import cache

from typing import List, Dict, Union, Callable

class SHAPGenerator():
    def __init__(self, numFactors: int, rangeMin: float, rangeMax: float, 
            contributions: Dict[str, int] = None):
        self.__numFactors    = numFactors
        self.__numSubsets    = 2 ** self.__numFactors
        self.__rangeMax      = rangeMax 
        self.__rangeMin      = rangeMin 

        if contributions is None:
            self.__contributions = \
                self.__generateContributions(rangeMin, rangeMax)
        else:
            self.__contributions = contributions

    def __generateContributions(self, rangeMin: int, rangeMax: int)\
            -> Dict[str, float]:
        #Declaring Variables
        contributions = {self.__numFactors*"0": 0}

        #For each subset
        for i in range(1,self.__numSubsets):
            bits = self.__subsetIndexToBin(i)

            #Map binary string to random contribution
            contributions[bits] = rand.randrange(rangeMin, rangeMax)

        return contributions

    def __subsetIndexToBin(self, subset: int, numbits: int = None):
        if numbits is None:
            numbits = self.__numFactors

        #Convert subset to binary string
        bits = ""
        for j in range(numbits):
            bits = str((subset&(1<<j))>>j) + bits
        return bits

    @staticmethod
    def __subsetSize(subset: str):
        size = 0
        for i in range(len(subset)):
            size += int(subset[i])
        return size

    @staticmethod
    def lambdaGenerateContributions(
            numFactors: int, generator: Callable[[int], float]
        ) -> dict[str, float]:
        return {f'{i:0{numFactors}b}': generator(i) for i in range(2**numFactors)}

    def computeShap(self, targetFactor: int) -> float:
        shap = 0

        subsetsExcludingT = self.__numSubsets>>1
        #For all subsets not including the target factor
        for i in range(subsetsExcludingT):
            #Convert to a binary string
            strI = self.__subsetIndexToBin(i, self.__numFactors-1)

            #Subset
            subset   = strI[:targetFactor] + "0" + strI[targetFactor:]
            #Subset union target
            subsetUt = strI[:targetFactor] + "1" + strI[targetFactor:]

            #Determine subset size
            sSize = self.__subsetSize(subset)

            #Calculate coefficient
            shapCoef =\
                1/(comb(self.__numFactors - 1, sSize)*(self.__numFactors))

            #Add to shapvalue
            shap += shapCoef *\
                (self.__contributions[subsetUt] - self.__contributions[subset])

        return shap

    @staticmethod
    @cache
    def computeLogShapApproxCoef(n:int,m:int,l:int):
        coef = 0
        binomDist = lambda p: ((p)**(m)) * ((1-p)**(n-m))
        sampleVal = lambda k: (np.cos((np.pi*(2*k+1))/(4*l)))**2
        sampleWei = lambda k: np.pi * np.sin(np.pi*(2*k+1)/(2*l)) / (2*l)
        #sampleWei = lambda k: abs(sampleVal(k+1)-sampleVal(k-1))/2

        for k in range(l):
            coef += sampleWei(k) * binomDist(sampleVal(k))
        return coef

    @staticmethod
    @cache
    def computeLogShapApproxCoefList(n:int, l:int) -> List[float]:
        coefs = (n+1)*[0]
        for m in range(n+1):
            coefs[m] = SHAPGenerator.computeLogShapApproxCoef(n, m, l)

        return coefs

    def computeLogBinomialApproxShap(self, targetFactor: int, l: int) -> float:
        shap = 0
        shapCoefList = self.computeLogShapApproxCoefList(
            self.getNumFactors()-1, l)

        subsetsExcludingT = self.__numSubsets>>1
        #For all subsets not including the target factor
        for i in range(subsetsExcludingT):
            #Convert to a binary string
            strI = self.__subsetIndexToBin(i, self.__numFactors-1)

            #Subset
            subset   = strI[:targetFactor] + "0" + strI[targetFactor:]
            #Subset union target
            subsetUt = strI[:targetFactor] + "1" + strI[targetFactor:]

            #Determine subset size
            sSize = self.__subsetSize(subset)

            #Caclulate coefficient
            shapCoef = shapCoefList[sSize]

            #Add to shapvalue
            shap += shapCoef *\
                (self.__contributions[subsetUt] - self.__contributions[subset])

        return shap
 
    @staticmethod
    @cache
    def computeShapApproxCoef(n:int,m:int,l:int):
        coef = 0
        for k in range(l):
            coef += (((2*k+1)/(2*l))**(m)) * ((1-(2*k+1)/(2*l))**(n-m))
        return coef / l

    @staticmethod
    @cache
    def computeShapApproxCoefList(n:int, l:int) -> List[float]:
        coefs = (n+1)*[0]
        for m in range(n+1):
            coefs[m] = SHAPGenerator.computeShapApproxCoef(n, m, l)
        return coefs

    def computeBinomialApproxShap(self, targetFactor: int, l: int) -> float:
        shap = 0
        shapCoefList = self.computeShapApproxCoefList(self.getNumFactors()-1, l)

        subsetsExcludingT = self.__numSubsets>>1
        #For all subsets not including the target factor
        for i in range(subsetsExcludingT):
            #Convert to a binary string
            strI = self.__subsetIndexToBin(i, self.__numFactors-1)

            #Subset
            subset   = strI[:targetFactor] + "0" + strI[targetFactor:]
            #Subset union target
            subsetUt = strI[:targetFactor] + "1" + strI[targetFactor:]

            #Determine subset size
            sSize = self.__subsetSize(subset)

            #Caclulate coefficient
            shapCoef = shapCoefList[sSize]

            #Add to shapvalue
            shap += shapCoef *\
                (self.__contributions[subsetUt] - self.__contributions[subset])

        return shap
    
    def shapToProbability(self, shap: float) -> float:
        shiftedShap = shap - self.__rangeMin
        return shiftedShap / self.getRangeLen()

    def probabilityToShap(self, prob) -> float:
        shiftedShap = prob * self.getRangeLen()
        return shiftedShap + self.__rangeMin

    def getTotalValue(self) -> float:
        return self.getContribution(self.__numSubsets - 1)

    def getContribution(self, key: Union[int, str]) -> float:
        if type(key) is int:
            key = f"{key:0{self.__numFactors}b}"
        return self.__contributions[key]

    def getNumFactors(self) -> int:
        return self.__numFactors
    def getNumSubsets(self) -> int:
        return self.__numSubsets
    def getContributions(self) -> Dict[str, float]:
        return self.__contributions.copy()
    def getRangeLen(self) -> float:
        return self.__rangeMax - self.__rangeMin

def main():
    rand.seed(0)
    numFactors = 10

    sg = SHAPGenerator(numFactors, -128, 127)
    #print(sg.getContributions())

    directShap      = []
    binomialShap    = []
    binomialLogShap = []
    for i in range(numFactors):
        directShap.append(sg.computeShap(i))
        binomialShap.append(sg.computeBinomialApproxShap(i, numFactors-1))
        binomialLogShap.append(
            sg.computeLogBinomialApproxShap(i, numFactors-1))
        print(f"Direct Shap({i}):      \t{directShap[-1]:.2f}")
        print(f"Binomial Shap({i}):    \t{binomialShap[-1]:.2f}\t|\tdelta: "
                + f"{abs(binomialShap[-1]-directShap[-1]):.4f}")
        print(f"Log Binomial Shap({i}):\t{binomialLogShap[-1]:.2f}\t|\tdelta: "
                + f"{abs(binomialLogShap[-1]-directShap[-1]):.4f}")
        print()

    print(f"Sum Direct Shap:      \t{sum(directShap):.2f}")
    print(f"Sum Binomial Shap:    \t{sum(binomialShap):.2f}")
    print(f"Sum Log Binomial Shap:\t{sum(binomialLogShap):.2f}")
    
    #Different test:
    #print("\nTest lambda contribution generator function:")
    #generator = lambda i: i
    #print(f"\t{sg.lambdaGenerateContributions(numFactors, generator) = }")

if __name__ == "__main__":
    main()
