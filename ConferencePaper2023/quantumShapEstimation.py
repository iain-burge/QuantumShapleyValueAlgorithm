
import numpy as np
import random
from qiskit import QuantumCircuit, transpile, Aer
from enum import Enum, auto
from shapRandomGates import getRandomShapGate
from shapExampleGenerator import SHAPGenerator

from tqdm import tqdm

class Registers(Enum):
    aux       = auto()
    factors   = auto()
    otherFac  = auto()
    targetFac = auto()
    output    = auto()
    full      = auto()
    classical = auto()

class QuantumShapleyWrapper():
    def __init__(self, gate, factorIndices: list[int], outputIndices: list[int],
        fullIndices: list[float]
    ) -> None:
        """This class acts as a wrapper which holds a gate representing the circuit
        of interest. This class is used to prepare shapley analysis on said gate.
  
        Args:
            gate (Gate): The gate being studied
            factorsIndices (list[int]): A list containing the index of each factor
            outputIndices (list[int]): A list of the output bits in order of significance
            fullIndices (list[int]): A list of all indices in the gate
        """

        self.__gate            = gate
        self.__outputIndices   = outputIndices
        self.__factorIndices   = factorIndices
        self.__fullIndices     = fullIndices
  
    def getShapCircuit(self, targetFactor: int, betaApproxBits: int = None,
        targetOn: int = True
    ) -> QuantumCircuit:

        #Init number of beta approximation bits
        if betaApproxBits is None:
            betaApproxBits = int(np.ceil(np.log2(len(self.__factorIndices))))
  
        #Get registers
        registers = self.getRegisters(targetFactor, betaApproxBits)
  
        #Construct the circuit
        circuit = self.initShapley(targetFactor, betaApproxBits, targetOn)
        circuit.append(self.__gate, registers[Registers.full])
  
        return circuit
  
    def getRegisters(self, targetFactor: int, betaApproxBits: int):
        #Create dictionary with all registers
        regDict = {}
  
        regDict[Registers.aux]       = [i for i in range(betaApproxBits)]
        regDict[Registers.full]      = [i+betaApproxBits for i 
                in self.__fullIndices]
        regDict[Registers.factors]   = [i+betaApproxBits for 
                i in self.__factorIndices]
        regDict[Registers.otherFac]  = [i for i in regDict[Registers.factors] 
                if i!=regDict[Registers.factors][targetFactor]]
        regDict[Registers.targetFac] = [regDict[Registers.factors][targetFactor]]
        regDict[Registers.output]    = [i+betaApproxBits for 
                i in self.__outputIndices]
        regDict[Registers.full]      = [i+betaApproxBits for 
                i in self.__fullIndices]
        regDict[Registers.classical] = [0]
        
        return regDict
  
    def initBetaApproxBits(self, auxReg: list[int], circuit: QuantumCircuit):

        L = 2**len(auxReg)
  
        auxWeights = np.arange(L)
        auxWeights = np.sin(np.pi*(2*auxWeights+1)/(2*L))
        auxWeights = auxWeights/np.sum(auxWeights)
        auxWeights = np.sqrt(auxWeights)
  
        circuit.initialize(auxWeights, auxReg)
  
    def getShapleyInitGate(self, betaApproxBits: int, numFactors: int, 
        target: int, targetOn: bool = True
    ) -> ...:

        #Define Subregisters
        auxReg       = [i for i in range(betaApproxBits)]
        factorsReg   = [i+betaApproxBits for i in range(numFactors)]
        otherFacReg  = [i for i in factorsReg if i!=factorsReg[target]]
        targetFacReg = [factorsReg[target]]
  
        circuit = QuantumCircuit(len(auxReg)+len(factorsReg), 
            name="Shap Factor Init")
  
        #Initializing Target Factor Qubit
        if targetOn:
            circuit.x(targetFacReg)
  
        #Initializing Other Factor Qubits
        for factorQubit in otherFacReg:
            circuit.ry(np.pi/(2**(betaApproxBits+1)), factorQubit)
            for k, auxQubit in enumerate(auxReg):
                circuit.cry(np.pi*2**-(betaApproxBits-k), auxQubit, factorQubit)

        #print(circuit.draw())
        
        return circuit.to_gate()

    def initShapley(self, targetFactor: int, betaApproxBits: int,
            targetOn: bool = True
        ) -> None:
        """This function returns a gate which prepares a circuit to have probability
        amplitudes corresponding to Shapley coefficients.

        Args:
            numFactors (int): The number of total factors of a system
            targetFactor (int): The factor of interest
            auxBits (int): The resolution of our beta function approximation, will 
                be forced nearest power of 2 greater than L
        """
        #Define Registers
        regDict = self.getRegisters(targetFactor, betaApproxBits)

        factorBits = len(regDict[Registers.factors])
        L          = 2**betaApproxBits

        #Create Circuit
        circuit = QuantumCircuit(betaApproxBits + self.__gate.num_qubits, 1,
            name="SHAP State Init")

        #Initializing Aux Qubits
        self.initBetaApproxBits(regDict[Registers.aux], circuit)

        #Initializing Factors
        factorInitGate = self.getShapleyInitGate(
            betaApproxBits, factorBits, targetFactor, targetOn)
        circuit.append(factorInitGate, 
            regDict[Registers.aux]+regDict[Registers.factors])
    
        #print(circuit.draw())
        return circuit

    def approxShap(self, targetFactor: int, numMeasurements: int = 2**12, 
        betaApproxBits: int = None, rangeMin: float = -128, 
        rangeMax: float = 127, directProb: bool = False
    ) -> float:
        ON  = 1
        OFF = 0
        
        if betaApproxBits is None:
            betaApproxBits = int(np.ceil(np.log2(len(self.__factorIndices))))
        reg = self.getRegisters(targetFactor, betaApproxBits)

        probs   = 2*[0]
        counts  = 2*[0]
        for toggle in [ON, OFF]:
            #Initialize Circuit
            circuit = self.getShapCircuit(
                targetFactor, betaApproxBits, toggle==ON)
            #Measure Circuit Output
            #if not directProb:
            #    circuit.measure(reg[Registers.output], reg[Registers.classical])

            #Use Aer
            sim = Aer.get_backend('aer_simulator')
            circuit.save_statevector()

            #Compile the circuit down to low-level QASM instructions
            compiled_circuit = transpile(circuit, sim)
            
            #Simulate circuit
            result = sim.run(compiled_circuit, 
                shots=numMeasurements #if not directProb else 1
            ).result()
            out_state = result.get_statevector(circuit, decimals=4)

            #Retrieve probs
            probs[toggle] = out_state.probabilities(reg[Registers.output])
            #Compensating for rounding errors:
            probs[toggle][1] = max(0, min(1, probs[toggle][1]))

            counts[toggle] = np.random.binomial(numMeasurements, probs[toggle][1])

        if directProb:
            shap = (rangeMax-rangeMin)*(probs[ON][1]-probs[OFF][1])
        else:
            shap = (rangeMax-rangeMin)*(counts[ON]-counts[OFF])/numMeasurements

        return shap

def main():
    ##temp
    #bits = 5
    #betaApproxBits = 2
    #targetFactor = 0
    #circuit = QuantumCircuit(bits, name="test")
    #for i in range(1,bits):
    #    circuit.cry(np.pi*i/(2**(bits-1)), i, 0)
    #gate = circuit.to_gate()
    ##temp

    #qsw = QuantumShapleyWrapper(
    #    gate, [i for i in range(1,bits)], [0], [i for i in range(bits)])
    #print(qsw.getRegisters(targetFactor, betaApproxBits))
    ## qsw.initShapley(targetFactor, betaApproxBits)
    #shapCircuit = qsw.getShapCircuit(targetFactor, betaApproxBits)
    #print(shapCircuit.draw())

    #random.seed(0)
    #Generate Shap Example
    numMeasurements = 2**16
    directProb = True
    betaApproxBits  = 4
    numFactors = 6
    rangeMin   =-128
    rangeMax   = 127
    #Contributions
    sg = SHAPGenerator(numFactors, rangeMin, rangeMax)
    contributions = sg.getContributions()
    normContributions = \
        {subset: (cont-rangeMin)/(rangeMax-rangeMin)
        for subset, cont in contributions.items()}
    gate = getRandomShapGate(numFactors, normContributions)

    #Put in wrapper class
    qsw = QuantumShapleyWrapper(
        gate, [i for i in range(numFactors)], 
        [numFactors], [i for i in range(numFactors+1)]
    )

    #Generate shap approx circuit
    classicalShaps = numFactors * [0]
    quantumShaps   = numFactors * [0]
    biApproxShaps  = numFactors * [0]
    deltas         = numFactors * [0]

    #Collect Results
    for i in tqdm(range(numFactors)):
        classicalShaps[i] = sg.computeShap(i)
        biApproxShaps[i]  = sg.computeLogBinomialApproxShap(i, 2**betaApproxBits)
        quantumShaps[i]   = qsw.approxShap(
            i, betaApproxBits=betaApproxBits, numMeasurements=numMeasurements,
            directProb=directProb, rangeMin=rangeMin, rangeMax=rangeMax)
        deltas[i]        = abs(classicalShaps[i] - quantumShaps[i])
    print(30*"=")

    #Print Results
    for i in range(numFactors):
        print(f"Factor {i}:")
        print(f"\tClassical Shap: {classicalShaps[i]:.2f}")
        print(f"\tBiApprox  Shap: {biApproxShaps[i]:.2f}")
        print(f"\tQuantum Shap:   {quantumShaps[i]:.2f}")
        print(f"\n\tDelta:          {deltas[i]:.2f}\n")

    #Print Stats
    print(30*"-")
    print(f"{np.average(deltas) = :.4f}\n{np.std(deltas)     = :.4f}")

if __name__ == "__main__":
    main()

