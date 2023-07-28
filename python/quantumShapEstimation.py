
import numpy as np
import matplotlib.pyplot as plt

from enum import Enum, auto
from tqdm import tqdm

from shapRandomGates import getRandomShapGate
from shapExampleGenerator import SHAPGenerator

from qiskit.circuit.library import StatePreparation, QFT, ZGate
from qiskit import QuantumCircuit, transpile, Aer
from qiskit.circuit.gate import Gate


class Registers(Enum):
    aux       = auto()
    factors   = auto()
    otherFac  = auto()
    targetFac = auto()
    output    = auto()
    ampEst    = auto()
    full      = auto()
    classical = auto()

class QuantumShapleyWrapper():
    #Constants
    ON  = 1
    OFF = 0


    def __init__(self, gate: Gate, factorIndices: list[int], 
        outputIndices: list[int], fullIndices: list[float]
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
  
    def getRegisters(
            self, targetFactor: int, betaApproxBits: int, ampEstBits: int = 0
        ) -> dict[Registers, list[int]]:
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
        regDict[Registers.ampEst]    = [i+betaApproxBits+len(self.__fullIndices) 
                for i in range(ampEstBits)]
        regDict[Registers.classical] = [0]
        
        return regDict
  
    @staticmethod
    def initBetaApproxBits(auxReg: list[int], circuit: QuantumCircuit) -> None:
        L = 2**len(auxReg)
  
        auxWeights = np.arange(L)
        auxWeights = np.sin(np.pi*(2*auxWeights+1)/(2*L))
        auxWeights = auxWeights/np.sum(auxWeights)
        auxWeights = np.sqrt(auxWeights)

        auxWeights = auxWeights.astype('complex128')
  
        prep = StatePreparation(
            params=auxWeights,
        )

        circuit.append(prep, auxReg)
  
    @staticmethod
    def getShapleyInitGate(betaApproxBits: int, numFactors: int, 
        target: int, targetOn: bool = True
    ) -> Gate:

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

        # print(circuit.draw())
        
        return circuit.to_gate()

    def initShapley(self, targetFactor: int, betaApproxBits: int,
            targetOn: bool = True
        ) -> QuantumCircuit:
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
        if betaApproxBits is None:
            betaApproxBits = int(np.ceil(np.log2(len(self.__factorIndices))))
        reg = self.getRegisters(targetFactor, betaApproxBits)

        probs   = 2*[0]
        counts  = 2*[0]
        for toggle in [self.ON, self.OFF]:
            #Initialize Circuit
            circuit = self.getShapCircuit(
                targetFactor, betaApproxBits, toggle==self.ON)
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
            shap = (rangeMax-rangeMin)*(probs[self.ON][1]-probs[self.OFF][1])
        else:
            shap = (rangeMax-rangeMin)*(counts[self.ON]-counts[self.OFF])\
                /numMeasurements

        return shap

    def constructAGate(
        self, auxReg: list[int], factorReg: list[int], targetFactor: int,
        targetOn: bool
    ) -> Gate:
        #Defining the A Gate
        ACircuit = QuantumCircuit(len(auxReg+factorReg))

        #Prepare auxiliary register
        self.initBetaApproxBits(auxReg, ACircuit)

        ACircuit.append(
            self.getShapleyInitGate(
                betaApproxBits=len(auxReg),
                numFactors=len(self.__factorIndices),
                target=targetFactor,
                targetOn=targetOn,
            ),
            auxReg + factorReg
        )

        return ACircuit.to_gate(label="A")
        
    def constructWGate(self, fullReg: list[int]) -> Gate:
        #Circuit on the full register + output
        WCircuit = QuantumCircuit(len(fullReg))
        WCircuit.append(self.__gate, range(len(fullReg)))

        return WCircuit.to_gate(label="W")

    def constructVGate(
        self, reg: dict[Registers, list[int]]
    ) -> Gate:
        VCircuit = QuantumCircuit(len(reg[Registers.aux]+reg[Registers.full]))
        VCircuit.z(reg[Registers.output])

        return VCircuit.to_gate(label="V")

    def constructS0Gate(
        self, reg: dict[Registers, list[int]]
    ) -> Gate:
        S0Circuit = QuantumCircuit(
            len(reg[Registers.aux]+reg[Registers.full]))
        mczGate = ZGate().control(
            num_ctrl_qubits=len(reg[Registers.aux]+reg[Registers.factors]),
            label="Z"
        )

        S0Circuit.x(reg[Registers.aux]+reg[Registers.full])
        S0Circuit.append(
            mczGate,
            reg[Registers.aux] + reg[Registers.full]
        )
        S0Circuit.x(reg[Registers.aux]+reg[Registers.full])
        
        return S0Circuit.to_gate(label="S_0")

    def constructQGate(
        self, reg: dict[Registers, list[int]], AGate: Gate, WGate: Gate, 
        VGate: Gate, S0Gate: Gate
    ) -> Gate:
        QCircuit = QuantumCircuit(len(reg[Registers.aux]+reg[Registers.full]))

        #Define a phase negation gate
        phaseCircuit = QuantumCircuit(1)
        phaseCircuit.p(np.pi, 0); phaseCircuit.x(0)
        phaseCircuit.p(np.pi, 0); phaseCircuit.x(0)
        phaseGate = phaseCircuit.to_gate(label="-I")

        #Define some inverse gates
        invAGate = AGate.inverse(); invAGate.label = "A^-1"
        invWGate = WGate.inverse(); invWGate.label = "W^-1"

        #Main body of Q gate
        QCircuit.append(phaseGate, [0])
        QCircuit.append(VGate,  reg[Registers.aux]+reg[Registers.full])
        QCircuit.append(invWGate,  reg[Registers.full])
        QCircuit.append(invAGate,  reg[Registers.aux]+reg[Registers.factors])
        QCircuit.append(S0Gate, reg[Registers.aux]+reg[Registers.full])
        QCircuit.append(AGate,  reg[Registers.aux]+reg[Registers.factors])
        QCircuit.append(WGate,  reg[Registers.full])

        print(QCircuit.draw())

        return QCircuit.to_gate(label="Q")
    
    @staticmethod
    def controlledGatePowers(gate: Gate, amplitudeEstBits: int):
        gatePowers = []

        for i in tqdm(range(amplitudeEstBits)):
            gatePowers.append(
                gate.repeat(2**i)
            )
            gatePowers[i].label = f"{gate.label}^{2**i}"
            gatePowers[i] = gatePowers[i].control()
        
        return gatePowers


    def approxShapAmpEst(self, targetFactor: int, numMeasurements: int = 2**4, 
        amplitudeEstBits: int = 3, betaApproxBits: int = None, 
        rangeMin: float = -128, rangeMax: float = 127
    ) -> float:
        """This function makes use of Brassard et al. 2000 and Ashley Montanaro 2017.
           It uses amplitude estimation to quickly estimate the expected value of
           measuring the output bit, from which we can approximate the shapley value.

        Args:
            targetFactor (int): The index of the target player
            numMeasurements (int, optional): The number of times the process 
                is repeated. Defaults to 2**4.
            amplitudeEstBits (int, optional): Number of bits used to perform 
                amplitude estimation. Defaults to None.
            betaApproxBits (int, optional): Number of bits to approximate the 
                beta function distribution. Defaults to None.
            rangeMin (float, optional): Upper bound for coalition value. 
                Defaults to -128.
            rangeMax (float, optional): Lower bound for coalition value. 
                Defaults to 127.

        Returns:
            float: Approximate Shapley value.
        """
  
        if betaApproxBits is None:
            betaApproxBits = int(np.ceil(np.log2(len(self.__factorIndices))))
        reg = self.getRegisters(targetFactor, betaApproxBits, amplitudeEstBits)       

        #Defining target toggle independent gates
        WGate  = self.constructWGate(reg[Registers.full])
        VGate  = self.constructVGate(reg)
        S0Gate = self.constructS0Gate(reg)

        outputProbs = 2*[0]
        for toggle in [self.OFF, self.ON]:
            #Target toggle dependent gates:
            AGate  = self.constructAGate(
                reg[Registers.aux], 
                reg[Registers.factors], 
                targetFactor, 
                targetOn=toggle
            )
            QGate = self.constructQGate(reg, AGate, WGate, VGate, S0Gate)
            QGatePowers = self.controlledGatePowers(QGate, amplitudeEstBits)

            circuit = QuantumCircuit(
                len(reg[Registers.aux]+reg[Registers.full]+reg[Registers.ampEst])
            )

            #Prepare Initial State
            circuit.append(AGate, reg[Registers.aux]+reg[Registers.factors])
            circuit.append(WGate, reg[Registers.full])

            #Prepare Amplitude Estimation Register
            circuit.h(reg[Registers.ampEst])
            
            #Apply Q Multiplicity
            for i in range(amplitudeEstBits):
                circuit.append(
                    QGatePowers[i], 
                    [reg[Registers.ampEst][i]]+reg[Registers.aux]+reg[Registers.full]
                )

            #Inverse Fourier transform on the amplitude estimation register
            qftGate = QFT(num_qubits=amplitudeEstBits, inverse=True).to_gate()
            qftGate.label = "QFT^-1"
            circuit.append(qftGate, reg[Registers.ampEst])
            # circuit.append(QFT(num_qubits=amplitudeEstBits).to_gate(), reg[Registers.ampEst])

            #Run the circuit
            outputProbs[toggle] = assessCircuit(circuit, reg[Registers.ampEst])
            print(circuit.draw())
        
        #Post processing results
        weights = np.sin(
            (np.pi/2**amplitudeEstBits) * np.arange(2**amplitudeEstBits)
        )**2

        shapleyEstimation = np.dot(weights, np.array(outputProbs[self.ON]))\
            - np.dot(weights, np.array(outputProbs[self.OFF]))
        shapleyEstimation *= rangeMax-rangeMin

        print(shapleyEstimation)



def assessCircuit(circuit: QuantumCircuit, register: list[int]):
    #Use Aer
    sim = Aer.get_backend('aer_simulator')
    circuit.save_statevector()

    #Compile the circuit down to low-level QASM instructions
    compiled_circuit = transpile(circuit, sim)

    #Simulate circuit
    result = sim.run(compiled_circuit).result()
    out_state = result.get_statevector(circuit, decimals=4)

    #Visualize output
    probs = out_state.probabilities(register)
    state = out_state.to_dict()

    angles = {
        key: np.sin(np.angle(value)/2)
        for key, value in state.items()
    }

    # print(state)
    # print(angles)
    colors = [f"#{int(180*angle)%256:02x}0080" for angle in angles.values()]
    # print(colors)

    x = np.linspace(0, 1-1/len(probs), len(probs))
    plt.bar(x, probs, align="edge", width=.9/(len(probs)), color=colors)
    plt.xticks(
        ticks=x+.9/(2*len(probs)), 
        labels=[f"{num:0{len(register)}b}" for num in range(len(probs))],
        rotation=90
    )
    plt.show()

    # print(probs)

    return probs


def main():
    #Generate Shap Example
    numMeasurements = 2**4
    directProb = True
    betaApproxBits  = 3
    amplitudeEstBits = 6
    numFactors = 4
    rangeMin   =-128
    rangeMax   = 127

    #Contributions
    sg = SHAPGenerator(
        numFactors, 
        rangeMin, 
        rangeMax, 
        # contributions={'000':0,'001':0,'010':0,'011':0,'100':0,'101':1,'110':1,'111':1},
    )
    contributions = sg.getContributions()
    normContributions = {
        subset: (cont-rangeMin)/(rangeMax-rangeMin)
        for subset, cont in contributions.items()
    }
    gate = getRandomShapGate(numFactors, normContributions)

    #Put in wrapper class
    qsw = QuantumShapleyWrapper(
        gate, [i for i in range(numFactors)], 
        [numFactors], [i for i in range(numFactors+1)]
    )

    for target in [0]:
        qsw.approxShapAmpEst(
            targetFactor=target,
            numMeasurements=numMeasurements,
            amplitudeEstBits=amplitudeEstBits,
            betaApproxBits=betaApproxBits,
            rangeMin=rangeMin,
            rangeMax=rangeMax,
        )
        print(sg.computeShap(target))
    

# def main():
    # #Generate Shap Example
    # numMeasurements = 2**16
    # directProb = True
    # betaApproxBits  = 3
    # numFactors = 5
    # rangeMin   =-128
    # rangeMax   = 127

    # #Contributions
    # sg = SHAPGenerator(numFactors, rangeMin, rangeMax)
    # contributions = sg.getContributions()
    # normContributions = \
    #     {subset: (cont-rangeMin)/(rangeMax-rangeMin)
    #     for subset, cont in contributions.items()}
    # gate = getRandomShapGate(numFactors, normContributions)

    # #Put in wrapper class
    # qsw = QuantumShapleyWrapper(
    #     gate, [i for i in range(numFactors)], 
    #     [numFactors], [i for i in range(numFactors+1)]
    # )

    # #Generate shap approx circuit
    # classicalShaps = numFactors * [0]
    # quantumShaps   = numFactors * [0]
    # biApproxShaps  = numFactors * [0]
    # deltas         = numFactors * [0]

    # #Collect Results
    # for i in tqdm(range(numFactors)):
    #     classicalShaps[i] = sg.computeShap(i)
    #     biApproxShaps[i]  = sg.computeLogBinomialApproxShap(i, 2**betaApproxBits)
    #     quantumShaps[i]   = qsw.approxShap(
    #         i, betaApproxBits=betaApproxBits, numMeasurements=numMeasurements,
    #         directProb=directProb, rangeMin=rangeMin, rangeMax=rangeMax)
    #     deltas[i]        = abs(classicalShaps[i] - quantumShaps[i])
    # print(30*"=")

    # #Print Results
    # for i in range(numFactors):
    #     print(f"Factor {i}:")
    #     print(f"\tClassical Shap: {classicalShaps[i]:.2f}")
    #     print(f"\tBiApprox  Shap: {biApproxShaps[i]:.2f}")
    #     print(f"\tQuantum Shap:   {quantumShaps[i]:.2f}")
    #     print(f"\n\tDelta:          {deltas[i]:.2f}\n")

    # #Print Stats
    # print(30*"-")
    # print(f"{np.average(deltas) = :.4f}\n{np.std(deltas)     = :.4f}")

if __name__ == "__main__":
    main()
