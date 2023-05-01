
import numpy as np
from shapExampleGenerator import SHAPGenerator
from qiskit import QuantumCircuit

def toYRotation(val: float=0):
    return 2*np.arccos(np.sqrt(val))

def binaryStringToXGates(circuit: QuantumCircuit, register: list[int],
        bits: str
    ) -> None:
    for i, bit in enumerate(bits):
        if bit == '0':
            circuit.x(register[i])

def getRandomShapGate(numFactors: int, contributions = dict[str, float], 
        negative: bool=False
    ) -> ...:
    circuit    = QuantumCircuit(numFactors + 1)
    factorsReg = [i for i in range(numFactors)]
    outputReg  = [numFactors]

    for subset, contribution in contributions.items():
        binaryStringToXGates(circuit, factorsReg, subset)
        circuit.mcry(
            (-1 if negative else 1) * toYRotation(1-contribution),
            factorsReg,
            outputReg
        )
        binaryStringToXGates(circuit, factorsReg, subset)

    return circuit.to_gate()

def main():
    #Generate Shap Example
    numFactors = 4
    rangeMin   =-128 
    rangeMax   = 127
    #Contributions
    sg = SHAPGenerator(numFactors, rangeMin, rangeMax)
    contributions = sg.getContributions()
    normContributions = \
        {subset: (cont-rangeMin)/(rangeMax-rangeMin)
        for subset, cont in contributions.items()}

    getRandomShapGate(numFactors, normContributions)

if __name__ == "__main__":
    main()

