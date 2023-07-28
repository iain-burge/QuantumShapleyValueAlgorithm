
import numpy as np
import matplotlib.pyplot as plt

from qiskit import QuantumCircuit, transpile, Aer
from qiskit.circuit.library import MCXGate
from shapExampleGenerator import SHAPGenerator
from quantumShapEstimation import QuantumShapleyWrapper as qsw

from tqdm import tqdm

def constructPlusOneGate(numBits: int, numControls: int = 0, name: str = None):
    """
        Gate Structure:
            [Controls] + [Input]

        Returned gate adds one mod 2**numBits to input qubits
    """
    #init circuit
    if name is None:
        name = "+1"
        if numControls > 0:
            name = f"c({numControls}){name}"
    circuit = QuantumCircuit(numBits, name=name)

    #build circuit
    for i in range(numBits-1):
        mcx = MCXGate(numBits-i-1)
        circuit.append(
            mcx, [j for j in range(numBits-1, i-1, -1)]
        )
    circuit.x(numBits-1)

    #convert circuit to gate and return
    if numControls > 0:
        return circuit.to_gate().control(numControls)
    else:
        return circuit.to_gate()
 
def constructFixedAdditionGate(numBits: int, increment: int, 
        numControls: int = 0, name: str = None
    ):
    if name is None:
        name = f"+{increment}"
        if numControls > 0:
            name = f"c({numControls}){name}"

    circuit = QuantumCircuit(numBits, name=name)

    for bit in range(numBits):
        if increment & (1<<bit) != 0:
            gate = constructPlusOneGate(numBits-bit)
            circuit.append(gate, [i for i in range(0, numBits-bit)])

    #print(circuit.draw())
    #convert circuit to gate and return
    if numControls > 0:
        return circuit.to_gate().control(numControls)
    else:
        return circuit.to_gate()

def randomVotingGame(numPlayers: int, thresholdBits: int, roughVariance: int=1):
    """Distributes votes among players such that the grand coalition
      has between threshold and 2*threshold votes.

    Args:
        numPlayers (int): The number of voters in the game
        thresholdBits (int):  The vote threshold to pass a vote, will be forced
            to be 2**(thresholdBits-1). This is so I can be lazy when designing
            the circuit and just check one bit.
        roughVariance (int): Larger means the players will have more
            voting variance
    """
    #Define players voting power arr 
    playerVotingPower  = np.zeros(numPlayers, dtype=int)
    #Define a random order to assign player points
    randomOrderPlayers = np.arange(numPlayers)
    np.random.shuffle(randomOrderPlayers)

    #Due to my being lazy, force threshold to be a nearby power of 2
    threshold  = 2**(thresholdBits-1)

    #Assign each player votes
    totalVotes = int(np.floor((1+np.random.rand()) * threshold))
    while True:
        for i in randomOrderPlayers:
            while np.random.rand() > 0.5:
                if totalVotes <= 0:
                    break
                variance = int(np.ceil(roughVariance * np.random.rand()))
                change = min(totalVotes, variance)
                playerVotingPower[i] += change
                totalVotes -= change

        if totalVotes <= 0:
            break
        np.random.shuffle(randomOrderPlayers)

    return list(playerVotingPower)

def randomVotingGameGate(thresholdBits: int, playerVal: list[float]):
    playerReg = np.arange(len(playerVal)).tolist()
    voteReg   = np.arange(
        len(playerVal), len(playerVal)+thresholdBits).tolist()
    allReg = playerReg + voteReg
    utilityReg = [len(playerVal)]
    
    circuit = QuantumCircuit(len(playerReg) + len(voteReg))

    for player in playerReg:
        circuit.append(
            constructFixedAdditionGate(len(voteReg), playerVal[player], 1),
            [player] + voteReg
        )

    #temp:
    print(circuit.draw())
    #:temp

    return circuit.to_gate(), playerReg, utilityReg, allReg

def classicalVotingShap(
        threshold: int, playerVals: list[int]
    )-> list[float]:
    #Sorry for the nested function, better than writing this as a lambda funciton

    #Start of nested funciton
    def intToCoalitionValue(coalition: int)->float:
        #A function which takes integer encodings of coalitions
        #and outputs their values in a voting game
        votes = 0
        for j in range(len(playerVals)):
            jthPlayerOn = (coalition & (1<<j)) > 0
            if jthPlayerOn:
                #Add player value, Note: endian-ness gets reversed 
                votes += playerVals[len(playerVals)-j-1]
        
        return 1 if votes >= threshold else 0
    #End of nested function

    #Calculate coalition values
    coalitionValues = SHAPGenerator.lambdaGenerateContributions(
        numFactors=len(playerVals),
        generator=intToCoalitionValue
    )
    
    #Calculate shapley values
    sg = SHAPGenerator(
        numFactors=len(playerVals),
        rangeMin=0,
        rangeMax=1,
        contributions=coalitionValues
    )
    shapleyValues = []
    for i in range(len(playerVals)):
        shapleyValues.append(sg.computeShap(i))

    return shapleyValues

def quantumVotingShap(
        threshold: int, playerVals: list[int], ell: int = 2,
    )-> list[float]:
    """Uses quantum circuit to estimate shapley values of voting game

    Args:
        threshold (int): Number of votes needed for a successful vote.
          Make sure threshold is a power of 2 - since I was lazy with the circuit.
        playerVals (list[int]): The number of votes attributed to each player.

    Returns:
        list[float]: List of Shapley values of players
    """
    thresholdBits = int(np.floor(np.log2(threshold))+1)

    gate, playerReg, utilityReg, allReg = randomVotingGameGate(
        thresholdBits=thresholdBits,
        playerVal=playerVals,
    )

    votingQShapWrapper = qsw(
        gate=gate,
        factorIndices=playerReg,
        outputIndices=utilityReg,
        fullIndices=allReg
    )

    numPlayers = len(playerVals)
    qshaps = numPlayers*[0]
    for i in range(numPlayers):
        qshaps[i] = votingQShapWrapper.approxShap(
            i, rangeMin=0, rangeMax=1, betaApproxBits=ell, directProb=True
        )
    
    return qshaps


def main():
    thresholdBits = 3
    numPlayers    = 3

    threshold = 2**(thresholdBits-1)

    rvgArray = randomVotingGame(
        numPlayers=numPlayers,
        thresholdBits=thresholdBits,
        roughVariance=2,
    )

    #quantum Shapley
    qshaps = quantumVotingShap(
        threshold=threshold,
        playerVals=rvgArray,
    )

    #classical Shapley
    cshaps = classicalVotingShap(
        threshold=threshold,
        playerVals=rvgArray,
    )

    print("Player Values:  ", rvgArray)
    print("Threshold:      ", 2**(thresholdBits-1))
    print(80*"=")
    print("Shapley Values: ")
    for i in range(numPlayers):
        print(f"\tPlayer {i}:")
        print(f"\t\tqshaps[{i}] = {qshaps[i]:.4f}")
        print(f"\t\tcshaps[{i}] = {cshaps[i]:.4f}")
        
    # plt.hist(rvgArray)
    # plt.show()

if __name__ == "__main__":
    main()

