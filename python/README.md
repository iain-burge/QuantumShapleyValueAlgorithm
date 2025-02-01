# Quantum Algorithms for Shapley Value Calculation

## 1. Authors

<a href="https://github.com/iain-burge/iain-burge">Iain Burge</a>

<a href="https://carleton.ca/scs/people/michel-barbeau/">Michel Barbeau</a>

<a href="http://www-public.imtbs-tsp.eu/~garcia_a/web/">Joaquin Garcia-Alfaro</a>

## 2. Available Resources

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/python/">Python Code</a>.

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">Matlab Code</a>.

## 3. Sample Code and Results

### 3.1 Random Voting Games & Quantum Shapley Values

The <a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/python/">python scripts under this folder</a>
address random voting games using our quantum algorithm to estimate the Shapley value of each player. The scripts also perform
some basic data analysis on the predictions. A summary of the code and the results is shown below.

#### Import Libraries

``` python
import quantumBasicVotingGame as vg
from quantumShapEstimation import QuantumShapleyWrapper as qsw

import numpy as np
import matplotlib.pyplot as plt
import pickle
from tqdm.auto import tqdm
```
#### Define Variables

``` python
numTrails = 32

maxEll    = 7

#Defining the different conditions
numPlayersCond    = [4,8,12]
thresholdBitCond  = [4,5,6]
roughVarianceCond = [1,2,2]
```

#### Run Simulations

``` python
simulations = {}

for trialNum in tqdm(range(numTrails), desc="Current Trial"):
    for ell in tqdm(range(1,maxEll), desc="Current Ell"):
        for n, thresholdBits, roughVariance in zip(
            numPlayersCond, thresholdBitCond, roughVarianceCond
        ):
            trial = (n,ell,trialNum)


            #New random game
            threshold = 2**(thresholdBits-1)

            playerVals = vg.randomVotingGame(
                numPlayers=n,
                thresholdBits=thresholdBits,
                roughVariance=roughVariance
            )

            #quantum Shapley
            qshaps = vg.quantumVotingShap(
                threshold=threshold,
                playerVals=playerVals,
                ell=ell
            )

            #classical Shapley
            cshaps = vg.classicalVotingShap(
                threshold=threshold,
                playerVals=playerVals,
            )

            #Store outcome
            simulations[trial] = (qshaps, cshaps)

```

#### Save Results

``` python
with open('shapleyVoteResults.pkl', 'wb') as f:
    pickle.dump(simulations, f)
```
#### Analyze Trials

``` python
def meanAbsError(qshaps, cshaps):
    err = 0
    for qshap, cshap in zip(qshaps, cshaps):
        err += abs(qshap-cshap)
    return err
```
``` python
plt.rcParams['figure.figsize'] = [12, 5]
fig, ax = plt.subplots(1, len(numPlayersCond))

#We're looking to find reciprocal mean abs error per trial
#For each trial with n players
for i, n in enumerate(numPlayersCond):
    #Orient data
    resultsX = []
    resultsY = []
    resultErr = []
    for ell in range(1, maxEll):
        trialOutcomes = []

        for trialNum in range(numTrails):
            qshaps, cshaps = simulations[(n,ell,trialNum)]
            trialOutcomes.append(
                meanAbsError(qshaps, cshaps)
            )

        trialOutcomes = np.array(trialOutcomes)
        resultsX.append(ell)
        resultsY.append(trialOutcomes.mean())
        resultErr.append(trialOutcomes.std())

        # resultsX += len(trialOutcomes) * [ell]
        # resultsY += trialOutcomes

    ax[i].set_title(f"{n} Players")#, Threshold: {2**thresholdBitCond[i]}")
    ax[i].bar(
        np.array(resultsX),
        1/np.array(resultsY),
        # yerr=resultErr,
        align='center',
        alpha=0.5,
        ecolor='black',
        capsize=10,
    )
    ax[i].set_xlabel(r"$\ell$")
    ax[i].set_ylabel(r"Reciprocal Mean Absolute Error")

plt.tight_layout()
plt.show()
```

![](figures/546.png)

## Minimal Requirements

The code has been tested using the following set of libraries and requirements:

- Python 3.9.2
- IPython 8.14.0
- ipykernel 6.25.0
- ipywidgets 8.1.0
- jupyter_client 8.3.0
- jupyter_core 5.3.1
- jupyter_server 2.7.0
- jupyterlab 4.0.4
- nbclient 0.8.0
- nbconvert 7.7.3
- nbformat 5.9.2
- notebook 7.0.1
- qtconsole 5.4.3
- traitlets 5.9.0
- numpy 1.25.2
- matplotlib 3.7.2
- qiskit 0.44.0
- tqdm 4.65.0
- qiskit aer-0.12.2

Additionally, the following tutorial for amplitude estimation was used as a framework: 
https://qiskit.org/ecosystem/finance/tutorials/00_amplitude_estimation.html 

## References

If using this code for research purposes, please cite:

Iain Burge, Michel Barbeau and Joaquin Garcia-Alfaro. Quantum Algorithms for Shapley Value Calculation. 2023 IEEE International Conference on Quantum Computing and Engineering (QCE 2023), Bellevue, WA, United States, September 17-22, 2023.

Iain Burge, Michel Barbeau and Joaquin Garcia-Alfaro. A Shapley Value Estimation Speedup for Efficient Explainable Quantum AI. arXiv 2412.14639.

```
@inproceedings{burge-barbeau-alfaro2023Shapley,
  title={Quantum Algorithms for Shapley Value Calculation},
  author={Burge, Iain and Barbeau, Michel and Garcia-Alfaro, Joaquin},
  booktitle={2023 IEEE International Conference on Quantum Computing and Engineering (QCE 2023), Bellevue, WA, United States, September 17-22, 2023},
  pages={1--9},
  year={2023},
  month={September},
}

@misc{burge2024shapleyvalueestimationspeedup,
  title={A Shapley Value Estimation Speedup for Efficient Explainable Quantum AI},
  author={Iain Burge and Michel Barbeau and Joaquin Garcia-Alfaro},
  year={2024},
  eprint={2412.14639},
  archivePrefix={arXiv},
  primaryClass={cs.CR},
  url={https://arxiv.org/abs/2412.14639},
}
```
