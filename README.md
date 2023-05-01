# Quantum Algorithms for Shapley Value Calculation

## Authors

<a href="https://github.com/iain-burge/iain-burge">Iain Burge</a>

<a href="https://carleton.ca/scs/people/michel-barbeau/">Michel Barbeau</a>

<a href="http://www-public.imtbs-tsp.eu/~garcia_a/web/">Joaquin Garcia-Alfaro</a>

## Resources

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/python/">Python Code</a>.

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">Matlab Code</a>.

## Summary of Code and Results 

### One equation model

#### Construct the two quantum systems (may be combined in one)

[![](data:image/png;base64,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){width="129" height="19"}]{texencoding="B^{\pm} (H^{\otimes n}\otimes I) \vert 0 \rangle^{\otimes n+1}
" style="vertical-align:-5px"}

** TODO **

### Random Voting Games & Quantum Shapley Values

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

## References

If using this code for research purposes, please cite:

Iain Burge, Michel Barbeau and Joaquin Garcia-Alfaro. Quantum Algorithms for Shapley Value Calculation. *To appear*. May 2023.

```
@inproceedings{burge-barbeau-alfaro2023Shapley,
  title={Quantum Algorithms for Shapley Value Calculation},
  author={Burge, Iain and Barbeau, Michel and Garcia-Alfaro, Joaquin},
  booktitle={To appear},
  pages={1--8},
  year={2023},
  month={May},
}
```


 
