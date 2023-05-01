# Quantum Algorithms for Shapley Value Calculation

## Authors

<a href="https://github.com/iain-burge/iain-burge">Iain Burge</a>

<a href="https://carleton.ca/scs/people/michel-barbeau/">Michel Barbeau</a>

<a href="http://www-public.imtbs-tsp.eu/~garcia_a/web/">Joaquin Garcia-Alfaro</a>

## Resources

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/python/">Python Code</a>.

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">Matlab Code</a>.

## Summary of Results 

### Random Voting Games & Quantum Shapley Values

The followinf results report random voting games and uses our quantum 
algorithm to estimate the Shapley value of each player. It then performs 
some basic data analysis on the predictions.
:::

::: {.cell .markdown}
### Import Libraries
:::

::: {.cell .code execution_count="46"}
``` python
import quantumBasicVotingGame as vg
from quantumShapEstimation import QuantumShapleyWrapper as qsw

import numpy as np
import matplotlib.pyplot as plt
import pickle
from tqdm.auto import tqdm
```
:::

::: {.cell .markdown}
### Define Variables
:::

::: {.cell .code execution_count="47"}
``` python
numTrails = 32

maxEll    = 7

#Defining the different conditions
numPlayersCond    = [4,8,12]
thresholdBitCond  = [4,5,6]
roughVarianceCond = [1,2,2]
```
:::

::: {.cell .markdown}
### Run Simulations
:::

::: {.cell .code execution_count="48"}
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

::: {.output .display_data}
``` json
{"model_id":"353ad19873484d85b191564c5d014cf5","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"0760b9d057a24e8390f8f8e0caebcd01","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"9befaa47222f48fb840ef40fa8033733","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"19bed865d4394f8eb24a31bc1013b1ce","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"36fcf2f9fa854bb3a8a4fc56be39395e","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"3f5317b05d10419e9a9952b7ec1a1e3e","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"67a529282747423b87fc6b88f3c8ddd4","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"16df4f5d03cd4c9d8c3f8af914640195","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"9949a71471964a91976aedf211b0d441","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"a8d169cf001e42c9ba3947d7aa71148c","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"94cbe7e36ef94467b673ac872e89ad01","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"d3e97c6ef17241ba84065f6858a6f4e4","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"78290f474dbe465c9f8dc32b61ac85e3","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"04af1e37e3ef445d9e9f166383df42a4","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"0c2d5faea455416b982ac2b1f4951f65","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"48ab54beca4844618ace3428c8ee2ca0","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"78f3e89d0f24451db9fe58153d675dc0","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"5f76714145ab466b8661a0d9393a4cd2","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"be97a6aef43249ba85d73247b7d9c558","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"3f9a9906dd66427080b36a2a8a508772","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"d18776efc9a343d2bff47659ba7d2456","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"96fb13f06f9044bea80d7b420aa8e755","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"34cd2f38fd034a39a8d56acddf586f94","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"63261b2ee3b9401ba4bfe4fadafa9872","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"119ab0e0a2fe4a6fa687e05c8d66fc00","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"5cef775c249c429bb06ed2738a095f51","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"bee397c2997a4308ac0bfe2f122f3462","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"b84c0dc39f42451eb44fbdc324b2ac9f","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"69849a0629b84278a3736faa0f110af6","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"6576855ac23f47919a0f3e2c8d026f70","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"ace449c98d24441cb5e72a0e6fb347ce","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"e6b3b988b5db461886c13c07a2578d92","version_major":2,"version_minor":0}
```
:::

::: {.output .display_data}
``` json
{"model_id":"e0e15151b8d546af9218dbfafa546e59","version_major":2,"version_minor":0}
```
:::
:::

::: {.cell .markdown}
### Save Results
:::

::: {.cell .code execution_count="49"}
``` python
with open('shapleyVoteResults.pkl', 'wb') as f:
    pickle.dump(simulations, f)
```
:::

::: {.cell .code execution_count="50"}
``` python
# with open('shapleyVoteResults_Apr30_11-55PM.pkl', 'rb') as fp:
#     simulations = pickle.load(fp)
```
:::

::: {.cell .markdown}
### Analyze Trials
:::

::: {.cell .code execution_count="51"}
``` python
def meanAbsError(qshaps, cshaps):
    err = 0
    for qshap, cshap in zip(qshaps, cshaps):
        err += abs(qshap-cshap)
    return err
```
:::

::: {.cell .code execution_count="52"}
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

::: {.output .display_data}
![](figures/546.png)
:::
:::

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


 
