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

**Figure 1.** We generated 32 random weighted voting games for each condition. For each case, we generated random weights $w_j\in \mathbb{N}$ such that $q\leq \sum_j w_j < 2q$. There were three primary scenarios: (1) Four players, voting threshold $q=8$; (2) Eight players, voting threshold $q=16$; and (3) 12 players, voting threshold $q=32$. For each scenario, we approximated every player's Shapley value with our quantum algorithm using various $\ell$'s, assuming $\epsilon=0$. Next, we compared each approximated Shapley value to it's true value using absolute error. Finally, we took the reciprocal of the mean for all Shapley value errors in each random game for each condition.

### 3.2 One-Equation Model

The <a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">matlab scripts under this folder</a> can
be used to construct the following quantum system: $B^{\pm} (H^{\otimes n}\otimes I) \vert 0 \rangle^{\otimes n+1}$.

Repeatedly measure the rightmost qubit. The difference of the expected values of the systems is $\frac{\Phi_i}{V_{max}}$.

#### Quantum Version of $\gamma(n,m)$ and $V(S)$

##### Unitaries general form

<img src="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/blob/main/figures/matlab1.png" width=50%>

where

* $\phi^{\pm}(i,n) = \gamma(n,c(i))  \cdot V^{\pm}(i)$
* $c(i)$ is the number of ones in the binary representation of $i$.

``` matlab
gamma = @(n,m) 1 / ( nchoosek(n,m)*(n+1) );
c = @(i) nnz(dec2bin(i)-'0');
```

Example assuming $n+1=3, F = {0, 1, 2}$ and $i=1$

``` matlab
V_empty = 0;
V_0 = 0;
V_0_1  = 1;
V_0_2  = 1;
V_1 = 0;
V_1_2 = 0;
V_2 = 0;
V_0_1_2 = 1;
```

``` {verbatim}
W(S)
```

``` matlab
W_0 = V_0 - V_empty;
W_0_1 = V_0_1 - V_1;
W_0_2 = V_0_2 - V_2
```

``` {verbatim}
W_0_2 = 1
```

``` matlab
W_0_1_2 = V_0_1_2 - V_1_2;
W = [ W_0 W_0_1 W_0_2 W_0_1_2 ];
Wmax = max(W);
phi = @(i,n) gamma(n,c(i)) * W(i+1)/Wmax;
n = 2;
B = [...
    sqrt(1-phi(0,n)) sqrt(phi(0,n))   0 0 0 0 0 0;
    sqrt(phi(0,n))  -sqrt(1-phi(0,n)) 0 0 0 0 0 0;
    0 0 sqrt(1-phi(1,n)) sqrt(phi(1,n)) 0 0 0 0;
    0 0 sqrt(phi(1,n)) -sqrt(1-phi(1,n)) 0 0 0 0;
    0 0 0 0 sqrt(1-phi(2,n))  sqrt(phi(2,n)) 0 0;
    0 0 0 0 sqrt(phi(2,n)) -sqrt(1-phi(2,n)) 0 0;
    0 0 0 0 0 0 sqrt(1-phi(3,n))  sqrt(phi(3,n));
    0 0 0 0 0 0 sqrt(phi(3,n)) -sqrt(1-phi(3,n))
    ]
```

``` {verbatim}
B = 8×8
    1.0000         0         0         0         0         0         0         0
         0   -1.0000         0         0         0         0         0         0
         0         0    0.9129    0.4082         0         0         0         0
         0         0    0.4082   -0.9129         0         0         0         0
         0         0         0         0    0.9129    0.4082         0         0
         0         0         0         0    0.4082   -0.9129         0         0
         0         0         0         0         0         0    0.8165    0.5774
         0         0         0         0         0         0    0.5774   -0.8165
```

``` matlab
% check if unitary
isequal(B*B', eye(size(B*B',1)))
```

``` {verbatim}
ans = logical
   1
```

``` matlab
% create the input state
H = gate.qft(2);
I = gate.id(2)
```

``` {verbatim}
   (1,1)        1
   (2,2)        1
```

``` matlab
% check if unitary
isequal(Bminus*Bminus', eye(size(Bminus*Bminus',1)))
```

``` {verbatim}
ans =
   1
```


``` matlab
% create the input state
H = gate.qft(2);
I = gate.id(2)
```

``` {verbatim}
   (1,1)        1
   (2,2)        1
```

``` matlab
In = u_propagate(state('000'),tensor(tensor(H,H),I))
```


``` {verbatim}
In = +0.5 |000> +0.5 |010> +0.5 |100> +0.5 |110>
```

``` matlab
u_propagate(In, B)
```


``` {verbatim}
ans = +0.5 |000> +0.456435 |010> +0.204124 |011> +0.456435 |100> +0.204124 |101> +0.408248 |110> +0.288675 |111>
```

``` matlab
% samples
u = [];
for k=1:10000
    % apply the Shappley gate
    Ou = u_propagate(In, B);
    % measure the 3rd qubit
    [~,b,~] = measure(Ou, 3);
    cbit = b - 1;
    u = [ u cbit ];
end
fprintf('Shapply value is: %2.4f', 2^n * Wmax * mean(u) );
```

``` {verbatim}
Shapply value is:    0.6800
```

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
  year={2025},
  eprint={2412.14639},
  archivePrefix={arXiv},
  primaryClass={cs.CR},
  url={https://arxiv.org/abs/2412.14639}, 
}
```



