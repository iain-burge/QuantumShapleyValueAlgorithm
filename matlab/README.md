# Quantum Algorithms for Shapley Value Calculation

## 1. Authors

<a href="https://github.com/iain-burge/iain-burge">Iain Burge</a>

<a href="https://carleton.ca/scs/people/michel-barbeau/">Michel Barbeau</a>

<a href="http://www-public.imtbs-tsp.eu/~garcia_a/web/">Joaquin Garcia-Alfaro</a>

## 2. Available Resources

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/python/">Python Code</a>.

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">Matlab Code</a>.

## 3. Sample Code and Results

### 3.1 One-Equation Model

The <a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">matlab scripts under this folder</a> can
be used to construct the following example, consisting of the two quantun systems (may be combined in one): $B^{\pm} (H^{\otimes n}\otimes I) \vert 0 \rangle^{\otimes n+1}$

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

Example assuming $n+1=3, F = {1, 2, 3}$ and $i=1$

``` matlab
V_0 = 0;
V_1 = 0;
V_1_2  = 1;
V_1_3  = 1;
V_2 = 0;
V_2_3 = 0;
V_3 = 0;
V_1_2_3 = 1;
```

$V^{+}()$ and $\phi^{-}()$

``` matlab
Vm = [V_0 V_3 V_2 V_2_3 ];
phim = @(i,n) gamma(n,c(i)) * Vm(i);
% value function upper bound
Vmax = max([Vp Vm]);
n = 2;
Bplus = 1/Vmax * [...
    sqrt(1-phip(1,n)) sqrt(phip(1,n)) 0 0 0 0 0 0;
    sqrt(phip(1,n)) -sqrt(1-phip(1,n)) 0 0 0 0 0 0;
    0 0 sqrt(1-phip(2,n)) sqrt(phip(2,n)) 0 0 0 0;
    0 0 sqrt(phip(2,n)) -sqrt(1-phip(2,n)) 0 0 0 0;
    0 0 0 0 sqrt(1-phip(2,n))  sqrt(phip(2,n)) 0 0;
    0 0 0 0 sqrt(phip(2,n)) -sqrt(1-phip(2,n)) 0 0;
    0 0 0 0 0 0 sqrt(1-phip(2,n))  sqrt(phip(2,n));
    0 0 0 0 0 0 sqrt(phip(2,n)) -sqrt(1-phip(2,n))
    ]
```

``` {verbatim}
Bplus = 8×8
    1.0000         0         0         0         0         0         0         0
         0   -1.0000         0         0         0         0         0         0
         0         0    0.9129    0.4082         0         0         0         0
         0         0    0.4082   -0.9129         0         0         0         0
         0         0         0         0    0.9129    0.4082         0         0
         0         0         0         0    0.4082   -0.9129         0         0
         0         0         0         0         0         0    0.9129    0.4082
         0         0         0         0         0         0    0.4082   -0.9129

```

``` matlab
% check if unitary
isequal(Bplus*Bplus', eye(size(Bplus*Bplus',1)))
```

``` {verbatim}
ans =
   1
```

``` matlab
Bminus = 1/Vmax * [...
    sqrt(1-phim(1,n)) sqrt(phim(1,n)) 0 0 0 0 0 0;
    sqrt(phim(1,n)) -sqrt(1-phim(1,n)) 0 0 0 0 0 0;
    0 0 sqrt(1-phim(2,n)) sqrt(phim(2,n)) 0 0 0 0;
    0 0 sqrt(phim(2,n)) -sqrt(1-phim(2,n)) 0 0 0 0;
    0 0 0 0 sqrt(1-phim(2,n))  sqrt(phim(2,n)) 0 0;
    0 0 0 0 sqrt(phim(2,n)) -sqrt(1-phim(2,n)) 0 0;
    0 0 0 0 0 0 sqrt(1-phim(2,n))  sqrt(phim(2,n));
    0 0 0 0 0 0 sqrt(phim(2,n)) -sqrt(1-phim(2,n))
    ]
```

``` {verbatim}
Bminus = 8×8
     1     0     0     0     0     0     0     0
     0    -1     0     0     0     0     0     0
     0     0     1     0     0     0     0     0
     0     0     0    -1     0     0     0     0
     0     0     0     0     1     0     0     0
     0     0     0     0     0    -1     0     0
     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0    -1

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
% samples
sp = []; sm = [];
for k=1:50
    % apply the Shappley gate
    Ou = u_propagate(In, Bplus);
    % measure the 3rd qubit
    [~,b,~] = measure(Ou, 3);
    cbit = b - 1;
    sp = [ sp cbit ];
    % apply the Shappley gate
    Ou = u_propagate(In, Bminus);
    % measure the 3rd qubit
    [~,b,~] = measure(Ou, 3);
    cbit = b - 1;
    sm = [ sm cbit ];
end
fprintf('Shapply value is: %6.1f', Vmax * (mean(sp) - mean(sm)));
```

``` {verbatim}
Shapply value is:    0.1
```

## References

If using this code for research purposes, please cite:

Iain Burge, Michel Barbeau and Joaquin Garcia-Alfaro. Quantum Algorithms for Shapley Value Calculation. *To appear*. May 2023.

```
@inproceedings{burge-barbeau-alfaro2023Shapley,
  title={Quantum Algorithms for Shapley Value Calculation},
  author={Burge, Iain and Barbeau, Michel and Garcia-Alfaro, Joaquin},
  booktitle={To appear},
  pages={1--9},
  year={2023},
  month={May},
}
```



