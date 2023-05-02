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
be used to construct the following quantum system: $B^{\pm} (H^{\otimes n}\otimes I) \vert 0 \rangle^{\otimes n+1}$.

Repeatedly measure the rightmost qubit. The difference of the expected values of the systems is $\frac{\Phi_i}{V_{max}}$.

#### Quantum Version of $\gamma(n,m)$ and $V(S)$

##### Unitaries general form

<img src="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/blob/main/figures/matlab1.png" width=35%>

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


### 3.2 Quantum Algorithm for Weighted Voting Games

The <a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">matlab scripts under this folder</a> can also
be used to compute Shapley values for weighted voting games.

#### 3.2.1 Register Preparation

Let us prepare the following register:

<img src="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/blob/main/figures/matlab2.png" width=35%>

``` matlab
%partition functions
t = @(ell,k) ( sin( (k*pi) / 2^(ell+1) ) )^2;
w = @(ell,k) t(ell,k+1) - t(ell,k);
ell = 1;
for k=0:2^ell-1
    w(ell,k)
end
```

``` {verbatim}
ans = 0.5000
ans = 0.5000
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



