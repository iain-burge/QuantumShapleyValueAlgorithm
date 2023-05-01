%% Quantum Shapley Value Calculation
% *Author*: Michel Barbeau
% 
% Version: April 26, 2023

% cd to QIT directory
clear
cd './qit'
% run the init script
init
%% *One-equation quantum model (non-efficient, made efficient with Iain's ideas)*
% Construct the two quantum systems (may be combined in one):
% $$B^{\pm} (H^{\otimes n}\otimes I) \vert 0 \rangle^{\otimes n+1}$$
% 
% Repeatedly measure the rightmost qubit. The difference of the expected values 
% of the systems is \frac{\Phi_i}{V_{max}}.
%% Quantum Version of $\gamma \left(n,m\right)$ and$V(S)$
% Unitaries general form
% 
% 
% $$B^\pm(n) = \bigoplus_{i=0}^{2^n-1} \pmatrix{ \sqrt{1-\phi^\pm(i,n)} &   
% \sqrt{\phi^\pm(i,n)}  \cr\sqrt{\phi^\pm(i,n)} & -\sqrt{1-\phi^\pm(i,n)}   }$$
% 
% where
%% 
% * $\phi^{\pm}(i,n)=\gamma(n,c(i)) \cdotV^{\pm}(i)$
% * $c\left(i\right)$ is the number of ones in the binary representation of 
% $i$.

gamma = @(n,m) 1 / ( nchoosek(n,m)*(n+1) );
c = @(i) nnz(dec2bin(i)-'0');
%% 
% Example assuming $n+1=3$, $F = \{0, 1, 2\}$ and $i=1$

V_empty = 0; 
V_0 = 0;
V_0_1  = 1;
V_0_2  = 1;
V_1 = 0; 
V_1_2 = 0;
V_2 = 0; 
V_0_1_2 = 1;
%% 
% $$W(S)$$

W_0 = V_0 - V_empty;
W_0_1 = V_0_1 - V_1;
W_0_2 = V_0_2 - V_2
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
% check if unitary
isequal(B*B', eye(size(B*B',1)))
% create the input state
H = gate.qft(2);
I = gate.id(2)
In = u_propagate(state('000'),tensor(tensor(H,H),I))
u_propagate(In, B)
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
%% Quantum Algorithm for Weighted Voting Games
% Prepartion of partition register
% 
% $$W_{\ell} H^{\otimes 2^{\ell}} \vert 0 \rangle^{\otimes 2^{\ell}}$$
% 
% $$W_\ell = \bigoplus_{i=0}^{2^n-1} \pmatrix{ \sqrt{W_{\ell}(i))} &   \sqrt{W_{\ell}(i+1)}  
% \cr\sqrt{W_{\ell}(i+1)} & -\sqrt{W_{\ell}(i)}   }$$

%partition functions
t = @(ell,k) ( sin( (k*pi) / 2^(ell+1) ) )^2;
w = @(ell,k) t(ell,k+1) - t(ell,k); 
ell = 1;
for k=0:2^ell-1
    w(ell,k)
end
%% 
% Controlled rotations
% 
% 


%% Reference
% Burge, I., Barbeau, M., & Garcia-Alfaro, J. (2023). A Quantum Algorithm for 
% Shapley Value Estimation. _arXiv preprint arXiv:2301.04727_.
% 
%