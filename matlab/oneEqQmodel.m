%% *One-equation quantum model* 
% Construct the two quantum systems (may be combined in one):
% $$B^{\pm} (H^{\otimes n}\otimes I) \vert 0 \rangle^{\otimes n+1}$$
% 
% Repeatedly measure the rightmost qubit. The difference of the expected values 
% of the systems is $\frac{\Phi_i}{V_{max}}$.
%% Quantum Version of $\gamma \left(n,m\right)$ and$V(S)$
% Unitaries general form
% 
% 
% $$B^\pm(n) = \bigoplus_{i=0}^{2^n-1} \frac{1}{V_{max}}\pmatrix{ \sqrt{1-\phi^\pm(i,n)} 
% &   \sqrt{\phi^\pm(i,n)}  \cr\sqrt{\phi^\pm(i,n)} & -\sqrt{1-\phi^\pm(i,n)}   
% }$$
% 
% where
%% 
% * $\phi^{\pm}(i,n)=\gamma(n,c(i)) \cdotV^{\pm}(i)$
% * $c\left(i\right)$ is the number of ones in the binary representation of 
% $i$.

gamma = @(n,m) 1 / ( nchoosek(n,m)*(n+1) );
c = @(i) nnz(dec2bin(i)-'0');
%% 
% Example assuming $n+1=3, F = {1, 2, 3}$and $i=1$

V_0 = 0; 
V_1 = 0;
V_1_2  = 1;
V_1_3  = 1;
V_2 = 0; 
V_2_3 = 0;
V_3 = 0; 
V_1_2_3 = 1;
%% 
% $V^{+}()$ and $\phi^{+}()$

Vp = [V_1 V_1_3 V_1_2 V_1_2_3 ];
phip = @(i,n) gamma(n,c(i)) * Vp(i);
%% 
% $V^{-}()$ and $\phi^{-}()$

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
% check if unitary
isequal(Bplus*Bplus', eye(size(Bplus*Bplus',1)))
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
% check if unitary
isequal(Bminus*Bminus', eye(size(Bminus*Bminus',1)))

% create the input state
H = gate.qft(2);
I = gate.id(2)
In = u_propagate(state('000'),tensor(tensor(H,H),I))
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
%%