% Calculation of branched decays (mother -> daughter1 + daughter2): 1 -> 2 + 3 
% using the Bateman equation.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
% K–Ar dating
% -------------------------------------------------------
%                                   /-(beta)-(br1_1(%): 89.28%) -> Ca-40 (stable)
% K-40 (t(1/2)=1.248*10^9 years)) --   
%                                   \-(ec)---(br1_2(%): 10.72%) -> Ar*-40 -(gamma)-> Ar-40(stable) 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% (Note that K-40 -(ec)-> Ar-40 (0.001%) branch is ignored.)
%
% l = lamda = decay constant
% t(1/2) = half-life (units (e.g. seconds, days, years), should be equal for all t(1/2)
% N0 is the initial quantity of substance that will decay.
% Since the script calculates relative values, changes of N0 have no effect on the result. 
% N(t) is the quantity that still remains and has not yet decayed after a time t.
% dN/dt = -l*N
% N(t) = N0*exp(-l*t)
% l = log(2)/t(1/2) (log means the natural logarithm)
% br = branching ratio

clear all;

syms N1(t) N2(t) N3(t) l1 l2 l3 N0 br1_1 br1_2;

eq_1 = diff(N1(t),t) == -l1*N1(t);                 % decay of 1
eq_2 = diff(N2(t),t) == -l2*N2(t)+br1_1*l1*N1(t);  % decay of 2 + formation of 2 from 1 times branching ratio
eq_3 = diff(N3(t),t) == -l3*N3(t)+br1_2*l1*N1(t);  % decay of 3 + formation of 3 from 1 times branching ratio

% solve the system of diff. equations
% conditions N1(t=0) = N0; N2(0) = 0; N3(0) = 0
sol = dsolve ([eq_1 eq_2 eq_3, N1(0) == N0 N2(0) == 0 N3(0) == 0]); 

% ----------------------------------------------------
% values to change
l1_v = log(2)/1.248e9;  % l for K-40;  l=log(2)/t(1/2)
l2_v = 0;               % 2 is a stable isotope
l3_v = 0;               % 3 is considered to be stable for the calculation
N0_v = 1000;            % number of active K-40 nuclei at t = 0
br1_1_v = 0.8928;       % branching ratio: K-40 -> Ca-40
br1_2_v = 0.1072;       % branching ratio: K-40 -> Ar*-40
% ----------------------------------------------------

% symbolic substitution of constants with "real" values
N1s=subs(sol.N1, [l1,l2,l3,N0], [l1_v,l2_v,l3_v,N0_v]);
N2s=subs(sol.N2, [l1,l2,l3,br1_1,br1_2,N0], [l1_v,l2_v,l3_v,br1_1_v,br1_2_v,N0_v]);
N3s=subs(sol.N3, [l1,l2,l3,br1_1,br1_2,N0], [l1_v,l2_v,l3_v,br1_1_v,br1_2_v,N0_v]);

%plot section
hold on;

n1=fplot(N1s);  % plot the decay of K-40
n2=fplot(N2s);  % plot the formation of Ca-40
n3=fplot(N3s);  % plot the formation of Ar*-40

n1.Color='r';
n2.Color='g';
n3.Color='b';

% The formula in text books for the amount of Ar is sometimes given as: 
% Ar-40 = 0.1072 * K-40 *(e^(l*t)-1).
% compare to N3s by uncommenting the following 3 lines
% K40(t) = N1s;
% Ar40(t) = 0.1072 * K40(t) *(exp((log(2)/1.248e9)*t)-1);
% n4=fplot(Ar40,[0 10e9]);

hold off;

%plot options
ax2=gca;
ax2.Title.String = {'K–Ar dating (branched decay)','K-40 \rightarrow Ca-40 (89.28%) + Ar-40 (10.72%)'};
%ax2.YScale='log';
ax2.YLim=[0 N0_v];
ax2.XLim=[0 10e9];
ax2.XLabel.String='t /a';
ax2.YLabel.String='N';
legend('N(K-40)','N(Ca-40)', 'N(Ar*-40)');
