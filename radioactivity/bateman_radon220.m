% Rn-220 (Thoron) decay chain
%
% For the basics have a look at bateman.m.
%
%
%                                           (r41=64%)-l4-> Po-212 
%                                          /                     \l5
% Rn-220 -l1-> Po-216 -l2-> Pb-212 -l3-> Bi-212                   ----> Pb-208 (stable)
%                                          \                     /l6
%                                           (r42=36%)-l4-> Tl-208 
% 1: Rn-220 t(1/2) = 55.6 s
% 2: Po-216 t(1/2) = 0.15 s
% 3: Pb-212 t(1/2) = 10.6 h
% 4: Bi-212 t(1/2) = 60.6 min
% 5: Po-212 t(1/2) = 299 ns
% 6: Tl-208 t(1/2) = 3.05 min
%
% l = lamda = decay constant
% t(1/2) = half-life 
% N0 is the initial quantity of substance that will decay.
% N(t) is the quantity that still remains and has not yet decayed after a time t.
% dN/dt = -l*N
% N(t) = N0*exp(-l*t)
% l = log(2)/t(1/2) (log means the natural logarithm)
% A = activity; is the number of decays per unit time of a radioactive sample
% A = l*N
% A0 = l1*N0
% rxy = branching ratio 
% et = elution time

clear all;

syms N1(t) N2(t) N3(t) N4(t) N5(t) N6(t) l1 l2 l3 l4 l5 l6 r41 r42 N0 A0;

eq_1 = diff(N1(t),t) == -l1*N1(t);              % decay of Rn-220        
eq_2 = diff(N2(t),t) == -l2*N2(t)+l1*N1(t);     % decay of Po-216 and formation from Rn-220
eq_3 = diff(N3(t),t) == -l3*N3(t)+l2*N2(t);     % decay of Pb-212 and formation from Po-216
eq_4 = diff(N4(t),t) == -l4*N4(t)+l3*N3(t);     % decay of Bi-212 and formation from Pb-212
eq_5 = diff(N5(t),t) == -l5*N5(t)+r41*l4*N4(t); % decay of Po-212 and formation from Bi-212 
eq_6 = diff(N6(t),t) == -l6*N6(t)+r42*l4*N4(t); % decay of Tl-208 and formation from Bi-212 

% solve the system of diff. equations
% conditions N1(t=0) = N0; N2 to N5(t=0) = 0

sol = dsolve ([eq_1 eq_2 eq_3 eq_4 eq_5 eq_6, N1(0)==N0 N2(0)==0 N3(0)==0 N4(0)==0 N5(0)==0 N6(0)==0]); 

A1 = sol.N1*l1; % transform N to A
A2 = sol.N2*l2; % transform N to A
A3 = sol.N3*l3; % transform N to A
A4 = sol.N4*l4; % transform N to A
A5 = sol.N5*l5; % transform N to A
A6 = sol.N6*l6; % transform N to A

A1 = subs(A1,[N0*l1],[A0]); % substitute N0*l1 with A0
A2 = subs(A2,[N0*l1],[A0]); % substitute N0*l1 with A0
A3 = subs(A3,[N0*l1],[A0]); % substitute N0*l1 with A0
A4 = subs(A4,[N0*l1],[A0]); % substitute N0*l1 with A0
A5 = subs(A5,[N0*l1],[A0]); % substitute N0*l1 with A0
A6 = subs(A6,[N0*l1],[A0]); % substitute N0*l1 with A0

A1 = A1 / A0; % divide by A0 to get relative activities
A2 = A2 / A0; % divide by A0 to get relative activities
A3 = A3 / A0; % divide by A0 to get relative activities
A4 = A4 / A0; % divide by A0 to get relative activities
A5 = A5 / A0; % divide by A0 to get relative activities
A6 = A6 / A0; % divide by A0 to get relative activities

A2 = simplify(A2,'Steps',40); % simplify and divide by A0 to get relative activities
A3 = simplify(A3,'Steps',40); % simplify and divide by A0 to get relative activities
A4 = simplify(A4,'Steps',40); % simplify and divide by A0 to get relative activities
A5 = simplify(A5,'Steps',40); % simplify and divide by A0 to get relative activities
A6 = simplify(A6,'Steps',40); % simplify and divide by A0 to get relative activities

% -----------------------------------------------------------------------
% values to change
l1_v = log(2)/55.6;         % l for Rn-220 /s;  l=log(2)/t(1/2)
l2_v = log(2)/0.15;         % l for Po-216 /s; l=log(2)/t(1/2)
l3_v = log(2)/(10.6*60*60); % l for Pb-212 /s; l=log(2)/t(1/2)
l4_v = log(2)/(60.6*60);    % l for Bi-212 /s;  l=log(2)/t(1/2)
l5_v = log(2)/2.99e-7;      % l for Po-212 /s; l=log(2)/t(1/2)
l6_v = log(2)/(3.05*60);    % l for Tl-208 /s; l=log(2)/t(1/2)

r41_v = 0.64;               % Bi-212 branching to Po-212
r42_v = 0.36;               % Bi-212 branching to Tl-208
% -----------------------------------------------------------------------

% symbolic substitution of constants with "real" values
A1s = subs(A1, [l1,l2,l3,l4,l5,l6,r41,r42], [l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,r41_v,r42_v]);
A2s = subs(A2, [l1,l2,l3,l4,l5,l6,r41,r42], [l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,r41_v,r42_v]);
A3s = subs(A3, [l1,l2,l3,l4,l5,l6,r41,r42], [l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,r41_v,r42_v]);
A4s = subs(A4, [l1,l2,l3,l4,l5,l6,r41,r42], [l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,r41_v,r42_v]);
A5s = subs(A5, [l1,l2,l3,l4,l5,l6,r41,r42], [l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,r41_v,r42_v]);
A6s = subs(A6, [l1,l2,l3,l4,l5,l6,r41,r42], [l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,r41_v,r42_v]);
A_sum = A1s + A2s + A3s + A4s + A5s + A6s; % sum of all activities

% plot section
hold on;
fplot(A1s,'Color','blue');
fplot(A2s,'Color','red');
fplot(A3s,'Color','green');
fplot(A4s,'Color','yellow');
fplot(A5s,'Color','magenta');
fplot(A6s,'Color','black');
fplot(A_sum,'Color','black','LineStyle','--');
hold off

% plot options
ax=gca;
ax.Title.String = {'Decay of Rn-220'};
ax.YLim=[1e-15 3];
ax.XLim=[1e-4 5e6];
ax.YScale='log';
ax.XScale='log';
ax.XLabel.String='t /s';
ax.YLabel.String='rel. activity A/A_0';
grid on;
grid minor;
legend('Rn-220','Po-216','Pb-212','Bi-212','Po-212','Tl-208','A_{sum}','Location','best');
