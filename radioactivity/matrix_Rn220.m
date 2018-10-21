% Matrix method
% A quick way to solve equations related to radioactive decay chains.
%
% For the mathematical background have a look at:
% L. Moral and A. F. Pacheco, Am. J. Phys. 71(7) 2003, 684 - 686. 
% DOI: 10.1119/1.1571834
%
% The Rn-220 (Thoron) decay chain:
% 
% For a classic solution with 'dsolve' have a look at bateman_radon220.m
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
% 
% In the manuscript the matrix containing all lambdas is named [A]. 
% A in a radioactive decay context usually refers to the
% activities. So [A] has been labeled [AM] to avoid confusion.

clear all;

% -----------------------------------------------------------------------
% values to change
l1 = log(2)/55.6;         % l for Rn-220 /s; l=log(2)/t(1/2)
l2 = log(2)/0.15;         % l for Po-216 /s; l=log(2)/t(1/2)
l3 = log(2)/(10.6*60*60); % l for Pb-212 /s; l=log(2)/t(1/2)
l4 = log(2)/(60.6*60);    % l for Bi-212 /s; l=log(2)/t(1/2)
l5 = log(2)/2.99e-7;      % l for Po-212 /s; l=log(2)/t(1/2)
l6 = log(2)/(3.05*60);    % l for Tl-208 /s; l=log(2)/t(1/2)

r41 = 0.64;               % Bi-212 branching to Po-212
r42 = 0.36;               % Bi-212 branching to Tl-208

N0 = 1;                   % The script calculates relative values. 
                          % So N0 can be any number > 0.
% -----------------------------------------------------------------------

% The lambda matrix is a NxN matrix (N = number of active species),
% which describes decay and formation of species in the decay chain.
% If you formulate the differential equations of the decay chain you can 
% easily see how the matrix is build up.
%
% rows of the lambda matrix [AM]
dN1 = [-l1   0   0   0       0   0]; % dN1(t)/dt = -l1*N1; decay of 1
dN2 = [ l1 -l2   0   0       0   0]; % dN2(t)/dt =  l1*N1(t)-l2*N2(t); formation of 2 from 1 and decay of 2
dN3 = [  0  l2 -l3   0       0   0]; % dN3(t)/dt =  l2*N2(t)-l3*N3(t); formation of 3 from 2 and decay of 3
dN4 = [  0   0  l3 -l4       0   0]; % dN4(t)/dt =  l3*N3(t)-l4*N4(t); formation of 4 from 3 and decay of 4
dN5 = [  0   0   0  r41*l4 -l5   0]; % dN5(t)/dt =  r41*l4*N4(t)-l5*N5(t); formation of 5 from 4 * branching and decay of 5
dN6 = [  0   0   0  r42*l4   0 -l6]; % dN6(t)/dt =  r42*l4*N4(t)-l6*N6(t); formation of 6 from 4 * branching and decay of 6
%
% the lambda matrix [AM]
AM = [dN1; dN2; dN3; dN4; dN5; dN6];

% the N0 vector
N0_vec = [N0; 0; 0; 0; 0; 0];

% AM_vec = eigenvector of [AM], AM_eval = eigenvalue of [AM]
[AM_vec,AM_eval]=eig(AM);

% create a diagonal matrix [L] with exp(AM_eval*t) which is exp(-lamda*t)
syms t; 
L = expm(AM_eval*t);

% [N] = eigenvector of [AM] * [L] * eigenvector of [AM]^-1^ * [N0]
N = AM_vec*L*inv(AM_vec)*N0_vec;

% multiply N * l for activities A; A = N * l 
% l is in the eigenvalue of AM, but in the wrong order and with the wrong sign
% take only values from the diagonale and reorder with 'flipud'
A = N.*flipud(diag(-AM_eval));

% divide by A0 (A0 = N0 * l1) to get relative activities
Arel = A/(N0*l1);

% sum of all activities
A_sum_rel=sum(Arel);

% plot section
hold on;
fplot(Arel(1),'Color','blue');
fplot(Arel(2),'Color','red');
fplot(Arel(3),'Color','green');
fplot(Arel(4),'Color','yellow');
fplot(Arel(5),'Color','magenta');
fplot(Arel(6),'Color','black');
fplot(A_sum_rel,'Color','black','LineStyle','--');
hold off;

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
