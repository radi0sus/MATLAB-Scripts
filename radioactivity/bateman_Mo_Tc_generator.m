% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Characteristics of a Tc-99m generator
% -------------------------------------------------------
% Mo-99 (t(1/2))=66 h) -> Tc-99m (t(1/2)=6 h)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
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
% r12 = branching ratio (only 86% of Mo-99 decays to Tc-99m)
% et = elution time
% A Tc-99m quantity will be eluted after "et".
%
% This script was made for autodidactic purposes and may contain errors. 
% DON'T USE IT FOR MEDICAL OR OTHER APPLICATIONS!

clear all;

syms N1(t) N2(t) l1 l2 N0 A0 et r12;

eq_1 = diff(N1(t),t) == -l1*N1(t);                 % decay of 1
eq_2 = diff(N2(t),t) == -l2*N2(t)+r12*l1*N1(t);    % decay of 2 + formation of 2 from 1 times branching ratio

% solve the system of diff. equations
% conditions N1(t=0) = N0; N2(t=et) = 0
sol = dsolve ([eq_1 eq_2, N1(0) == N0 N2(et) == 0]); 

A1 = sol.N1*l1; % transform N to A
A2 = sol.N2*l2; % transform N to A

A1 = subs(A1,[N0*l1],[A0]); % substitute N0*l1 with A0
A2 = subs(A2,[N0*l1],[A0]); % substitute N0*l1 with A0

% -----------------------------------------------------------------------
% values to change
l1_v = log(2)/65.94; % l for Mo-99;  l=log(2)/t(1/2)
l2_v = log(2)/6.01;  % l for Tc-99m; l=log(2)/t(1/2)
A0_v = 10;           % Mo-99 activity at t=0 in MBq or KBq or Bq....
r12_v = 0.86;        % branching ratio: only 86% Mo-99 decays to Tc-99m

et_max = 240;        % elution time maximum /h
et_int = 12;         % elution time interval /h
et_s = 0;            % start first elution after ... h 
et_p = 0;            % pause /h
% -----------------------------------------------------------------------

% symbolic substitution of constants with "real" values
A1s = subs(A1, [l1,l2,A0], [l1_v,l2_v,A0_v]);
A2s = subs(A2, [l1,l2,A0,r12], [l1_v,l2_v,A0_v,r12_v]);

% transform A2s to plot functions
f1 = symfun(A2s,et);      % plot activity from et to t
f2 = symfun(A2s,[et t]);  % obtain maximum Tc-99m activity at et + t

% plot section
hold on;
n1=fplot(A1s);   % plot the decay of Mo-99
n2=fplot(f1(0)); % plot the decay of Tc-99m without elution, et = 0

n1.Color='b';
n2.Color='g';

% plot decay of Tc-99m after each elution interval
if et_int+et_p == 0 
    % elution starts after et_s with no further elution of Tc-99m
    n3=fplot(f1(et_s));
    n3.Color='r';
else
    % elution starts after et_s with further elution of Tc-99m
    for et=et_s-et_p:et_int+et_p:et_max
        n3=fplot(f1(et+et_p),[0 et+et_int+et_p]);
        n3.Color='r';
        
        st=stem(et+et_int+et_p,f2(et+et_p,et+et_int+et_p),'-');
        st.Marker='none';
        st.Color='r';
    end
end
hold off;

% plot options
ax2=gca;
ax2.Title.String = {'Tc-99m generator','Mo-99 \rightarrow Tc-99m'};
%ax2.YScale='log';
ax2.YLim=[0 A0_v];
ax2.YTick=[0:1:A0_v];
ax2.XLim=[0 240];
ax2.XTick=[0:12:240];
ax2.XTickLabelRotation = 90;
ax2.XMinorTick='on';
%grid on;
%grid minor;
ax2.XLabel.String='t /h';
ax2.YLabel.String='Activity';
legend('Mo-99','Tc-99m (no elution)', ['Tc-99m, start elution after ' num2str(et_s) ' h, elution every ' num2str(et_int) ' h, pause ' num2str(et_p) ' h']);
