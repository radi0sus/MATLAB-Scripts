% Calculation of radioactive decay chains (mother -> daugther -> granddaughter .....): 
% 1 -> 2 -> 3 -> 4 -> .....using the Bateman equation.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Transient equilibrium: t(1/2)_mother > t(1/2)_daugther
% -------------------------------------------------------
% Mo-99 (t(1/2))=66 h) -> Tc-99m (t(1/2)=6 h)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% l = lamda = decay constant
% t(1/2) = half-life (units (e.g. seconds, days, years), should be equal for all t(1/2)
% N0 is the initial quantity of substance that will decay.
% Since the script calculates relative values, changes of N0 have no effect on the result. 
% N(t) is the quantity that still remains and has not yet decayed after a time t.
% dN/dt = -l*N
% N(t) = N0*exp(-l*t)
% l = log(2)/t(1/2) (log means the natural logarithm)
% A = activity; is the number of decays per unit time of a radioactive sample
% A = l*N
% A0 = l1*N0
%
% normal scale and log scale plots 

clear all;

syms N1(t) N2(t) l1 l2 N0;

eq_1 = diff(N1(t),t) == -l1*N1(t);                  % decay of 1
eq_2 = diff(N2(t),t) == -l2*N2(t)+l1*N1(t);         % decay of 2 + formation of 2 from 1

% solve the system of diff. equations
% conditions N1(t=0) = N0; N2(t=0) = 0
sol = dsolve ([eq_1 eq_2, N1(0) == N0 N2(0) == 0]); 

sol.N1 = sol.N1*l1/(N0*l1); % transform N to A and make it relative "/(N0*l1)"
sol.N2 = sol.N2*l2/(N0*l1); % transform N to A and make it relative "/(N0*l1)"

sum_act = sol.N1+sol.N2;    % sum of all activities A

% disp(sol.N1);
% disp(sol.N2);
% disp(sum_act);

l1_v = log(2)/66; % insert t(1/2) for 1, l=log(2)/t(1/2)
l2_v = log(2)/6;  % insert t(1/2) for 2, l=log(2)/t(1/2)
N0_v = 1;         % Since the script calculates relative values, 
                  % changes of N0 have no effect on the result. 

% symbolic substitution of constants with "real" values
sol.N1=subs(sol.N1, [l1,l2,N0], [l1_v,l2_v,N0_v]);
sol.N2=subs(sol.N2, [l1,l2,N0], [l1_v,l2_v,N0_v]);
sum_act=subs(sum_act, [l1,l2,N0], [l1_v,l2_v,N0_v]);

max_n2=solve(diff(sol.N2),t);   % time to reach maximum activity; 1st derivate
max_sum=solve(diff(sum_act),t); % time to reach maximum overall activity; 1st derivate

%plot section

fig=figure;
fig.Name = 'Bateman equation';

% normal scale plot
ax1=subplot(2,1,1);

hold on;

n1=fplot(sol.N1);   % rel. activitity for 1
n2=fplot(sol.N2);   % rel. activitity for 2
sum=fplot(sum_act); % sum of the 1 + 2

n1.Color='r';
n2.Color='g';
sum.Color='b';

% set mark for time to reach maximum activity of 2
st=stem(max_n2,subs(sol.N2,t,max_n2),'--');
st.Marker='none';
st.Color='g';

% set mark for time to reach maximum activity of 1 + 2
st=stem(max_sum,subs(sum_act,t,max_sum),'--');
st.Marker='none';
st.Color='b';

%plot options
ax1.Title.String = {'Transient equilibrium','Mo-99 \rightarrow Tc-99m'};
ax1.YLim=[0.01 2];
ax1.XLim=[0 200];
ax1.XLabel.String='t /h';
ax1.YLabel.String='rel. Activity (A/A0)';
legend('Mo (t(1/2) = 66 h)','Tc (t(1/2) = 6 h)','sum',[num2str(double(max_n2)) ' h'],[num2str(double(max_sum)) ' h']);

hold off;

% log scale plot
ax2=subplot(2,1,2);

hold on;

n1=fplot(sol.N1);   % rel. activitity for 1
n2=fplot(sol.N2);   % rel. activitity for 2
sum=fplot(sum_act); % sum of the 1 + 2

n1.Color='r';
n2.Color='g';
sum.Color='b';

% set mark for time to reach maximum activity of 2
st=stem(max_n2,subs(sol.N2,t,max_n2),'--');
st.Marker='none';
st.Color='g';

% set mark for time to reach maximum activity of 1 + 2
st=stem(max_sum,subs(sum_act,t,max_sum),'--');
st.Marker='none';
st.Color='b';

%plot options
ax2.Title.String = {'Transient equilibrium','Mo-99 \rightarrow Tc-99m'};
ax2.YScale='log';
ax2.YLim=[0.01 2];
ax2.XLim=[0 200];
ax2.XLabel.String='t /h';
ax2.YLabel.String='rel. Activity (A/A0)';
legend('Mo (t(1/2) = 66 h)','Tc (t(1/2) = 6 h)','sum',[num2str(double(max_n2)) ' h'],[num2str(double(max_sum)) ' h']);

hold off;
