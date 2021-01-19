%Code containing all methods of imposing inner boundary conditions

%The SBP operators used,d2_2,d2_4,d2_6,d2_8 and d2_10, are available in:
%https://github.com/scicompuu/sbplib

syms a c x t mu
u = -a*tanh(a*(x-c*t)/(2*mu)) + c; %Analytical solution used for outer boundary data
u = matlabFunction(u);
ux = matlabFunction(diff(u,x));

clear a c x t mu;

mu1 = 0.2; % Viscosity in first domain
mu2 = 0.1; % Viscosity in second domain
mu = (mu1+mu2)/2; %Average viscosity across both domains

%Constants for analytical solution
c = 1; 
a = 1;

u = @(t,x) u(a,c,mu,t,x);
ux = @(t,x) ux(a,c,mu,t,x);

%Left and right boundary
x_l = -1;
x_r = 1;


%Number of data points
m = 21;
% m = 41;
% m = 81;
% m = 161;
% m = 321;

%Skew symmetric splitting gives alpha
alpha = 1/3;

%Dicretizing domain
x_i = (x_r+x_l)/2;
x1 = linspace(x_l,x_i,m);
x2 = linspace(x_i,x_r,m);
x = [x1 x2];
h = x(2)-x(1);
T = 0.2; %Simulation time
t = 0; %Starting time
dt = 0.5*h^2; %Time step
count = T/dt; %Steps taken 

[H, HI, D1, D2, e_1, e_m, M, Q, S_1, S_m] = d2_6(m,h); %SBP Operators

w_sol = u(0,x');%Stacked solutionvector 

alpha1 = H(m,m);
alpha2 = H(1,1); 

%Penalty terms
%Outer boundary, Robin
taul = -1;
taur = 1;

%Inner boundary, Dirichlet 
taud1 = @(w_sol) 1/3*w_sol(m);
taud2= @(w_sol) -1/3*w_sol(m+1);

sigd1 = -1/(4*(alpha1+alpha2)); 
sigd2 = sigd1;

%Inner boundary, Neumann
beta1 = -alpha2/(alpha1+alpha2); 
beta2= beta1 + 1;

%Stacked operators 
DD1 = sparse([D1 zeros(m);zeros(m) D1]);
DD2 = sparse([mu1*D2 zeros(m); zeros(m) mu2*D2]);
HItot = sparse([HI zeros(m);zeros(m) HI]);
I = sparse(diag(ones(m,1)));
Itot = sparse([I zeros(m); zeros(m) I]);

a_l = @(u,t) (u + abs(u))/6;
a_r = @(u,t) (u - abs(u))/6;

g_l = @(t) a_l(u(t,x_l),t)*u(t,x_l) - mu1*ux(t,x_l); %Left outer boundary data
g_r = @(t) a_r(u(t,x_r),t)*u(t,x_r) - mu2*ux(t,x_r); %Right outer boundary data


%SAT Outer boundary
SATY = @(w_sol,t) [HI*taul*e_1*((a_l(w_sol(1),t)*e_1'-mu1*S_1')*w_sol(1:m) - g_l(t));
              HI*taur*e_m*((a_r(w_sol(end),t)*e_m'-mu2*S_m')*w_sol(m+1:end)- g_r(t))];

%SAT Inner Dirichlet condition          
SATIDir1 = @(w_sol,t)[HI*taud1(w_sol)*e_m*(w_sol(m)- w_sol(m+1));
                     HI*taud2(w_sol)*e_1*(w_sol(m+1) - w_sol(m))];
SATIDir2 = @(w_sol,t)[HI*sigd1*e_m*mu1*(w_sol(m) - w_sol(m+1)); 
                      HI*sigd2*e_1*mu2*(w_sol(m+1)- w_sol(m))];                 
%SAT Inner Neumann condition 
SATINeu = @(w_sol,t) [HI*beta1*e_m*(mu1*S_m'*w_sol(1:m) - mu2*S_1'*w_sol(m+1:end)) ;
                      HI*beta2*e_1*(mu2*S_1'*w_sol(m+1:end) - mu1*S_m'*w_sol(1:m))];  

%Projection Inre Neumann randvillkor + Dirichlet                  
L = [e_m' -e_1'; mu1*S_m' -mu2*S_1']; % Use if using projection method

%Projection Inre Dirichlet
% L = [e_m' -e_1']; %Use if using hybrid method

P = Itot - HItot*L'*((L*HItot*L')^-1)*L; %Projection operator

%Projection Damping term
sigP = @(w_sol) -norm(w_sol,inf)./h ;

%Burgers equation                
A = @(w_sol,t) -(alpha)*diag(w_sol)*DD1 - ((1-alpha)/2)*DD1*diag(w_sol) + DD2;

%Choose Method of imposing inner boundary conditions
%Projection
wt = @(w_sol,t,x) P*A(P*w_sol,t)*P*w_sol + P*SATY(P*w_sol,t)+...
    sigP(P*w_sol)*(Itot-P)*w_sol;

%SAT
% wt = @(w_sol,t,x) A(w_sol,t)*w_sol + SATY(w_sol,t)+...
%     SATIDir1(w_sol,t)+SATIDir2(w_sol,t)+ SATINeu(w_sol,t);

%Hybrid
% wt = @(w_sol,t,x) P*A(P*w_sol,t)*P*w_sol + P*SATY(P*w_sol,t)+...
%      P*SATINeu(P*w_sol,t) + sigP(P*w_sol)*(Itot-P)*w_sol;

%Runge Kutta used to time step
for f=1:count
    
    k1 = dt * wt(w_sol,t);
    k2 = dt * wt(w_sol+0.5*k1,t+0.5*dt);
    k3 = dt * wt(w_sol+0.5*k2,t+0.5*dt);
    k4 = dt * wt(w_sol+k3,t+dt);
    
    w_sol = w_sol + 1/6*(k1+2*k2+2*k3+k4);
    
    t = t + dt;
    %To plot solution in real time
    plot(x,w_sol,'*'); 
    plot(x,w_sol,'LineWidth', 2)
    hold on
    xline(0,'--r');
    ylim([0 2]);
    hold off
    drawnow;
end

%Plotting final solution
plot(x,w_sol,'LineWidth', 2)
hold on
xline(0,'--r');
ylim([0 2]);
hold off
