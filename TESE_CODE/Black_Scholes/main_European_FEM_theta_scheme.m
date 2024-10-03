clear all
close all
clc;

figure_increment = 0;

format shortE

%% Constantes do problema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.06; % Interest rate
sigma=0.3; % Volatility of the underlying
%d=0.3; % Continuous dividend yield
d=0; 
T=1; % Maturation (expiry)of contract
K=10; % Exercise price of the underlying (K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Espaço (Price), e a sua descretização

Smax=200; % Maximum share price considered
Smin=0; % Minimum share price considered

%Msi=400; % Number of share price points (interior)
%hS = (Smax-Smin)/(Msi+1);


hS=0.5;
Msi = (Smax-Smin)/(hS*0.5^5) - 1;
hS = (Smax-Smin)/(Msi+1);

mesh1D = mesh_line1D(Smin,Smax,Msi);

S = Smin:hS:Smax; % S = Smin+dS*(0:1:Ms);
mm = 0:1:Msi+1;

%% Time
%Nt=100; % Number of time points
%ht = T/(Nt+1);

ht=0.1;
Nt = round((T-0)/(ht*0.5^0) - 1);
ht = (T-0)/(Nt+1);


t = 0:ht:T; %time vector
nn = 0:1:Nt+1;



%% Initializing the matrix of the option values
%V = initial_BS_CallOption_European(r,sigma,d,T,K,S,t);
V = initial_BS_PutOption_European(r,sigma,d,T,K,S,t);

%% Matrizes (Space/Price)
% assemble constant matrices/vectors
xpow2_ux_vx_fun=@(u,v,ux,vx,x,h) x.^2.*ux.*vx;
%B_xpow2_ux_vy=bilinear_assembly(xpow2_ux_vx_fun,mesh);
B_xpow2_ux_vy=Bilinear_assembly_xpow2_ux_vx(S);

x_ux_v_fun=@(u,v,ux,vx,x,h) x.*ux.*v;
%B_x_ux_v=bilinear_assembly(x_ux_v_fun,mesh);
B_x_ux_v=Bilinear_assembly_x_ux_v(S);

u_v_fun=@(u,v,ux,vx,x,h) u.*v;
%B_u_v=bilinear_assembly(u_v_fun,mesh);
B_u_v=Bilinear_assembly_u_v(S);

ux_vx_fun=@(u,v,ux,vx,x,h) ux.*vx;
%B_ux_vx=bilinear_assembly(ux_vx_fun,mesh);
B_ux_vx=Bilinear_assembly_ux_vx(S);

h0 = S(2)-S(1);
hMp1 = S(end)-S(end-1);

% integral na fronteira:
B_front = sparse(Msi+2,Msi+2);
B_front(1,1) = S(1)^2/h0; B_front(1,2) = -S(1)^2/h0;
B_front(end,end-1) = -S(end)^2/hMp1; B_front(end,end) = S(end)^2/hMp1;

% matriz do operador bilinear:
B = -1/2*sigma^2 * B_xpow2_ux_vy +... 
    + (r-sigma^2)* B_x_ux_v +...
    -r* B_u_v +...
    + 0.5*sigma^2*B_front;


% matriz do produto do dual:
Pd = B_u_v;

%% Determining the matrix coeficients of the algorithm: theta-scheme
% Implementing the algorithm
%U'(t) = (Un_plus_1-Un)/ht

%theta = 0; % Euler Explicito (Forward-Euler): condição para o esquema explicito: sigma^2*S.^2*ht/(hS^2)<=1
theta = 1; % Euler Implicito (Backward-Euler)
%theta = 0.5; % Crank-Nicolson

cond_euler_exp = max(sigma^2*S.^2*ht/(hS^2))




E = theta*B - (1/ht)*Pd;
H = (1-theta)*B + (1/ht)*Pd;

E_R = E(2:end-1,:);
E_r = E(2:end-1,2:end-1);

H_R = H(2:end-1,:);



for n=Nt+1:-1:1
    V(2:Msi+1,n)= (-E_r)\(E_R*V(:,n) + H_R*V(:,n+1));
end


%% Analitic Solution for European Put Option with Discretization

% Initializing the matrix of the option value
Vp_ana = [];
Vc_ana = [];

for i=1:length(S);

d1 = (log(S(i)/K) + (r+0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
d2 = (log(S(i)/K) + (r-0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));

%N_d1 = normcdf(d1);
%N_d2 = normcdf(d2);

Vc_ana(i,:)= S(i) .* normcdf(d1)  - K * exp( -r.*(T-t) ) .* normcdf(d2);
Vp_ana(i,:)=- S(i) .* normcdf(-d1) + K * exp( -r.*(T-t) ) .* normcdf(-d2);


end


%% Erro
% haverá sempre um valor NaN na solução analitica, temos de descartar tal
% ocorrencia para calcular o erro
V_min_Vp_ana = reshape(Vp_ana,1,[]) - reshape(full(V),1,[]);
V_min_Vp_ana = V_min_Vp_ana(~isnan(V_min_Vp_ana));

V_min_Vc_ana = Vc_ana - full(V);
V_min_Vc_ana = V_min_Vc_ana(~isnan(V_min_Vc_ana));

e_p = norm(V_min_Vp_ana,"inf")
e_c = norm(V_min_Vc_ana,"inf")


S_0=8;

[mm,ii] = min(abs(S-S_0))

Valor_real_S0=Vp_ana(ii,1)

erro=abs(full(V(ii,1))-Vp_ana(ii,1))

%{

%% Plots

% Figure of the Value of the European Option, V(S,t)
figure_increment =figure_increment+1;
figure(figure_increment)
%contourf(S,t,V')
hh = surf(S,t,full(V)');
set(hh,'LineStyle','none')
%mesh(S,t,full(V(1:length(S),1:length(t))'))
colormap jet
colorbar
title('European Option value surface, V(S,t), within the FEM Method')
xlabel('S')
ylabel('t')


% Figure of the value of the Europena Option, V(S,t), as a function of S
% at three different times:t=0, T/2 and T (expiry).
figure_increment =figure_increment+1;
figure(figure_increment)
plot(S,V(:,1)','r-',S,V(:,round(Nt/2))','g-',S,V(:,end)','b-','LineWidth',1.5);
xlabel('S');
ylabel('V(S,t)');
title('European Option within the FEM Method, at three different times: t=0, T/2 and T (expiry)');
legend('t=0','t=T/2','t=T');
grid;

% Solucao analitica
figure_increment =figure_increment+1;
figure(figure_increment)
hh = surf(S,t,Vp_ana');
set(hh,'LineStyle','none')
title('Analitic European Option solution, V(S,t)')
xlabel('S')
ylabel('t')

figure_increment =figure_increment+1;
figure(figure_increment)
plot(S,Vp_ana(:,1)','r-',S,Vp_ana(:,round(Nt/2))','g-',S,Vp_ana(:,end)','b-','LineWidth',1.5);
xlabel('S');
ylabel('V(S,t)');
title('Analitic European Option solution, V(S,t), at three different times: t=0, T/2 and T (expiry)');
legend('t=0','t=T/2','t=T');
grid;


%% Plot Erro

err_dif = abs(full(V)-Vp_ana);
figure_increment =figure_increment+1;
figure(figure_increment)
mesh(S,t,err_dif')
title('Erro')
xlabel('S')
ylabel('t')

%}
