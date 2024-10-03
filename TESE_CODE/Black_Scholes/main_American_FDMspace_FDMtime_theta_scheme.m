clear all
close all
clc;

figure_increment = 0;

format long

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
Msi = round((Smax-Smin)/(hS*0.5^0) - 1);
hS = (Smax-Smin)/(Msi+1);

mesh1D = mesh_line1D(Smin,Smax,Msi);

S = Smin:hS:Smax; % S = Smin+dS*(0:1:Ms);
mm = 0:1:Msi+1;

%% Time
%Nt=100; % Number of time points
%ht = T/(Nt+1);

ht=0.1;
Nt = round((T-0)/(ht*(0.5^0)) - 1);
ht = (T-0)/(Nt+1);


t = 0:ht:T; %time vector
nn = 0:1:Nt+1;


%% Initializing the matrix of the option values
%V = initial_BS_CallOption_American(r,sigma,d,T,K,S,t);
V = initial_BS_PutOption_American(r,sigma,d,T,K,S,t);


%% Obstaculo: Valor final de V
% Final conditions prescribed by the Call payoff at expiry: V(S,T)=max(E-S,0)
g=V(:,end);
gr=g(2:end-1);


G = sparse(Msi+2,Nt+2);
for n=Nt+2:-1:1
    G(:,n) = g;
end





%% Matriz da discretização no espaço do operador eliptico

L = sparse(Msi,Msi+2);

Am = S/(2*hS).*(sigma^2*S./hS - r);
Bm =  -r - sigma^2*(S.^2)./(hS^2);
Cm = S/(2*hS).*(sigma^2*S./hS + r);


for m = 1:Msi
    L(m,m) = Am(m+1);
    L(m,m+1) = Bm(m+1);
    L(m,m+2) = Cm(m+1);
end

Lr = L(:,2:end-1);
R = sparse(Msi,Msi+2);

R(1,1) = L(1,1); R(end,end) = L(end,end);

I = speye(Msi,Msi);


%% Determining the matrix coeficients of the algorithm: theta-scheme

%theta = 0; % Euler Explicito (Forward-Euler): condição para o esquema explicito: sigma^2*S.^2*ht/(hS^2)<=1
theta = 1; % Euler Implicito (Backward-Euler)
%theta = 0.5; % Crank-Nicolson

A = theta*Lr - (1/ht)*I;
B = (1-theta)*Lr + (1/ht)*I;


tic
for n=Nt+1:-1:1
    c = (1-theta)*R*V(:,n+1) + theta*R*V(:,n);
    %V(2:Msi+1,n)= (-A)\(B*V(2:Msi+1,n+1)+c);
    %V(2:Msi+1,n)= max( (-A)\(B*V(2:Msi+1,n+1)+c) , gr);
    V(2:Msi+1,n) = SSNM_iterative_min(A,-(B*V(2:Msi+1,n+1)+c),gr); 
    n
end
toc
%% Plots

% Figure of the Value of the European Option, V(S,t)
figure_increment =figure_increment+1;
figure(figure_increment)
%contourf(S,t,V')
%ss = surf(S,t,full(V)')
%set(ss,'LineStyle','none')
mesh(S,t,full(V'))
%colormap jet
%colorbar
%title('American Option value surface, V(S,t), within the FEM Method')
xlabel('S','FontSize',20);
ylabel('t','FontSize',20);
zlabel('V(S,t)','FontSize',20);
xlim([0 15])

% Figure of the value of the Europena Option, V(S,t), as a function of S
% at three different times:t=0, T/2 and T (expiry).
figure_increment =figure_increment+1;
figure(figure_increment)
plot(S,V(:,1)','r-',S,V(:,round(Nt/2))','g-',S,V(:,end)','b-','LineWidth',3);
xlabel('S','FontSize',20);
ylabel('V(S,t)','FontSize',20);
%title('American Option within the FEM Method, at three different times: t=0, T/2 and T (expiry)');
legend('t=0','t=T/2','t=T');
xlim([0 15])
grid;


ii = find(V-G==0); % indices onde a diferença é zero, ou seja V=G
UU =  sparse(Msi+2,Nt+2);
UU(ii)=1; % marcar esses pontos
WW = sparse(Msi+2,Nt+2);
WW(2:end,:) = UU(1:end-1,:) - UU(2:end,:); % diferenciar para obter uma curva
WW(end,:) = 0; % eliminar a fronteira
WW(end-1,:) = 0; % eliminar a fronteira
WW(end-2,:) = 0; % eliminar a fronteira

S_0=K;
[mm,ii] = min(abs(S-K));
for n=1:length(t)
    WW(ii+1:end,n) = 0; % eliminar a fronteira
end 

[px,py] = find(WW==1);
[iipp] = find(WW==1);

figure_increment = figure_increment + 1;
figure(figure_increment)

% Criar a superfície com cor mais transparente
ss = surf(S, t, full(UU)');
set(ss, 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Ajustar a transparência e remover as bordas
colormap jet

hold on

% Traçar a linha
p = plot3(S(px), t(py), full(V(iipp)), 'black', 'LineWidth', 5);

% Adicionar marcadores de cor para a legenda usando scatter3
h1 = scatter3(NaN, NaN, NaN, 'or', 'filled', 'MarkerFaceColor', 'red'); % Para a legenda do valor 1
h2 = scatter3(NaN, NaN, NaN, 'ob', 'filled', 'MarkerFaceColor', 'blue'); % Para a legenda do valor 0

% Adicionar a legenda
legend([p, h1, h2], {'S_f(t)', 'Stop', 'Hold'}, 'Location', 'Best','FontSize',20);

% Configurar os eixos e os rótulos
xlabel('S', 'FontSize', 20);
ylabel('t', 'FontSize', 20);
xlim([0 15]);

hold off




%{
figure_increment =figure_increment+1;
figure(figure_increment)
%contourf(S,t,V')
%ss = surf(S,t,full(UU)')
%set(ss,'LineStyle','none')
%colormap jet
%colorbar
mesh(S,t,full(UU)')
%mesh(S,t,full(WW)')
hold on
plot3(S(px),t(py),full(V(iipp)),'black','LineWidth',5);
%title('Regiao onde V=g')
xlabel('S','FontSize',20);
ylabel('t','FontSize',20);
xlim([0 15])

figure_increment =figure_increment+1;
figure(figure_increment)
%contourf(S,t,V')
%ss = surf(S,t,full(UU)')
%set(ss,'LineStyle','none')
colormap jet
colorbar
%mesh(S,t,full(UU)')
mesh(S,t,full(WW)')
%title('Regiao onde V=g')
xlabel('S','FontSize',20);
ylabel('t','FontSize',20);
xlim([0 20])


% Figure of the Value of the European Option, V(S,t)
figure_increment =figure_increment+1;
figure(figure_increment)
%contourf(S,t,V')
%ss = surf(S,t,full(V)')
%set(ss,'LineStyle','none')
mesh(S,t,full(V'))
hold on
plot3(S(px),t(py),full(V(iipp)),'o','LineWidth',1);
%colormap jet
%colorbar
%title('American Option value surface, V(S,t), within the FEM Method')
xlabel('S','FontSize',20);
ylabel('t','FontSize',20);
zlabel('V(S,t)','FontSize',20);
xlim([0 15])

figure_increment =figure_increment+1;
figure(figure_increment)
plot3(S(px),t(py),full(V(iipp)),'LineWidth',1);
xlabel('S','FontSize',20);
ylabel('t','FontSize',20);
zlabel('V(S,t)','FontSize',20);
xlim([0 15])
%}



%% VALORES V(S_0,T_0)



S_0=6;
[mm,ii] = min(abs(S-S_0));
V_S0_6 = full(V(ii,1))

S_0=7;
[mm,ii] = min(abs(S-S_0));
V_S0_7 = full(V(ii,1))

S_0=8;
[mm,ii] = min(abs(S-S_0));
V_S0_8 = full(V(ii,1))

S_0=10;
[mm,ii] = min(abs(S-S_0));
V_S0_10 = full(V(ii,1))


S_0=12;
[mm,ii] = min(abs(S-S_0));
V_S0_12 = full(V(ii,1))

S_0=14;
[mm,ii] = min(abs(S-S_0));
V_S0_14 = full(V(ii,1))

S_0=16;
[mm,ii] = min(abs(S-S_0));
V_S0_16 = full(V(ii,1))
