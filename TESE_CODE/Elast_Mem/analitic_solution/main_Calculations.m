close all;
clear all;
clear;
clc;

figure_increment = 0;

%% calculo dos parâmetros da função g(r)
% r = sqrt(x.^2+y.^2)
bg = 0.9; % picewise r point 
c1 = -bg/sqrt(1-bg^2);
c2 = sqrt(1-bg^2)-c1*bg;

%% funções (anonimas) de g [obstaculo]

gr_anony_fun=@(r) sqrt(1-r.^2).*(r<bg)+ ...
    + (c1.*r + c2).*(r>=bg);
gr_line_anony_fun=@(r) (-r)./sqrt(1-r.^2).*(r<bg)+ ...
    + (c1).*(r>=bg);
g_anony_fun=@(x,y) sqrt(1-sqrt(x.^2+y.^2).^2).*(sqrt(x.^2+y.^2)<bg)+ ...
    + (c1.*sqrt(x.^2+y.^2) + c2).*(sqrt(x.^2+y.^2)>=bg);
% derivadas em x e y:
gx_anony_fun=@(x,y) (-2*x.*(sqrt(1-sqrt(x.^2+y.^2).^2))).*(sqrt(x.^2+y.^2)<bg)+ ...
    + (c1.*x./sqrt(x.^2+y.^2) ).*(sqrt(x.^2+y.^2)>=bg);
gy_anony_fun=@(x,y) (-2*y.*(sqrt(1-sqrt(x.^2+y.^2).^2))).*(sqrt(x.^2+y.^2)<bg)+ ...
    + (c1.*y./sqrt(x.^2+y.^2) ).*(sqrt(x.^2+y.^2)>=bg);

%% Descobrir ponto 'a'
% assumindo que existe uma região r in (a,2), onde (g-u)<0 <=> g<u
% se (g-u)<0 então temos Lamb=0
% e ficamos com uma ODE para resolver: u_rr-(1/r)*u_r=1
% cuja solução: u(r) = A+B*log(r)+r^2/4  ; u'(r)=B*(1/r)+r/2

% onde temos as seguintes condições (para decobrir A e B [em função de 'a']):
% u(a)=g(a) <=> A+B*log(a)+a^2/4 = g(a)
% u(2)=0  <=> A+B*log(2)+a^2/4 = 0

% onde se fica com um sistema:
%A+B*log(a) = g(a)-a^2/4
%A+B*log(2) = -1

% agora tiramos A e B, como solução de um sistema:

syms a g_a g_lin_a log2 r

M = [1 , log(a);
     1 , log2 ];
b = [g_a-a^2/4 ; -1];

A_and_B = M^-1*b
A = log(a)/(log2 - log(a)) + (log2*(- a^2/4 + g_a))/(log2 - log(a));
B = - 1/(log2 - log(a)) - (- a^2/4 + g_a)/(log2 - log(a));

u = A + B*log(r) + r^2/4;
u_a_to_2 = simplify(u)

% Finalmente para descobrir o valor de a, devemos calcular a raiz da
% seguinte igualdade: u'(a)=g'(a)
% tiremos o valor de 'a', onde  u'(a)-g'(a)=0

ff =@(a) (2*a*(log(2) - log(a))+ (a^2-4-4*gr_anony_fun(a))*(1/a))/(4*(log(2) - log(a))) - gr_line_anony_fun(a);
a = fzero(ff,0.95) % x = fzero(fun,x0) encontra x onde fun=0, fzero começa 
% em x0e tenta localizar um ponto x1 onde fun(x1) tenha o sinal oposto de 
% fun(x0). Em seguida fzero, reduz iterativamente o intervalo onde o sinal
% de fun muda para chegar a uma solução.

%% funções (anonimas) de u [solucao]
ur_anony_fun=@(r) gr_anony_fun(r).*(r<=a)+...
    + ((-a^2 + 4*gr_anony_fun(a) + r.^2)*log(2) + (4-r.^2)*log(a)+(a^2-4-4*gr_anony_fun(a)).*log(r))./(4*(log(2) - log(a))).*(r>a);

ur_line_anony_fun=@(r) gr_line_anony_fun(r).*(r<=a)+...
    + ((2*r)*log(2) + (-2*r)*log(a)+(a^2-4-4*gr_anony_fun(a)).*(1./r))./(4*(log(2) - log(a))).*(r>a);
%{
epsilon = 1e-10;
ur_line_anony_fun = @(r) (ur_anony_fun(r+epsilon) - ur_anony_fun(r)) / epsilon;
%}

u_ana_anony_fun=@(x,y) gr_anony_fun(sqrt(x.^2+y.^2)).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((-a^2 + 4*gr_anony_fun(a) + sqrt(x.^2+y.^2).^2)*log(2) + (4-sqrt(x.^2+y.^2).^2)*log(a)+(a^2-4-4*gr_anony_fun(a)).*log(sqrt(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);
ux_ana_anony_fun=@(x,y) gx_anony_fun(x,y).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((2*x)*log(2) + (-2*x)*log(a)+(a^2-4-4*gr_anony_fun(a)).*(x./(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);
uy_ana_anony_fun=@(x,y) gy_anony_fun(x,y).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((2*y)*log(2) + (-2*y)*log(a)+(a^2-4-4*gr_anony_fun(a)).*(y./(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);

%% Dominio e Mesh
Xc = [0;0];
R = 2;
NN = 7;
mesh = Circular_Mesh_RUI(NN,R,Xc);

figure_increment = figure_increment+1;
figure(figure_increment)
TR = triangulation(mesh.t',mesh.p')
triplot(TR)
title('Mesh triangular')
xlabel('x')
ylabel('y')

figure_increment = figure_increment+1;
figure(figure_increment)
%plot boundary:
triplot(TR)
hold on
plot(mesh.p(1,mesh.b)',mesh.p(2,mesh.b)','or','MarkerSize',8);
hold on
%plot interior:
plot(mesh.p(1,mesh.in)',mesh.p(2,mesh.in)','.k','MarkerSize',16);
xlabel('x')
ylabel('y')
legend('Mesh Connections','Boundary Points','Interior Points')


[xs,ys]=meshgrid(-2:0.1:2);


%% Plot do Obstaculo

r = linspace(0,2,10000);
gr = gr_anony_fun(r);
gr_line = gr_line_anony_fun(r);
ur = ur_anony_fun(r);
ur_line = ur_line_anony_fun(r);

figure_increment = figure_increment+1;
figure(figure_increment)
plot(r,gr)
grid 
title('g(r)')
xlabel('r')

figure_increment = figure_increment+1;
figure(figure_increment)
plot(r,gr_line)
grid 
title('g\_line(r)')
xlabel('r')


figure_increment = figure_increment+1;
figure(figure_increment)
surf(xs,ys,g_anony_fun(xs,ys));
title('g(x,y)')
xlabel('x')
ylabel('y')

figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',g_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('Obstacle: g(x,y)');
xlabel('x')
ylabel('y')

figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',gx_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('Obstacle: g_x(x,y)');
xlabel('x')
ylabel('y')

figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',gy_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('Obstacle: g_y(x,y)');
xlabel('x')
ylabel('y')

%% Plot do Obstaculo

figure_increment = figure_increment+1;
figure(figure_increment)
plot(r,ur)
grid 
title('u(r)')
xlabel('r')

figure_increment = figure_increment+1;
figure(figure_increment)
plot(r,ur_line)
grid 
title('u\_line(r)')
xlabel('r')


figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',u_ana_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('Analitic Solution: u(x,y)');
xlabel('x')
ylabel('y')

figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',ux_ana_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('u_x(x,y)');
xlabel('x')
ylabel('y')

figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',uy_ana_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('u_y(x,y)');
xlabel('x')
