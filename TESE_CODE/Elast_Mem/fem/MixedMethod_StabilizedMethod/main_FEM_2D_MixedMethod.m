close all;
clear all;
clear;
clc;

figure_increment = 0;

format short
format shortE

%% funções (anonimas) 
f_anony_fun=@(x,y)(-1);
%obstaculo
bg = 0.9;
c1 = -bg/sqrt(1-bg^2);
c2 = sqrt(1-bg^2)-c1*bg;
a = 0.8294;
gr_anony_fun=@(r) sqrt(1-r.^2).*(r<bg)+ ...
     (c1.*r + c2).*(r>=bg);
g_anony_fun=@(x,y) sqrt(1-sqrt(x.^2+y.^2).^2).*(sqrt(x.^2+y.^2)<bg)+ ...
    + (c1.*sqrt(x.^2+y.^2) + c2).*(sqrt(x.^2+y.^2)>=bg);
gx_anony_fun=@(x,y) (-2*x.*(sqrt(1-sqrt(x.^2+y.^2).^2))).*(sqrt(x.^2+y.^2)<bg)+ ...
    + (c1.*x./sqrt(x.^2+y.^2) ).*(sqrt(x.^2+y.^2)>=bg);
gy_anony_fun=@(x,y) (-2*y.*(sqrt(1-sqrt(x.^2+y.^2).^2))).*(sqrt(x.^2+y.^2)<bg)+ ...
    + (c1.*y./sqrt(x.^2+y.^2) ).*(sqrt(x.^2+y.^2)>=bg);

ur_anony_fun=@(r) gr_anony_fun(r).*(r<=a)+((-a^2 + 4*gr_anony_fun(a) + r.^2)*log(2) + (4-r.^2)*log(a)+(a^2-4-4*gr_anony_func(a)).*log(r))./(4*(log(2) - log(a))).*(r>a);
u_ana_anony_fun=@(x,y) gr_anony_fun(sqrt(x.^2+y.^2)).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((-a^2 + 4*gr_anony_fun(a) + sqrt(x.^2+y.^2).^2)*log(2) + (4-sqrt(x.^2+y.^2).^2)*log(a)+(a^2-4-4*gr_anony_fun(a)).*log(sqrt(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);

ux_ana_anony_fun=@(x,y) gx_anony_fun(x,y).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((2*x)*log(2) + (-2*x)*log(a)+(a^2-4-4*gr_anony_fun(a)).*(x./(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);
uy_ana_anony_fun=@(x,y) gy_anony_fun(x,y).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((2*y)*log(2) + (-2*y)*log(a)+(a^2-4-4*gr_anony_fun(a)).*(y./(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);




%% Dominio e Mesh
Xc = [0;0];
R = 2;
%NN = 2.3; %h=0.8
NN = 5; %h=0.4
%NN = 10; %h=0.2
%NN = 22; %h=0.1
%NN = 45; %h=0.05
%NN = 95; %h=0.025
mesh = Circular_Mesh_with_baricenter(NN,R,Xc);
%mesh = Mesh_by_Tom();


Np = length(mesh.p);
Nv = length(mesh.p) - length(mesh.ba);
Nt = max(size(mesh.t));
 

figure_increment = figure_increment+1;
figure(figure_increment)
TR = triangulation(mesh.t',mesh.p(:,1:Nv)');
triplot(TR)
title('Mesh triangular')
xlabel('x')
ylabel('y')


%Plot apenas dos pontos da mesh
figure_increment = figure_increment+1;
figure(figure_increment)
%plot boundary:
triplot(TR)
hold on
plot(mesh.p(1,mesh.b)',mesh.p(2,mesh.b)','or','MarkerSize',8);
hold on
%plot interior:
plot(mesh.p(1,mesh.in)',mesh.p(2,mesh.in)','.k','MarkerSize',16);
hold on
%plot baricenter:
plot(mesh.p(1,mesh.ba)',mesh.p(2,mesh.ba)','.g','MarkerSize',8);
xlabel('x')
ylabel('y')
legend('Mesh Connections','Boundary Points','Interior Points','baricenter')



figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',g_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('Obstacle');


%% Matrizes 
% assemble constant matrices/vectors

grad_u_dot_grad_v=@(u,v,ux,uy,vx,vy,x,y,lamb,h) ux.*vx+uy.*vy;
A=bilininear_assembly_VhVh(grad_u_dot_grad_v,mesh);

spy(A)

lamb_v=@(u,v,ux,uy,vx,vy,x,y,lamb,h) lamb.*v;
B=bilininear_assembly_VhP0(lamb_v,mesh);

fv=@(u,v,ux,uy,vx,vy,x,y,lamb,h) f_anony_fun(x,y).*v;
f=linear_assembly_Vh(fv,mesh);

g_lamb=@(u,v,ux,uy,vx,vy,x,y,lamb,h) g_anony_fun(x,y).*lamb;
g=linear_assembly_P0(g_lamb,mesh);



%% Algoritmo

U_ana = u_ana_anony_fun(mesh.p(1,[mesh.b , mesh.in])',mesh.p(2,[mesh.b , mesh.in])');

% constantes do algoritmo:
alpha=1e-1;
%alpha = 5
c=1e-2;

% Com este algoritmo queremos descobrir 'Lambda' e 'U', para isso então imaginemos um vector 'Y=[U , Lambda]'.
% 'U' terá Np pontos e Lamb terá Nt pontos:
indx_U=1:Np;
indx_U_b=mesh.b; %boundary/fronteira
indx_U_in=setdiff(1:Np,mesh.b); %interior
indx_Lamb=(1:Nt)+Np; % indx_Lamb=Np+1:Nt+Np; indx_Lamb=indx_U+Np;


U=zeros(Np,1);
Lamb=zeros(Nt,1);
Y=zeros(Np+Nt,1);

Lamb=zeros(Nt,1);
indx_Lamb_pos = find(Lamb>0) % inicialmente vai ser o conjunto vazio
indx_Y=[indx_U_in Np+indx_Lamb_pos];

% -> Consrução do sistema matricial
% Temos o seguinte sistema: A*U - B^T*Lamb = f
% Com a seguinte condição (desigualdade): -B*U - alpah*C*Lamb + g_alpah <= 0 , onde g_alpah = g-alpha*gC
% Lamb^T*(B*U + C*Lamb - g_alpha) = 0
% Lamb >= 0

% Assumindo que a condição: -B*U - alpah*C*Lamb + g_alpah <= 0, está ativa
% ou seja existe algumas linhas desse sistema de desigualdade que se
% verifica a igualdade: 
% -B*U - alpah*C*Lamb + g_alpah = 0 , e Lamb >= 0
% ou está inativa: 
% -B*U - alpah*C*Lamb + g_alpah < 0 , e Lamb = 0

AA=[A -B'; -B sparse(Nt,Nt)];
bb=[f;-g];
tic
for itr=1:100
%display(['Iteration: ',num2str(itr)]);

% 1º Quais linhas de Y podemos contar, para resolver o nosso sistema (de igualdades). 
% queremos calcular o interior de U (SEMPRE), e as entradas que Lambda>0 (QUE VARIARÃO), ou seja:
indx_Y=[indx_U_in Np+indx_Lamb_pos];
% 2º Calcular Y a partir do sistema (de onde se tira U):
Y(indx_Y)=AA(indx_Y,indx_Y)\bb(indx_Y);
U_prev = U;
U=Y(indx_U); % guardar a solução U:

% ou, U pode ser ainda de outra maneira:
%a = indx_Lamb_pos;
%U(indx_U_in) = ( A(indx_U_in,indx_U_in)+B(a,indx_U_in)'*[(alpha*C(a,a))\B(a,indx_U_in)] )\( f(indx_U_in)+B(a,indx_U_in)'*[(alpha*C(a,a))\(g(a)-alpha*gC(a))] );
    
%Achar o novo Lambda (O novo Lambda NÂO É A SOLUÇÂO DO SISTEMA!!!!)
Lamb_prev=Lamb; % guardar Lambda anterior
indx_Lamb_pos_prev = indx_Lamb_pos; % guardar index's positivos do Lambda anterior

% Semismooth Newton method:
Lamb=Y(indx_Lamb);
%Lamb=max(Lamb+c*(g-alpha*gC-(B*U+alpha*C*Lamb)),zeros(Nt,1));
% or:
%Lamb = max((alpha*C)\(g-alpha*gC - B*U),zeros(Nt,1));

%indx_Lamb_pos=find(Lamb+c*(g-alpha*gC-(B*U+alpha*C*Lamb))>0)'; % Semismooth Newton method
indx_Lamb_pos=find(Lamb+c*(g - B*U)>0)';

if isequal(sort(indx_Lamb_pos_prev),sort(indx_Lamb_pos))
    if norm(Lamb_prev-Lamb)<1e-32
    %if norm(U-U_prev,inf)<1e-32
        display('Primal-dual strategy converged!');
        display(['Iteration: ',num2str(itr)]);
        break
    end
end
end
toc

%Agora o verdadeiro Lambda:
Lamb = max(Lamb,zeros(Nt,1));


%% Erro
uh = full(U());

[e_H1,e_L2,e_inf,H1K,L2K] = Compute_error_Mixed(uh,mesh,u_ana_anony_fun,ux_ana_anony_fun,uy_ana_anony_fun);

mesh
DOF = length(indx_U_in) + 2*length(mesh.t)
e_L2
e_H1
%e_inf
% Agora 


U(mesh.ba) = ( U(mesh.t(1,:)) + U(mesh.t(2,:)) + U(mesh.t(3,:)) )/3;
e_inf = norm(full(U([mesh.b , mesh.in]))-U_ana,Inf)



%% Plots
figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',u_ana_anony_fun(mesh.p(1,:)',mesh.p(2,:)'));
%zlim([-0.3 1])
title('Analitic Solution');

figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',full(U)); % 'full' para sair de 'sparse'
%zlim([-0.3 1])
title('Numeric Solution');


%{
figure_increment = figure_increment+1;
figure(figure_increment)
trisurf(mesh.t',mesh.p(1,:)',mesh.p(2,:)',full(U)); % 'full' para sair de 'sparse'
%zlim([-0.3 1])
hold on
scatter3(mesh.p(1,mesh.ba)',mesh.p(2,mesh.ba)',full(U(mesh.ba)),'.r')
%}


%{
figure_increment = figure_increment+1;
figure(figure_increment)
for i=1:length(Lamb)

    % Coordenadas dos vértices do triângulo    
    x = [mesh.p(1,mesh.t(1,i)) , mesh.p(1,mesh.t(2,i)) , mesh.p(1,mesh.t(3,i))];
    y = [mesh.p(2,mesh.t(1,i)) , mesh.p(2,mesh.t(2,i)) , mesh.p(2,mesh.t(3,i))];
    z = [Lamb(i) , Lamb(i) , Lamb(i)];


    % Plotar o triângulo preenchido
    fill3(x, y, z, 'r'); % 'r' indica a cor vermelha para o preenchimento
    hold on;
    % Plotar os vértices do triângulo
    %scatter3(x, y, z, 'filled', 'k'); % 'k' indica a cor preta para os pontos
end 
xlabel('x')
ylabel('y')
zlabel('u')
%}




figure_increment = figure_increment + 1;
figure(figure_increment)
colormap(jet); % Escolher um colormap, por exemplo, 'jet'
cmap = colormap; % Obter o colormap atual
num_colors = size(cmap, 1); % Número de cores no colormap

% Normalizar Lamb para usar como índice no colormap
Lamb_norm = (Lamb - min(Lamb)) / (max(Lamb) - min(Lamb));

for i = 1:length(Lamb)
    % Coordenadas dos vértices do triângulo    
    x = [mesh.p(1, mesh.t(1, i)), mesh.p(1, mesh.t(2, i)), mesh.p(1, mesh.t(3, i))];
    y = [mesh.p(2, mesh.t(1, i)), mesh.p(2, mesh.t(2, i)), mesh.p(2, mesh.t(3, i))];
    z = [Lamb(i), Lamb(i), Lamb(i)];
    
    % Determinar a cor baseada no valor normalizado de Lambda
    color_idx = round(Lamb_norm(i) * (num_colors - 1)) + 1;
    color = cmap(color_idx, :);
    
    % Plotar o triângulo preenchido
    fill3(x, y, z, color, 'FaceColor', color, 'EdgeColor', 'none'); % Modificação importante
    hold on;
end

teta=-pi:0.01:pi;
x=a*cos(teta);
y=a*sin(teta);
plot3(x,y,zeros(1,numel(x)),'r','linewidth',1.1)

xlabel('x', 'FontSize', 40);
ylabel('y', 'FontSize', 40);
%zlabel('u', 'FontSize', 40);

% Adicionar uma barra de cores para referência e remover os números
c = colorbar;
c.Ticks = [];

hold off;
