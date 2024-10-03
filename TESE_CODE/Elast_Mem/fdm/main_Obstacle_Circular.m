close all;
clear 
clc

figure_increment = 0;

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

ur_anony_fun=@(r) gr_anony_fun(r).*(r<=a)+((-a^2 + 4*gr_anony_fun(a) + r.^2)*log(2) + (4-r.^2)*log(a)+(a^2-4-4*gr_anony_func(a)).*log(r))./(4*(log(2) - log(a))).*(r>a);

u_ana_anony_fun=@(x,y) gr_anony_fun(sqrt(x.^2+y.^2)).*(sqrt(x.^2+y.^2)<=a)+ ...
    + ((-a^2 + 4*gr_anony_fun(a) + sqrt(x.^2+y.^2).^2)*log(2) + (4-sqrt(x.^2+y.^2).^2)*log(a)+(a^2-4-4*gr_anony_fun(a)).*log(sqrt(x.^2+y.^2)))./(4*(log(2) - log(a))).*(sqrt(x.^2+y.^2)>a);

u_bound_anony_fun=@(x,y) 0;

% plot do obstaculo
figure_increment = figure_increment+1;
figure(figure_increment)
[xs,ys]=meshgrid(-2:0.1:2);
surf(xs,ys,g_anony_fun(xs,ys));

%% Domain and Grid
Xc = [0;0];
R = 2;
N = 20;

grid = Grid2D_Circular(N,R,Xc)

% plot grid
figure_increment = figure_increment+1;
figure(figure_increment)
%scatter(p_int(:,1),p_int(:,2),'.')
scatter(grid.p(grid.idx_p_ip,1),grid.p(grid.idx_p_ip,2),'.')
hold on
scatter(grid.p(grid.idx_p_bp,1),grid.p(grid.idx_p_bp,2),'o')
hold on
%scatter(GGrid(:,1),GGrid(:,2),'x')
hold off
title('Uniform Orthogonal Grid')

%% Analitic Solution plot
u_ana = u_ana_anony_fun(grid.p(:,1),grid.p(:,2));
figure_increment = figure_increment+1;
figure(figure_increment)
scatter3(grid.p(:,1),grid.p(:,2),u_ana)
title('Analitic Solution')

figure_increment = figure_increment+1;
figure(figure_increment)
scatter3(grid.p(:,1),grid.p(:,2),u_ana_anony_fun(grid.p(:,1),grid.p(:,2)))
hold on
scatter3(grid.p(:,1),grid.p(:,2),g_anony_fun(grid.p(:,1),grid.p(:,2)))
title('Analitic Solution and obstacle')

%% Matrixs and Vectors assembly

%[A , b] = Laplacian_Matrix_assembly_VERSION1(u_bound_anony_fun,grid);
A = Laplacian_fdm_Matrix_assembly(grid);

f=sparse(length(grid.p),1);
f(:) = f_anony_fun(grid.p(:,1),grid.p(:,2));

g=sparse(length(grid.p),1); 
g(:) = g_anony_fun(grid.p(:,1),grid.p(:,2));

% plot matrix disposition
figure_increment = figure_increment+1;
figure(figure_increment)
%spy(A(grid.idx_p_ip,grid.idx_p_ip))
spy(A)
title('Matrix Spy')

% plot matrix disposition
figure_increment = figure_increment+1;
figure(figure_increment)
spy(A(grid.idx_p_ip,grid.idx_p_ip))
%spy(A)
title('Matrix Spy')

%% Algorythm

Np = length(grid.p); % total number of points
Nip = length(grid.idx_p_ip); % number of interior points
u = sparse(Np,1);
u(grid.idx_p_bp) = u_bound_anony_fun(grid.p(grid.idx_p_bp,1),grid.p(grid.idx_p_bp,2));
a = A*u; % vetor das extremidades da matrix, relativamente as condições do bordo

% Au + a + f <= 0 

% Para escrever na forma: Ax-b<=0
AA = A(grid.idx_p_ip,grid.idx_p_ip);
bb = -f(grid.idx_p_ip)+ a(grid.idx_p_ip);
gg = g(grid.idx_p_ip);



Dv=diag(AA); % vetor da diagonal
D = diag(diag(AA));   % D = diag(diag(A)); % matrix diagonal
L = tril(AA)- D;      % L =-tril(A,-1);
U = triu(AA)- D;      % U = -triu(A,1);

II=speye(length(AA)); %matrix identidade

AAminDv = AA - D; % tirar a diagonal



%u(grid.idx_p_ip) = max(-[A(grid.idx_p_ip,grid.idx_p_ip)]\(f(grid.idx_p_ip)+ a(grid.idx_p_ip)) , g(grid.idx_p_ip));
%u(grid.idx_p_ip) = [-A(grid.idx_p_ip,grid.idx_p_ip)]\(f(grid.idx_p_ip)+ a(grid.idx_p_ip)) ;
%u(grid.idx_p_ip) = PSOR([-A(grid.idx_p_ip,grid.idx_p_ip)],(f(grid.idx_p_ip)+ a(grid.idx_p_ip)),g(grid.idx_p_ip));
 
%uk = max(AA\bb,gg); % initial approach (NAO FUNCIONA COMO SULUCAO)
%uk = max(-D\((L+U)*uk-bb),gg);
uk = sparse(rand(length(AA),1));


G_jac = -D\(L+U);

% spectral radius condition
rho_jac = max(abs(eigs(G_jac)))
if rho_jac >= 1
    %error('no convergence')
    fprintf('PSOR no converge')
end
w=0.009
w = 2/(1+sqrt(1-rho_jac^2))
%w=1.3


eps=1e-64
rho=1/eps;



for iter = 1:100000
ukm1 = uk;

% "Projected" Jacobi
%uk = max(-(AAminDv*ukm1-bb)./Dv,gg);
%uk = max(-D\((L+U)*ukm1-bb),gg);

% "Projected" Gaus-Seidel
%uk = max(-(L+D)\(U*ukm1-bb),gg);

% Projected SOR (PSOR)
%uk = max(   (II+w*D\L)\[((1-w)*II-w*D\U)*ukm1 + w*D\bb] ,gg);
%uk = max(   (D+w*L)\[((1-w)*D-w*U)*ukm1 + w*bb] ,gg);
%uk = (D+w*L)\[((1-w)*D-w*U)*ukm1 + w*bb];

% Semi-Smooth Newtons Method
[phi,idx] = min([-AA*ukm1+bb,ukm1-gg],[],2);
da = 2-idx;
Da = sparse(diag(da));
db = idx-1;
Db = sparse(diag(db));

uk = ukm1 - [-Da*AA + Db]\sparse(phi);

% Penalty Method
[~,idx] = max([gg-ukm1,zeros(length(AA),1)],[],2);
P = sparse(rho*diag(2-idx));
%uk = ukm1 - [-AA+P]\(-AA*ukm1+bb-P*(gg-ukm1));

u(grid.idx_p_ip) = uk;

figure(figure_increment)
scatter3(grid.p(:,1),grid.p(:,2),full(u))
zlim([-0.6 1.1])
title('Numeric Solution each iteration')

if norm(uk-ukm1,Inf)< 1e-9
    iter
    break
end 
end

u(grid.idx_p_ip) = uk;



%% Plot Solution

figure_increment = figure_increment+1;
figure(figure_increment)
scatter3(grid.p(:,1),grid.p(:,2),full(u))
hold on
scatter3(grid.p(:,1),grid.p(:,2),g_anony_fun(grid.p(:,1),grid.p(:,2)))
hold on
aa = 0.8294;
teta=-pi:0.01:pi;
x=aa*cos(teta);
y=aa*sin(teta);
plot3(x,y,zeros(1,numel(x)),'k','linewidth',2)
title('Final Numeric Solution and Obstacle')



normlmax = norm(u_ana_anony_fun(grid.p(:,1),grid.p(:,2))-u,Inf)