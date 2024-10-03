%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descrição da rotina:
% cria uma mesh uniforme (ou regular)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% -> 1º ARGUMENTO:
%   Irá receber um real 'ax'
% -> 2º ARGUMENTO:
%   Irá receber um real 'bx'
% -> 3º ARGUMENTO:
%   Irá receber um inteiro 'Nx' que ditará o número de pontos em que o eixo 'x'
%   está discretizado, mais especificamente no intervalo [ax,bx]
% -> 4º ARGUMENTO:
%   Irá receber um real 'ay'
% -> 5º ARGUMENTO:
%   Irá receber um real 'by'
% -> 6º ARGUMENTO:
%   Irá receber um inteiro 'Ny' que ditará o número de pontos em que o eixo 'y'
%   está discretizado, mais especificamente no intervalo [ay,by]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% -> TRIMESH STRUCTURE: 'mesh' (mesh.p, mesh.t, mesh.edges, mesh.t2e, mesh.e2t)
% 
%   p       = 'N' nós (ou "nodes") numa 2xN-matrix, sendo cada coluna 
%              associada a um ponto, 1ª linha a coordenada em x, 2ª linha 
%              associada a coordenada em y.
%
%   t       = Triangles numa 3xNt-matrix (cada coluna associada a um 
%             triangulo, possui os 3 nós do triangulo [n1; n2; n3], e cada 
%             nó n_i [i=1,2,3 linha] é um indice da p-matrix).
%             Os 3 nós de cada coluna, estão na seguinte ordem: n1<n2<n3.
%             'Nt' é o numero total de triangulos.
%
%   b        = boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [mesh_r] = Retangular_Mesh_RUI(ax,bx,Nx,ay,by,Ny)

mesh_r.p = []; % Matriz de Pontos 

% Discretização do eixo "x" e "y", em pontos igualmente distanciados:
x = linspace(ax,bx,Nx);
y = linspace(ay,by,Ny);

% Gerar Todos os pontos do feixo do dominio, e guarda-los na matrix 'Points_Mat'
for i=1:Nx
    for j=1:Ny
        mesh_r.p = [mesh_r.p; x(i), y(j)];
    end 
end 

% este comando gera automaticamente 
% a matriz 'Triangles_Mat', através da "triangulação de Delaunay", que é
% uma forma de triangularizar pontos muito particular

%Lista de triangulos nós(nodes) dos triangles como indices da p-matrix:            
mesh_r.t = delaunay(mesh_r.p(:,1)',mesh_r.p(:,2)');% este comando gera automaticamente 
% a matriz 'Triangles_Mat', através da "triangulação de Delaunay", que é
% uma forma de triangularizar pontos muito particular

mesh_r.t = sort(mesh_r.t,2); %ordena em todas as linhas, coluna a coluna

Nv = max(size(mesh_r.p)); % Numero total de pontos (ou vertices) da mesh 
Nt = max(size(mesh_r.t)); % Numero total de triangulos


% Vertices (ou cantos) do nosso poligno, boundary, 
% Para dar alguma ordem, decidimos Ordenar no sentido 
% contrario dos ponteiros do relogio, a comecar no angulo 0º até 360º,
% portanto o ponto incial têm de ser igual ao ponto final:
pv = [bx,by ; ax,by ; ax,ay ; bx,ay ; bx,by];


[in,on]= inpolygon(mesh_r.p(:,1),mesh_r.p(:,2),pv(:,1),pv(:,2));
% 'in' que estão em todo o fecho do dominio Omega (interior + fronteira)
% 'on' (on the bourder) que estão na fronteira de Omega
% onde o 'in' e o 'on' são arrays de "logicals", ou seja são valores de '0' 
% ou '1',

% Com o intuito de fazer uma distinção de pontos da fronteira e pontos do
% interior do dominio, criou-se:
%Lista de Indices dos pontos interiores do array mesh_r.p:
mesh_r.in = find(in&~on == 1); %Lista de indices
%Lista de Indices dos pontos de fronteira do array mesh_r.p:
mesh_r.b = find(in&~on == 0); %Lista de indices


%% Meter nas dimensões certas:
if find(size(mesh_r.p) == max(size(mesh_r.p)))==1
    mesh_r.p = mesh_r.p';
    mesh_r.t = mesh_r.t';
    mesh_r.b = mesh_r.b';
    mesh_r.in = mesh_r.in';
end


end

