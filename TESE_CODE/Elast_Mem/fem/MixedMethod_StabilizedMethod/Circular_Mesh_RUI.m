%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descrição da rotina:
% cria uma mesh dentro de uma circunferência
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% -> 1º ARGUMENTO:
%   Irá receber um numero inteiro 'N', corresponderá ao número de
%   discretizações, em intervalos iguais, do perimetro da circunferência.
% -> 2º ARGUMENTO:
%   Irá receber um real 'R', que é o raio da circuferencia.
% -> 3º ARGUMENTO:
%   Irá receber um vetor real 'Xc', com uma componente no eixo 'x' e outra
%   em 'y'. Será o ponto onde a circuferência estará centrada.
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
%   edges   = Uma matrix de todas as "edges" (ou arestas/bordas) existentes
%             na mesh. Cada column corresponde a uma 'edge', possuindo 
%             2 nós: [n1; n2], e cada nó n_i [i=1,2 linha] é um indice da 
%             p-matrix. 
%             Os 2 nós de cada coluna, estão na seguinte ordem: n1<n2.
%
%   t2e     = Uma matrix que conecta 'triangles' and 'edges'. 
%             Cada coluna corresponde a um triangle tendo portanto 'Nt' 
%             colunas, com 3 entradas (linhas) cada 't2e(:,j)=[e1;e2;e3]', que são indices para
%             a matrix 'edges', que remetem para as 'triangle's edges' na ordem: n1->n2, n2->n3, n1->n3. 
%             edges(t2e(1,j)) = [n1 n2], edges(t2e(2,j)) = [n2 n3],edges(t2e(3,j)) = [n1 n3], j=1:Nt
%             onde cada coluna j=1:Nt corresponde por essa mesma ordem
%             à coluna j da matriz 'mesh.t' e portanto a um triangulo (j-ésimo)
%             
%   e2t     = Inverse of t2e. 
%             Cada coluna corresponde a uma aresta/edge da mesh, possuindo o 
%             mesmo numero de colunas que a matriz 'edges', e cada coluna
%             terá 2 entradas (linhas) 'e2t(:,j)=[t1;t2]', 
%             correspondendo aos indices dos dois triangulos que por essa
%             aresta/edge é partilhada. Para j=1:Ne onde cada coluna j=1:Nt 
%             corresponde por essa mesma ordem à coluna j da matriz 'mesh.edge' 
%             e portanto a uma aresta (j-ésima aresta).
%             Para "boundary edges" (arestas da fronteira) a 2ª linha é zero.
%
%   b       = boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mesh] = Circular_Mesh_RUI(N,R,Xc)

%% p t
%{
% Para defenir a fronteira (boundary), que é uma circunferência de raio 'R'
% iremos discretizar o seu perimetro em 'N' intervalos iguais.

% Portanto divida-se o intervalo angular [0,2pi](rads) em 'N' bocados
% diferentes, e coloque-se dentro do vetor 'ang_vec':
ang_vec = linspace(0,2*pi-2*pi/N,N-1);
% Agora vamos defennir um polinomio que simulará o nosso dominio circular,
% com os seus vertices, coincidindo com a circunferencia em questão:
pb_x = R*cos(ang_vec) + Xc(1);
pb_y = R*sin(ang_vec) + Xc(2);

pb = [pb_x', pb_y']; % Vertices do polinomio, que simula a circunferencia.


% De modo a gerar os pontos do dominio, dentro da circunferência, iremos
% dividir o raio da circunferência em 'M' (numero de divisoes do raio) 
% bocados, tal que os comprimentos das divisoes tenham
% mais ou menos a mesma proporcionalidade, das divisoes do perimetro maximo 
% ou seja 2*pi*R/N approx R/M:
M = round(R * 1/(2*pi*R/N)); % 2*pi*R/N = R/M
% Tudo isto para que os triangulos tenham todos mais ou menos a mesma
% proporcionalidade

p = [];
for i=0:M-1
    % com i=0, temos r=R, o raio original, e o raio 'r', vai diminuindo a
    % medida que 'i' sobe:
    r = R * (M-i)/M; % raio
    % numero de divisoes do perimetro do perimetro em consideração:
    n = floor(2*pi*r * 1/(r/(M-i)));

    ang_vec = linspace(0,2*pi -2*pi/n,n);
    p = [p ; [r*cos(ang_vec)' + Xc(1) , r*sin(ang_vec)' + Xc(2)] ];

end
mesh.p = p;

%Lista de triangulos nós(nodes) dos triangles como indices da p-matrix:            
t = delaunay(mesh.p(:,1)',mesh.p(:,2)'); 
mesh.t = sort(t,2); %ordena em todas as linhas, coluna a coluna

Np = max(size(mesh.p)); % Numero total de pontos da mesh
Nt = max(size(mesh.t)); % Numero total de triangulos

%}

%% Outra maneira para 'p' e 't'

gd=[1 0 0 2]';
G=decsg(gd);

[p,ee,t]=initmesh(G,'Hmax',R/N); %R/N = 2/5 = 0.4
%[p,ee,t]=initmesh(G,'Hmax',0.4); %R/N = 2/5 = 0.4
t=sort(t(1:3,:));
mesh.p=p';
mesh.t=t';
Np = max(size(mesh.p)); % Numero total de pontos da mesh
Nt = max(size(mesh.t)); % Numero total de triangulos

%% 't2e' 'edge'

% A cada triangulo, ligacao dos 3 nós: (n1->n2) (n2->n3) (n1->n3):
ligacao=[1 2; 
         2 3; 
         1 3];

edges=[];
for i=1:3
  edges=[edges ; sort([mesh.t(:,ligacao(i,1)),mesh.t(:,ligacao(i,2))],2)]; % empilha as ligacoes ordenadas (n1<n2)
end
% Neste ponto do algoritmo, 'edges' terá 3xNt linhas e 2 colunas.
% Muitas das suas linhas podem estar repetidas, porque afinal de contas, 
% existem triangulos que partilham a mesma edge.
% Para resolver esse assunto sacamos da matrix 'edges' as linhas que são
% repetidas, com o seguinte comando: 
% Exemplo: A = [9; 2; 9; 5]; [C, ia, ic] = unique(A), tal que temos C(ic,:)=A e A(ia,:)=C

[mesh.edges,~, mesh.t2e]= unique(edges,'rows');

% Repare-se que como "mesh.edges(mesh.t2e,:) = edges", ou seja 'mesh.t2e'
% neste momento são os indices de 'mesh.edges' que fazem retornar a 'edges'
% original, e portanto o comprimento de 'edges' e 'mesh.t2e' é igual, com
% um total de 3xNt linhas (a matriz 'mesh.edges' tem menos).

% Devido ao vetor 'ligacao=[1 2;2 3;1 3]', o vetor 'mesh.t2e' terá o
% seguinte significado:

% -> O 1º conjunto de 1:Nt linhas (em que cada linha i corresponde por essa mesma ordem
% à linha i da matriz 'mesh.t' e portanto a um triangulo) são os indices da 
% matriz 'mesh.edges', responsavel pela ligação 
% (n1->n2) = mesh.edges(mesh.t2e(1:Nt),1)->mesh.edges(mesh.t2e(1:Nt),2)

% -> O 2º conjunto de i=1+Nt:2*Nt linhas (em que cada linha i corresponde por essa mesma ordem
% à linha i-Nt da matriz 'mesh.t' e portanto a um triangulo) são os indices da 
% matriz 'mesh.edges', responsavel pela ligação 
% (n2->n3) = mesh.edges(mesh.t2e(1+Nt:2*Nt),1)->mesh.edges(mesh.t2e(1+Nt:2*Nt),2)

% -> O 3º conjunto de i=1+2*Nt:3*Nt linhas (em que cada linha i corresponde por essa mesma ordem
% à linha i-2*Nt da matriz 'mesh.t' e portanto a um triangulo) são os indices da 
% matriz 'mesh.edges', responsavel pela ligação 
% (n1->n3) = mesh.edges(mesh.t2e(1+2*Nt:3*Nt),1)->mesh.edges(mesh.t2e(1+2*Nt:3*Nt),2)

% agora empilhe-se estes 3 conjuntos, em 3 colunas diferentes, com o
% reshape:
mesh.t2e=reshape(mesh.t2e,Nt,3);
%ficando com:
% (n1->n2) = mesh.edges(mesh.t2e(:,1),1)->mesh.edges(mesh.t2e(:,1),2) = mesh.t(:,1) -> mesh.t(:,2)
% (n2->n3) = mesh.edges(mesh.t2e(:,2),1)->mesh.edges(mesh.t2e(:,2),2) = mesh.t(:,2) -> mesh.t(:,3)
% (n1->n3) = mesh.edges(mesh.t2e(:,3),1)->mesh.edges(mesh.t2e(:,3),2) = mesh.t(:,1) -> mesh.t(:,3)

% e portanto é verdade que: 
% mesh.edges(mesh.t2e(:,1),1) = mesh.edges(mesh.t2e(:,3),1) =  mesh.t(:,1)
% mesh.edges(mesh.t2e(:,1),2) = mesh.edges(mesh.t2e(:,2),1) = mesh.t(:,2)
% mesh.edges(mesh.t2e(:,2),2) = mesh.edges(mesh.t2e(:,3),2) = mesh.t(:,3)

% O mambo é confuso mesmo!! (esta matrix só serve para calcular a proxima)

%% 'e2t'
%vamos manipular indices de edges/arestas:
t2eee = [mesh.t2e(:,1); mesh.t2e(:,2); mesh.t2e(:,3)]; % empilha a matriz 'mesh.t2e' numa unica coluna
indx_t = 1:Nt; % têm o formato de linha
indx_ttt = repmat(indx_t',3,1); % empilha o vetor indx_t',replicado 3 vezes, numa unica coluna

% Repare-se agora que 't2eee' e 'indx_ttt', têm a mesma dimensão 3xNt
% linhas e 1 coluna. Para i=1:3xNt, temos que t2eee(i) (indice da mesh.edges) é uma aresta do
% triangulo indx_ttt(i) (indice da mesh.t).


% O comando "[C,ia,ic] = unique(A,occurrence)", onde 'occurrence' especifica 
% quais índices retornar caso de valores repetidos. a ocorrência pode ser 
% occurrence='first' (default) ou occurrence='last'.

[e_first,I_first,~]=unique(t2eee,'first'); % e_first = t2eee(I_first)
[e_last,I_last,~]=unique(t2eee,'last'); % e_last = t2eee(I_last)

% Guarda essa informação:
mesh.e2t(e_first,1)=indx_ttt(I_first);
mesh.e2t(e_last,2)=indx_ttt(I_last);

% Um aresta/edge, aparece no maximo 2 vezes 't2eee', e apenas uma vez caso
% seja uma fronteira. Caso só apareca uma vez (seja fronteira) então ai
% teremos uma linha com entradas iguais, meta-se a 2a entrada a zero nesse
% caso:
mesh.e2t((mesh.e2t(:,1)-mesh.e2t(:,2))==0,2)=0;

%% 'b'
% indices da matriz mesh.p que são a fronteira/boundary do conjunto

b = mesh.edges(mesh.e2t(:,2)==0,:);
mesh.b = unique(b(:));


%% outra maneira para 'b' 
% não funciona tão bem
%mesh.b = boundary(mesh.p(:,1),mesh.p(:,2));

%% 'in'
% indices da matriz mesh.p que são no interior do conjunto:
mesh.in=setdiff(1:Np,mesh.b);


%% Meter nas dimensões certas:
if find(size(mesh.p) == max(size(mesh.p)))==1
    mesh.p = mesh.p';
    mesh.t = mesh.t';
    mesh.edges = mesh.edges';
    mesh.t2e = mesh.t2e';
    mesh.e2t = mesh.e2t';
    mesh.b = mesh.b';
    %mesh.in = mesh.in';
end
%% parametro h
% affine mappings to all triangles
A=cell(2,2); %dimensao do numero de triangulos
A{1,1}=mesh.p(1,mesh.t(2,indx_t))'-mesh.p(1,mesh.t(1,indx_t))';
A{1,2}=mesh.p(1,mesh.t(3,indx_t))'-mesh.p(1,mesh.t(1,indx_t))';
A{2,1}=mesh.p(2,mesh.t(2,indx_t))'-mesh.p(2,mesh.t(1,indx_t))';
A{2,2}=mesh.p(2,mesh.t(3,indx_t))'-mesh.p(2,mesh.t(1,indx_t))';

% determinants of all affine mappings
detA=A{1,1}.*A{2,2}-A{1,2}.*A{2,1};

% mesh parameter 'h' for each element
mesh.h=sqrt(abs(detA))';

mesh.hmax = max(mesh.h);


end

