% N é o numero de pontos do interior

function [mesh] = mesh_line1D(xi,xf,N)

% uniform grid
mesh.p = xi:(xf-xi)/(N+1):xf;

% se tem N+1 pontos, então têm N elementos finitos

% A cada elemento finito, temos a ligacao dos 2 nós: (n1->n2) :
mesh.ele = [];

for i=1:length(mesh.p)-1
    mesh.ele = [mesh.ele; i i+1];
end 

%% 'b' (boundary) 
mesh.b = [1; length(mesh.p)];

%% 'in' 
% indices da matriz mesh.p que são no interior do conjunto:
mesh.in=setdiff(1:N+2,mesh.b);

end