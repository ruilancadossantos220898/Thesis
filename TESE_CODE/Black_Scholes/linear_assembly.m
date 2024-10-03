function [F] = linear_assembly(fun,mesh)
Np = max(size(mesh.p));
Ne = max(size(mesh.ele));

% initialize an empty matrix
%B=sparse(Np,Np);
% initialize an empty matrix
F=sparse(Np,1);

[Xq,Wq] = Quad_Gauss_Legendre();

Psi_1_fun = @(x) (-x+1)/2;
Psi_2_fun = @(x) (x+1)/2;

% P1 basis functions at quadrature points
Psi_cell{1}=Psi_1_fun(Xq(:));
Psi_cell{2}=Psi_2_fun(Xq(:));

% P1 basis function gradients at quadrature points
gradhat_Psi_cell{1}=repmat(-1/2,1,length(Wq));
gradhat_Psi_cell{2}=repmat(1/2,1,length(Wq));

% affine mappings to all elements
xi = mesh.p(mesh.ele(:,1))';
xf = mesh.p(mesh.ele(:,2))'; 

A = (xf-xi)/2;
b = (xf+xi)/2;

% mesh parameter 'h' for each element
h=repmat(xf-xi,size(Wq));

% map quadrature points to global coordinates
x = A.*Xq + b;

% determinants of all affine mappings
detA=A;

% all inverse transposes of affine mappings
invAt = 1./A;  


    
    for i=1:2
        % test function values
        v=repmat(Psi_cell{i}',Ne,1);
        vx=invAt.*gradhat_Psi_cell{i};
        

        F=F+sparse(mesh.ele(:,i),ones(size(mesh.ele(:,i))),(fun(v,vx,x,h)*Wq').*abs(detA),Np,1);
        
    end






end