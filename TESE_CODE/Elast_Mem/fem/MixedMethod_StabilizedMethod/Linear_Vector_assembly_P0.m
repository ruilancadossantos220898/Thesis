function [f] = Linear_Vector_assembly_P0(lin_fun,mesh)
Np = max(size(mesh.p));
Nt = max(size(mesh.t));

% initialize an empty matrix
f=sparse(Nt,1);

%
tind = 1:Nt;

[Xq,Wq] = TrianQuad_Gauss_Legendre();

Psi_1_fun = @(x,y) 1-x-y;
Psi_2_fun = @(x,y) x;
Psi_3_fun = @(x,y) y;

% P1 basis functions at quadrature points
Psi_cell{1}=Psi_1_fun(Xq(1,:),Xq(2,:));
Psi_cell{2}=Psi_2_fun(Xq(1,:),Xq(2,:));
Psi_cell{3}=Psi_3_fun(Xq(1,:),Xq(2,:));

% P1 basis function gradients at quadrature points
gradhat_Psi_cell{1}=repmat([-1;-1],1,length(Wq));
gradhat_Psi_cell{2}=repmat([1;0],1,length(Wq));
gradhat_Psi_cell{3}=repmat([0;1],1,length(Wq));

% affine mappings to all triangles
A=cell(2,2); %dimensao do numero de triangulos
A{1,1}=mesh.p(1,mesh.t(2,tind))'-mesh.p(1,mesh.t(1,tind))';
A{1,2}=mesh.p(1,mesh.t(3,tind))'-mesh.p(1,mesh.t(1,tind))';
A{2,1}=mesh.p(2,mesh.t(2,tind))'-mesh.p(2,mesh.t(1,tind))';
A{2,2}=mesh.p(2,mesh.t(3,tind))'-mesh.p(2,mesh.t(1,tind))';

b=cell(2,1);
b{1,1}=mesh.p(1,mesh.t(1,tind))';
b{2,1}=mesh.p(2,mesh.t(1,tind))';

% determinants of all affine mappings
detA=A{1,1}.*A{2,2}-A{1,2}.*A{2,1};

% all inverse transposes of affine mappings
invAt=cell(2,2);
invAt{1,1}=A{2,2}./detA;
invAt{2,1}=-A{1,2}./detA;
invAt{1,2}=-A{2,1}./detA;
invAt{2,2}=A{1,1}./detA;

% map quadrature points to global coordinates
% X = [A]X_hat + (b)
x=bsxfun(@plus,A{1,1}*Xq(1,:)+A{1,2}*Xq(2,:),b{1,1});
y=bsxfun(@plus,A{2,1}*Xq(1,:)+A{2,2}*Xq(2,:),b{2,1});


% mesh parameter 'h' for each element
h=repmat(sqrt(abs(detA)),1,size(Wq,2));

% loop over basis functions and integrate lin
v=repmat(ones(1,length(Wq)),Nt,1);

f=f+sparse(1:Nt,ones(1,Nt),(lin_fun(v,x,y,h)*Wq').*abs(detA),Nt,1);
end

