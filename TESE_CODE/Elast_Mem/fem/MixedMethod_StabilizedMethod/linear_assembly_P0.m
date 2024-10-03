function [G] = linear_assembly_P0(fun,mesh)
Np = max(size(mesh.p));
Nt = max(size(mesh.t));
Nv = Np-Nt;

ele_ind = 1:Nt;

% initialize an empty matrix
G=sparse(Nt,1);

% Quadrature points
[Xq,Wq] = TrianQuad_Gauss_Legendre();

% Basis functions
Psi_1_fun = @(x,y) 1-x-y;
Psi_2_fun = @(x,y) x;
Psi_3_fun = @(x,y) y;
Psi_b_fun = @(x,y) 27*(1-x-y).*x.*y;

Psi_b_fun = @(x,y) 27*(1-x-y).*x.*y; % Psi_b_fun = @(x,y) 27*(x.*y-x.^2*y-y.^2*x);
Psi_bx_fun = @(x,y) 27*(y-2*x.*y-y.^2);
Psi_by_fun = @(x,y) 27*(x-2*x.*y-x.^2);

% Basis functions at quadrature points
Psi_cell{1}=Psi_1_fun(Xq(1,:),Xq(2,:));
Psi_cell{2}=Psi_2_fun(Xq(1,:),Xq(2,:));
Psi_cell{3}=Psi_3_fun(Xq(1,:),Xq(2,:));
Psi_cell{4}=Psi_b_fun(Xq(1,:),Xq(2,:));

% P1 basis function gradients at quadrature points
gradhat_Psi_cell{1}=repmat([-1;-1],1,length(Wq));
gradhat_Psi_cell{2}=repmat([1;0],1,length(Wq));
gradhat_Psi_cell{3}=repmat([0;1],1,length(Wq));
gradhat_Psi_cell{4}=[Psi_bx_fun(Xq(1,:),Xq(2,:)) ; Psi_by_fun(Xq(1,:),Xq(2,:))];


%gradhat_Psi_b = [Psi_bx_fun(Xq(1,:),Xq(2,:)) ; Psi_by_fun(Xq(1,:),Xq(2,:))];


% affine mappings to all triangles
A=cell(2,2); %dimensao do numero de triangulos
A{1,1}=mesh.p(1,mesh.t(2,ele_ind))'-mesh.p(1,mesh.t(1,ele_ind))';
A{1,2}=mesh.p(1,mesh.t(3,ele_ind))'-mesh.p(1,mesh.t(1,ele_ind))';
A{2,1}=mesh.p(2,mesh.t(2,ele_ind))'-mesh.p(2,mesh.t(1,ele_ind))';
A{2,2}=mesh.p(2,mesh.t(3,ele_ind))'-mesh.p(2,mesh.t(1,ele_ind))';

b=cell(2,1);
b{1,1}=mesh.p(1,mesh.t(1,ele_ind))';
b{2,1}=mesh.p(2,mesh.t(1,ele_ind))';

% determinants of all affine mappings
detA=A{1,1}.*A{2,2}-A{1,2}.*A{2,1};

% all inverse transposes of affine mappings
invAt=cell(2,2);
invAt{1,1}=A{2,2}./detA;
invAt{2,1}=-A{1,2}./detA;
invAt{1,2}=-A{2,1}./detA;
invAt{2,2}=A{1,1}./detA;

% map quadrature points to global coordinates
x=bsxfun(@plus,A{1,1}*Xq(1,:)+A{1,2}*Xq(2,:),b{1,1});
y=bsxfun(@plus,A{2,1}*Xq(1,:)+A{2,2}*Xq(2,:),b{2,1});


% mesh parameter 'h' for each element
h=repmat(sqrt(abs(detA)),1,size(Wq,2));


lamb=repmat(ones(1,length(Wq)),Nt,1);

for i=1:4

    % function values
    u=repmat(Psi_cell{i},Nt,1);
    % grad_u = invA^T*gradhat_Psi
    ux=invAt{1,1}*gradhat_Psi_cell{i}(1,:)+invAt{1,2}*gradhat_Psi_cell{i}(2,:);
    uy=invAt{2,1}*gradhat_Psi_cell{i}(1,:)+invAt{2,2}*gradhat_Psi_cell{i}(2,:);

    % test function values
    v=repmat(Psi_cell{i},Nt,1);
    vx=invAt{1,1}*gradhat_Psi_cell{i}(1,:)+invAt{1,2}*gradhat_Psi_cell{i}(2,:);
    vy=invAt{2,1}*gradhat_Psi_cell{i}(1,:)+invAt{2,2}*gradhat_Psi_cell{i}(2,:);
    
end

G = G+sparse(1:Nt,ones(1,Nt),(fun(u,v,ux,uy,vx,vy,x,y,lamb,h)*Wq').*abs(detA),Nt,1);





end

