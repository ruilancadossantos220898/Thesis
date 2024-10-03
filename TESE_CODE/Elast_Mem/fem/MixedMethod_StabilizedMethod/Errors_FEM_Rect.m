
function [e_H1,e_L2,e_inf] = Errors_FEM_Rect(uh,mesh,u,ux,uy)

Np = max(size(mesh.p));
Nt = max(size(mesh.t));

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

% mesh parameter 'h' for each element
h=repmat(sqrt(abs(detA)),1,size(Wq,2));

% map quadrature points to global coordinates
x=bsxfun(@plus,A{1,1}*Xq(1,:)+A{1,2}*Xq(2,:),b{1,1});
y=bsxfun(@plus,A{2,1}*Xq(1,:)+A{2,2}*Xq(2,:),b{2,1});

% get exact solution values at quadrature points
uval=u(x,y);
% get exact solution x derivative values at quadrature points
uxval=ux(x,y);
% get exact solution y derivative values at quadrature points
uyval=uy(x,y);


% uh evaluated at quadrature points
uhval=zeros(Nt,length(Wq));
% x derivative of uh evaluated at quadrature points
uhxval=zeros(Nt,length(Wq));
% y derivative of uh evaluated at quadrature points
uhyval=zeros(Nt,length(Wq));


for i=1:3
    % basis function values at quadrature points
    %phival=repmat(phi{i},nt,1);
    Psival=repmat(Psi_cell{i},Nt,1);
    % basis function x derivative values at quadrature points
    %phixval=bsxfun(@times,invAt{1,1},gradhat_phi{i}(1,:))+bsxfun(@times,invAt{1,2},gradhat_phi{i}(2,:));
    Psixval=bsxfun(@times,invAt{1,1},gradhat_Psi_cell{i}(1,:))+bsxfun(@times,invAt{1,2},gradhat_Psi_cell{i}(2,:));

    % basis function y derivative values at quadrature points
    %phiyval=bsxfun(@times,invAt{2,1},gradhat_phi{i}(1,:))+bsxfun(@times,invAt{2,2},gradhat_phi{i}(2,:));
    Psiyval=bsxfun(@times,invAt{2,1},gradhat_Psi_cell{i}(1,:))+bsxfun(@times,invAt{2,2},gradhat_Psi_cell{i}(2,:));


    % add contribution to uh
    uhval=uhval+bsxfun(@times,uh(mesh.t(i,:)),Psival);
    % add contribution to uhx
    uhxval=uhxval+bsxfun(@times,uh(mesh.t(i,:)),Psixval);
    % add contribution to uhy
    uhyval=uhyval+bsxfun(@times,uh(mesh.t(i,:)),Psiyval);
end

% compute the L2 norm in each element
L2K=sqrt( (uhval-uval).^2*Wq'.*abs(detA) );

% compute the H1 seminorm in each element
H1K=sqrt( ((uhxval-uxval).^2+(uhyval-uyval).^2)*Wq'.*abs(detA) );

% global L2 norm
e_L2=sqrt( sum(L2K.^2) );

% global H1 seminorm
e_H1=sqrt( sum(H1K.^2) );

% global inf norm
e_inf=norm(uh([mesh.b , mesh.in])-u(mesh.p(1,[mesh.b , mesh.in])',mesh.p(2,[mesh.b , mesh.in])'),inf);
 