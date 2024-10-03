
% uniform orthogonal grid 

function [grid] = Grid2D_Circular(N,R,Xc)

%% Rectangular Universal Grid
% Limits
xmin=Xc(1)-R;
xmax=Xc(1)+R;
ymin=Xc(2)-R;
ymax=Xc(2)+R;

% axis lengths
Lx=xmax-xmin
Ly=ymax-ymin;

% uniform grid
Nx = N;
Ny = Nx;

% space step
hx = Lx/(Nx+1)
hy = Ly/(Ny+1);

u_grid=[];
for i=0:Ny+1
    for j=0:Nx+1
        u_grid=[u_grid ; xmin+j*hx , ymin+i*hy];
        %u_grid=[u_grid ; xmin+j*(xmax-xmin)/(Nx-1) , ymin+i*(ymax-ymin)/(Ny-1)];
    end
end

%% Boundary Points

grid.bp=[];
for j=0:Ny+1
    %linha: y = ymin+j*hy
    % quantos pontos da boundary, e em que coordenadas tocam nessa linha?
    py = ymin+j*hy;
    if sqrt((py-Xc(2))^2)<=R
        grid.bp=[grid.bp ; sqrt(R^2-(py-Xc(2))^2) py];
        grid.bp=[grid.bp ; -sqrt(R^2-(py-Xc(2))^2) py];
    end           
end

for i=0:Nx+1
    %linha: x = xmin+i*hx
    % quantos pontos da boundary, e em que coordenadas tocam nessa linha?
    px = xmin+i*hx;
    if sqrt((px-Xc(1))^2)<=R
        grid.bp=[grid.bp ; px sqrt(R^2-(px-Xc(1))^2)];
        grid.bp=[grid.bp ; px -sqrt(R^2-(px-Xc(1))^2)];
    end           
end

[grid.bp,~,~] = unique(grid.bp,'rows');

% Neste ponto do codigo, temos novos pontos

u_grid = [u_grid ; grid.bp];


%% Interior Points
idxin=find( sqrt((u_grid(:,1)-Xc(1)).^2 + (u_grid(:,2)-Xc(2)).^2)<R-0.01*max(hx,hy)  );
grid.ip = u_grid(idxin,:);

if length(u_grid) == length(idxin)
    8888
end 

%% All points
grid.p = [];
grid.p = [grid.ip ; grid.bp];

%% Indexs Interior and Bound

%grid.idx_p_ip=find(sqrt( sqrt((grid.p(:,1)-Xc(1)).^2 + (grid.p(:,2)-Xc(2)).^2))<R-0.1*max(hx,hy)  );
%grid.idx_p_bp = setdiff(1:length(grid.p),grid.idx_p_ip);
grid.idx_p_ip=1:length(grid.ip);
grid.idx_p_bp=length(grid.ip)+(1:length(grid.bp));

grid.hy = hy;
grid.hx = hx;



end

