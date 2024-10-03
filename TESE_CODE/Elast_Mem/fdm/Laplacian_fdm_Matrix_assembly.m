function [A] = Laplacian_fdm_Matrix_assembly(grid)

% initialize an empty matrix
Np = length(grid.p); % total number of points
Nip = length(grid.idx_p_ip); % number of interior points
A=sparse(Np,Np);

p_ipp_lig = zeros(Nip,5);
h = zeros(Nip,5);
cx = zeros(Nip,3);
cy = zeros(Nip,3);

for i=1:Nip

    %p_ipip_dif = grid.p(grid.idx_p_ip,:) - grid.p(grid.idx_p_ip(i),:);
    %p_ipbp_dif = grid.p(grid.idx_p_bp,:) - grid.p(grid.idx_p_ip(i),:);
    p_ipp_dif = grid.p - grid.p(grid.idx_p_ip(i),:);

    %South
    I = find(p_ipp_dif(:,1) == 0 & p_ipp_dif(:,2)<0);
    [hym1,id] = min(abs(p_ipp_dif(I,2)));
    p_ipp_lig(i,1) = I(id);
    h(i,1) = hym1;

    %West
    I = find(p_ipp_dif(:,2) == 0 & p_ipp_dif(:,1)<0);
    [hxm1,id] = min(abs(p_ipp_dif(I,1)));
    p_ipp_lig(i,2) = I(id);
    h(i,2) = hxm1;

    %P
    p_ipp_lig(i,3) = grid.idx_p_ip(i);

    %E
    I = find(p_ipp_dif(:,2) == 0 & p_ipp_dif(:,1)>0);
    [hx,id] = min(p_ipp_dif(I,1));
    p_ipp_lig(i,4) = I(id);
    h(i,4) = hx;

    %N
    I = find(p_ipp_dif(:,1) == 0 & p_ipp_dif(:,2)>0);
    [hy,id] = min(p_ipp_dif(I,2));
    p_ipp_lig(i,5) = I(id);
    h(i,5) = hy;

    cx(i,:) = [2*hx, -2*(hx+hxm1), 2*hxm1]/(hx*hxm1*(hx+hxm1));
    cy(i,:) = [2*hy, -2*(hy+hym1), 2*hym1]/(hy*hym1*(hy+hym1));
    c(i,:) = [cy(i,1),cx(i,1),cx(i,2)+cy(i,2),cx(i,3),cy(i,3)];

    for d=1:5
        A(grid.idx_p_ip(i),p_ipp_lig(i,d)) = c(i,d);     
    end 

end 


end

