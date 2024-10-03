function [B] = Bilinear_assembly_u_v(x)
Mx = length(x)-2;
B = zeros(Mx+2,Mx+2);

%B(1,1) = (x(2)-x(1))/3;

for m = 1:Mx+2
    % (0,0)
    if m==1
        B(1,1) = (x(m+1)-x(m))/3;
    end 
    % (m,m)
    if m>1 && m < Mx+2 
        B(m,m) = (x(m+1)-x(m-1))/3;
    end
    % (m,m+1)
    if m < Mx+2
        B(m,m+1) = (x(m+1)-x(m))/6;
    end 
    % (m,m-1)
    if m > 1
        B(m,m-1) = (x(m)-x(m-1))/6;
    end 
    % (M+1,M+1)
    if m==Mx+2
        B(m,m) = (x(m)-x(m-1))/3;
    end
    
end

%B(end,end) = (x(end)-x(end-1))/3;
B=sparse(B);
end

