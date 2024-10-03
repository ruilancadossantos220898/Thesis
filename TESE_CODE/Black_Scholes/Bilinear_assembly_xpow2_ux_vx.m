function [B] = Bilinear_assembly_xpow2_ux_vx(x)
Mx = length(x)-2;
B = zeros(Mx+2,Mx+2);
B(1,1) = (x(2)^2+x(1)*x(2)+x(1)^2)/(3*(x(2)-x(1)));

for m = 1:Mx+2

    if m>1 && m < Mx+2 
        B(m,m) = x(m)^2*(1/(x(m)-x(m-1))+1/(x(m+1)-x(m))) +(x(m+1)-x(m-1))/3; 
    end

    if m < Mx+2
        B(m,m+1) = -(x(m)^2/(x(m+1)-x(m)) + x(m) + (x(m+1)-x(m))/3 );
        
    end 

    if m > 1
        B(m,m-1) = -(x(m)^2/(x(m)-x(m-1))- x(m) + (x(m)-x(m-1))/3 );
    end 
    
end

B(end,end) = (x(end)^2+x(end-1)*x(end)+x(end-1)^2)/(3*(x(end)-x(end-1)));

B = sparse(B);

end

