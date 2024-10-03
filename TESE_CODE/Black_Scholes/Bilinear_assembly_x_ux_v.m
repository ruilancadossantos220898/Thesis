function [B] = Bilinear_assembly_x_ux_v(x)
    Mx = length(x)-2;

B = zeros(Mx+2,Mx+2);


B(1,1) = -(x(1)/3 + x(2)/6);


for m = 1:Mx+2
    if m>1 && m < Mx+2 
        B(m,m) = -(x(m+1) - x(m-1))/6;
        
    end
    if m < Mx+2
        B(m,m+1) = x(m+1)/6 + x(m)/3;
    end 

    if m > 1
        B(m,m-1) = -(x(m-1)/6 + x(m)/3);
    end 
    
end

B(end,end) = (x(end)/3 + x(end-1)/6);

B = sparse(B);
end

