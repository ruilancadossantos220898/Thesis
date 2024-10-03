function [Xq,Wq] = TrianQuad_Gauss_Legendre()

% Numero de nós (raizes dos polinomios Gauss-Legendre):
avv = 5;
if avv ==1
    %Pontos
    ee = [0];
    %Pesos
    w=[2];
end
if avv ==2
    %Pontos
    ee = [-sqrt(1/3) ; sqrt(1/3)];
    %Pesos
    w=[1; 1];
end
if avv ==3
    %Pontos
    ee = [-sqrt(3/5); 0 ;sqrt(3/5)];
    %Pesos
    w=[5/9; 8/9 ;5/9];
end
if avv ==4
    %Pontos
    ee = [-sqrt(3/7+(2/7)*sqrt(6/5)),-sqrt(3/7-(2/7)*sqrt(6/5)),sqrt(3/7-(2/7)*sqrt(6/5)),sqrt(3/7+(2/7)*sqrt(6/5))];
    %Pesos
    w = [(18-sqrt(30))/36, (18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
end
if avv ==5
    %Pontos
    ee = [-(1/3)*sqrt(5+2*sqrt(10/7)),-(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))];
    %Pesos
    w = [(322-13*sqrt(70))/(900),(322+13*sqrt(70))/(900),128/225,(322+13*sqrt(70))/(900),(322-13*sqrt(70))/(900)];
end

%Int_val = 0;
Xq = [];
Wq = [];

%% Transformação para triângulo


for i=1:avv
    for j = 1:avv
        k = i*j;
        x_k = (1-(1+ee(j))/4)*(1+ee(i))/2;
        y_k = (1-(1+ee(i))/4)*(1+ee(j))/2;
        Xq = [Xq, [x_k;y_k]];

        c_k = (1-(2+ee(i)+ee(j))/4)/4*w(i)*w(j);
        Wq = [Wq c_k];
        %G = At*[x_k;y_k] + bt;

        %Int_val = Int_val + c_k * f_anony_func( G(1), G(2) ) ;

    end 
end 
%Int_val = Int_val * abs( det(At) ) ;


end

