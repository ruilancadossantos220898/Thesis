function [Xq,Wq] = Quad_Gauss_Legendre()
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
Xq = ee;
Wq = w;
end

