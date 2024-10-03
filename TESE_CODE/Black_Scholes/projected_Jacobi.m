% Ax-b <= 0
% g-x <= 0
% (Ax+b)(g-x) = 0

% max{Ax+b , g-x}=0

function [x] = projected_Jacobi(A,b,g)

%A = L + D + U;
Dv=diag(A); % vetor da diagonal
D = diag(diag(A));   % D = diag(diag(A)); % matrix diagonal
L = tril(A)- D;      % L =-tril(A,-1);
U = triu(A)- D;      % U = -triu(A,1);
II=speye(length(A)); %matrix identidade
AAminDv = A - D; % tirar a diagonal

end

