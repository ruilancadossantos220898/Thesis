% Ax-b <= 0
% g-x <= 0
% (Ax-b)(g-x) = 0

% max{Ax-b , g-x}=0
% min{-Ax+b , x-g}=0

function [x] = SSNM_iterative_min(A,b,g)

%A = L + D + U;
Dv=diag(A); % vetor da diagonal
D = diag(diag(A));   % D = diag(diag(A)); % matrix diagonal
L = tril(A)- D;      % L =-tril(A,-1);
U = triu(A)- D;      % U = -triu(A,1);
II=speye(length(A)); %matrix identidade
AAminDv = A - D; % tirar a diagonal

% initial approach
%uk = max(A\b,g); % (NAO FUNCIONA COMO SULUCAO)
%uk = max(-D\((L+U)*uk-b),g);
uk = sparse(rand(length(A),1)); % algo aleatorio

for iter = 1:100000
ukm1 = uk;

% Semi-Smooth Newtons Method
[phi,idx] = min([-A*ukm1+b,ukm1-g],[],2);
da = 2-idx;
Da = sparse(diag(da));
db = idx-1;
Db = sparse(diag(db));

uk = ukm1 - [-Da*A + Db]\sparse(phi);

if norm(uk-ukm1,Inf)< 1e-12
    %iter
    break
end 

end 

x = uk;


end

