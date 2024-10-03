function [V] = initial_BS_CallOption(r,sigma,d,T,E,S,t)

Ms = length(S)-2;
Nt = length(t)-2;

%% Initializing the matrix of the option values
V= zeros(Ms+2,Nt+2); % V = zeros(Msi+2,Nti+2);

%% Final conditions prescribed by the Call payoff at expiry: V(S,t=T)=max(E-S,0)
V(:,end)=max(S-E, zeros(1,Ms+2) )';
%% Space Boundary conditions prescribed by Call Options:
% as S=0, we have: V(S=0,t)=0
V(1,:)=zeros(1,Nt+2); 

% as S -> infininty (S = Smax), we have: 
% 1) V(S,t)=Se^(-d*(T-t))-Ee^(-r(T-t)) with dividends; 
% 2) V(S,t)=S-Ee^(-r(T-t)) without dividends, (d=0)
% 3) ou simplesmente V(S,t) = S

Smax =S(end);

V(end,:)=Smax*exp(-d*t)-E*exp(-r*(T-t));
%V(end,:)=Smax;

V = sparse(V);

end

