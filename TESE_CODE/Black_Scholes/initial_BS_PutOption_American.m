function [V] = initial_BS_CallOption_American(r,sigma,d,T,K,S,t)
Ms = length(S)-2;
Nt = length(t)-2;

%% Initializing the matrix of the option values
V= sparse(Ms+2,Nt+2); % V = zeros(Msi+2,Nti+2);

%% Final conditions prescribed by the Put payoff at expiry: V(S,T)=max(S-E,0)
% Obstacle:
g = max(K-S, zeros(1,Ms+2) );
%V(:,end)=max(E-S, zeros(1,Ms+2) )';
V(:,end)=g;

%% Space Boundary conditions prescribed by Put Options:
% as S=0, we have: V(S=0,t)=K*e^(-r*(T-t))
h_Smin = K*exp(-r*(T-t));
%V(1,:)=K*exp(-r*(T-t)); 
V(1,:)=max(h_Smin,g(1)); 

% as S -> infininty (S = Smax), we have: 
%V(end,:)=0;
h_Smax = 0;
V(end,:)=max(h_Smax , g(end));

end

