function [V] = initial_BS_PutOption_European(r,sigma,d,T,K,S,t)
Ms = length(S)-2;
Nt = length(t)-2;

%% Initializing the matrix of the option values
V= zeros(Ms+2,Nt+2); % V = zeros(Msi+2,Nti+2);

%% Final conditions prescribed by the Put payoff at expiry: V(S,T)=max(S-E,0)
V(:,end)=max(K-S, zeros(1,Ms+2) )';

%% Space Boundary conditions prescribed by Put Options:
% as S=0, we have: V(S=0,t)=K*e^(-r*(T-t))
V(1,:)=K*exp(-r*(T-t)); 
%V(1,:)=max(K*exp(-r*(T-t)) , max(K-S(1), zeros(1,Nt+2) ) ); 

% as S -> infininty (S = Smax), we have: 
V(end,:)=0;
%V(end,:)=max(0 , max(K-S(end),0) );

V = sparse(V);
end

