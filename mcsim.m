% main program for Monte-Carlo simulation of the XY model with the
% Metropolis algorithm

global D L h theta beta

%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_therm = 1000;  % thermalization sweeps
N_prod  = 20000; % production sweeps

beta = 0.6;
L = 24;

D    = 2;         % dimension

%%%%%%%%%%%%%%%%%%%%%% geometry and initialization %%%%%%%%%%%%%%%%%%%%%%%%

vol = L^D;
h = hop();
theta = rand(vol,1)*2*pi; % hot start


%%%%%%%%%%%%%%%%%%%%%% thermalization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:N_therm
   [~, ~] = sweep();
end
%plot_cnfg(theta);

%%%%%%%%%%%%%%%%%%%%%% production %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize observables
H   = zeros(N_prod,1); % internal energy
M   = zeros(N_prod,1); % abs of magnetization
chi = zeros(N_prod,1); % mag. susceptibility
acc = zeros(N_prod,1); % acceptance rate
DH  = zeros(N_prod,1); % Proposed energy changes
Iam1= zeros(N_prod,1); % <Iam1> = 1
G   = zeros(N_prod,L); % slice-correlation

for k=1:N_prod
   [acc(k), Iam1(k)] = sweep();
   H(k)   = energy(theta);
   M(k)   = magnetization(theta);
   chi(k) = susceptibility(theta);
   G(k,:) = corr(theta);
end

%%%%%%%%%%%%%%%%%%%%%% data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mH, dH, ddH, tauH, dtauH] = UWerr(H);
[mM, dM, ddM, tauM, dtauM] = UWerr(M);
[mchi, dchi, ddchi, tauchi, dtauchi] = UWerr(chi);
macc = sum(acc)/vol/N_prod;
[mIam1, dIam1, ddIam1, tauIam1, dtauIam1] = UWerr(Iam1);
for t=1:L
   [mG(t), dG(t), ddG(t), tauG(t), dtauG(t)] = UWerr(G(:,t));
end

f = @(x) 1 / (log(x(1)/x(2))); % function to extract xi
xidat = [G(:,ceil(L/8)) G(:,ceil(L/8)+1)];
[xi, dxi, ddxi, tauxi, dtauxi] = UWerr(xidat,[],[],'',f);
   
fprintf('beta = %f\n',beta);
fprintf('energy          <H>   = %f +/- %f\n',mH,dH);
fprintf('magnetization   <M>   = %f +/- %f\n',mM,dM);
fprintf('susceptibility  <chi> = %f +/- %f\n',mchi,dchi);
fprintf('   tau_susc           = %f +/- %f\n',tauchi,dtauchi);
fprintf('acceptance rate: %f %%\n',macc*100);
fprintf('Iam1  <exp(-beta*DH)> = %f +/- %f\n',mIam1,dIam1);
fprintf('xi                    = %f +/- %f\n',xi,dxi);
fprintf('G(L/2)                = %f +/- %f\n',mG(L/2+1),dG(L/2+1));
fprintf('   tau_G              = %f +/- %f\n',tauG(L/2+1),dtauG(L/2+1));
figure();
errorbar([0:L-1],mG,dG,'ko');