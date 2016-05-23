%% 
function [F0,Fn,Zn,F] = SPGR_EPG_sim(theta, phi, TR, T1, T2,Npulse,varargin)

% [F0,Fn,Zn,F] = SPGR_EPG_sim(theta, phi, TR, T1, T2,Npulse,varargin)
%    ***** Required Arguments ******
%       theta = flip/rad
%       phi   = spoil phase increment/rad
%       T1/T2 = relaxation times/ms
%       TR    = TR /ms
%       Npulse= Number of TR periods to simulate
%
%    ***** Optional Arguments ******
%       To specify these use text argument 'diff' or 'kmax' followed by
%       the structure/value, as next argument.
%
%       diff (structure with following fields):
%            D     = diffusion coeff in m^2 sec^-1 (i.e. order of 10^-9 expected)
%            G     = gradient strength in mT/m
%            tau   = gradient duration in ms 
%            NOTE: The b-value function calculates this based on G and tau.
%            These can either be single values or time resolved, in which
%            case tau(i) is the duration for which gradient G(i) is played.
%            It is possible to include zero gradients here.
%       kmax
%             = Maximum order of configuration state to compute. 
%               Possible values:
%                               - not set: Variable number of states used to
%                                 guarantee F0 will be accurate
%                               - A low number: the calculation might be faster but it could suffer 
%                                 inacccuracies.
%                               - inf: No states are dropped
%
% Function returns predicted signal (F0), Fn, Zn and also whole EPG (F). The 
% layout of F is rows: [F0 Z0 F1 F-1 Z1 F2 F-2 Z2 ... etc]
%                cols: Number of TR periods
%
%
% Shaihan Malik Oct.2015

%% Additional arguments
for ii=1:length(varargin)
  
    % User defined maximum k-value
    if strcmpi(varargin{ii},'kmax')
        kmax_userdef = varargin{ii+1};
    end
    
    % Gradient params: needs fields grad.A grad.tau grad.D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
end


if ~exist('kmax_userdef','var')
    % Set maximum k-value
    kmax = ceil(Npulse/2 + 1);
    
    % For efficiency change the number of states evaluated at each pulse
    % Exclude evaluation of states that won't contribute to any echo amplitudes
    kloopmax = zeros([Npulse 1]);
    kloopmax(1:kmax)=1:kmax;
    tmp=kmax:-1:1;
    kloopmax((kmax+1):end)=tmp(1:(Npulse-kmax));
else
    if isinf(kmax_userdef)
        % this flags not dropping any states
        kloopmax = (1:Npulse)-1;
        kmax=Npulse-1;
    else
        % user has set a maximum value
        kmax_theoretical = ceil(Npulse/2 + 1);
        if kmax_userdef<kmax_theoretical
            kmax = kmax_userdef;
        else
            kmax = kmax_theoretical;
        end
        
        kloopmax = zeros([Npulse 1]);
        kloopmax(1:kmax_theoretical)=1:kmax_theoretical;
        tmp=kmax_theoretical:-1:1;
        kloopmax((kmax_theoretical+1):end)=tmp(1:(Npulse-kmax_theoretical));
        
        kloopmax(kloopmax>kmax)=kmax;
    end
end


% Store everything in one array, [F0 Z0 F1 F-1 Z1 F2 F-2 Z2 ... etc]
F = zeros([3*kmax+2 Npulse]); %<- state after pulse
Fr = zeros([3*kmax+2 Npulse]);%<- state after relaxation

% Parameters
E1=exp(-TR/T1);
E2=exp(-TR/T2);


%%% Diffusion
if exist('diff','var')
    
    %%% We are going to include diffusion, so prepare some values
    gmT = 42.58e6 * 1e-3 * 2*pi; % sec^-1 mT^-1
    tau = diff.tau(:)*1e-3; % convert to sec
    dur = sum(tau); % total duration of dephasing period to consider
    G = diff.G(:)';
    dk = gmT*G*tau; %<--- dk is total dephasing between two EPG states
 
   
    % longitudinal is easy to calculate
    bLong = @(n)((n*dk).^2*dur);
    
    % transverse from the helper function
    bTrans = @(n)(bfactors(n));
    
    D = diff.D;
    
    %%% More efficient to pre-compute these now.
    dvals = zeros([kmax+1 2]);
    for ii=0:kmax
        dvals(ii+1,1)=bLong(ii)*D;
        dvals(ii+1,2)=bTrans(ii)*D;
    end
else
    dvals = zeros([kmax+1 2]);
end


%% simulation

% equilibrium at start
Fr(2,1)=1;

% Phase for later demodulation
phi0 = zeros([Npulse 1]);

for p=1:Npulse    
    p0=((p-1)*(p-2)/2)*phi; % sum of quadratic sequence
    % save this phase
    phi0(p)=p0;
    
    % RF pulse transition matrix
    A = Trot(theta,p0);
    
    % 1. RF pulse: treat zero order state separately
    X = [Fr(1,p);conj(Fr(1,p));Fr(2,p)];
    X_p = A*X;
    F(1:2,p)=X_p([1 3]);% exclude F0*
    
    for ii=1:kloopmax(p) 
        kidx = (1:3)+2+3*(ii-1);
        F(kidx,p) = A*Fr(kidx,p);
    end
    
    %%% If we're on the last pulse, break here
    if p==Npulse
        break
    end
    
    % 2. Relaxation and dephasing
    
    % First apply relevant relaxation & diffusion attenuation 
    % Zero Order
    %bt = bTrans(0);
    %Fr(1,p+1) = F(1,p) * E2 * exp(-bt*D);
    Fr(1,p+1) = F(1,p) * E2 * exp(-dvals(1,2));
    Fr(2,p+1) = F(2,p) * E1 + 1-E1;
    % Now higher orders
    for ii=1:kloopmax(p)
        
        %bt = bTrans(ii);
        %bl = bLong(ii);
        
        %e1 = E1*exp(-bl*D);
        %e2 = E2*exp(-bt*D); 
        e1 = E1*exp(-dvals(ii+1,1));%<-- add 1 to index because 1st value is N=0
        e2 = E2*exp(-dvals(ii+1,2)); 
        fidx = (1:2)+2+3*(ii-1);
        zidx = (3)+2+3*(ii-1);
        Fr(fidx,p+1) = F(fidx,p) * e2;
        Fr(zidx,p+1) = F(zidx,p) * e1;
    end
       
    % Now shift transverse states
    fidx = [fliplr(4:3:(kmax*3+1)) 1 3:3:(kmax*3+1)];
    Fr(fidx,p+1)=circshift(Fr(fidx,p+1),[1 0]);
    
    % conjugate for 0th order
    Fr(1,p+1)=conj(Fr(1,p+1));
end
F0=F(1,:);

%%% phase demodulate
F0 = F0(:) .* exp(-1i*phi0(:)) *1i;

%%% Reconstruct Fn and Zn here
fidx = [fliplr(4:3:(kmax*3+1)) 1 3:3:(kmax*3+1)];
Fn=F(fidx,:);
Fn(1:kmax,:)=conj(Fn(1:kmax,:));% F-k is actually F-k*
Zn = F([2 5:3:(kmax*3+1)],:);



    %%% Helper function to define EPG transition matrix
    % As per Weigel et al JMR 2010 276-285
    function T = Trot(a,p)
        T = zeros([3 3]);
        T(1) = cos(a/2).^2;
        T(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*1i*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));
        T(5) = T(1);
        T(6) = 0.5*1i*exp(1i*p)*sin(a);
        T(7) = -1i*exp(1i*p)*sin(a);
        T(8) = 1i*exp(-1i*p)*sin(a);
        T(9) = cos(a);
    end



    %%% Helper function to compute b-value as function of N (order of EPG state)
    % This is a nested function: it sees arguments from above. See Weigel
    % review paper 2015 (JMRI)
    function bT = bfactors(N)
        
        k0 = N.*dk; %<--- initial k-space is the order of the state * dk
        %%% for each segment we need the initial k-value ki and final value kf
        %%% initialize with ki=k0.
        ki = k0;
        bT=0;
        for jj=1:length(G)
            kf = ki+gmT*G(jj)*tau(jj);
            bT = bT + (tau(jj)/3)*(ki.^2+kf.^2+ki.*kf);
            
            % now update ki for next iteration
            ki=kf;
        end
    end


end