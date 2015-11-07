function [s,mxy,mz] = SPGR_isochromat_sim(theta,phi_0,TR,T1,T2,Npulse,varargin) 

% [s,mxy,mz] = SPGR_isochromat_sim(theta, phi, TR, T1, T2,Npulse,varargin)
%    ***** Required Arguments ******
%       theta = flip/rad
%       phi   = spoil phase increment/rad
%       T1/T2 = relaxation times/ms
%       TR    = TR /ms
%       Npulse= Number of TR periods to simulate
%
%    ***** Optional Arguments ******
%       To specify these use text argument 'diff', 'psi', or 'Niso' followed by
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
%       
%       psi: set of dephasing angles (radians)
%
%       Niso: number of isochromats, in this case computed as 2*pi*(0:fix(Niso)-1)/fix(Niso)
%   Note, only one of psi/Niso should be set
% Shaihan Malik Oct.2015

%%% Set up dephasing gradients
Niso=2*(Npulse-1)+1;
psifun = @(n)(2*pi*(0:fix(n)-1)/fix(n));
psi = psifun(Niso);

for ii=1:length(varargin)
  
    % Range of dephasing angles
    if strcmpi(varargin{ii},'psi')
        psi = varargin{ii+1};
        Niso= length(psi);
    end
    
    % Gradient params: needs fields diff.A diff.tau diff.D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
    
    % User defined number of isochromats
    if strcmpi(varargin{ii},'Niso')
        Niso = varargin{ii+1};
        % redefine phase range
        psi = psifun(Niso);
    end
    
end

%%% Diffusion: compute b-factors
if exist('diff','var')
    gamma = 42.58e6 * 1e-3 * 2*pi; % sec^-1 mT^-1
    A = diff.G*diff.tau'*1e-3; %mT/m * sec
    tau = sum(diff.tau)*1e-3; % sec
    beff = (gamma*A).^2*tau;
    % combine diffusion effects
    d = beff*diff.D;
    diffusion_calc = true;
else
    diffusion_calc = false;
end

%%% Number of states
N = 3*Niso;

%%% Gradient dephasing matrices
rg={};
for jj=1:Niso
    rg{jj} = rotmat([0 0 psi(jj)]);
end
Rg = blkdiag(rg{:});
    
    
%%% Relaxation matrices 
E1=exp(-TR/T1);
E2=exp(-TR/T2);
E=eye(N);
ii = 1:(N/3);
E(3*N*(ii-1)+3*(ii-1)+1)=E2;
E(3*N*ii-2*N+3*(ii-1)+2)=E2;
E(3*N*ii-N+3*(ii-1)+3)=E1;

%%% composite matrix
Reg = E*Rg;
Reg=sparse(Reg);

%%% Now run the sim
M = zeros([N Npulse]);

% Initialize
M0 = zeros([N 1]);
M0(3:3:end)=1;

%%% Initialize RF rotation matrix, which is modified but not re-declared
%%% each time
T=sparse(zeros([N N]));

%%% Function to get quadratic phase
phi0 = @(p)(((p-1).*(p-2)/2)*phi_0);
build_T_matrix_sub_implicit(rotmat(theta*[cos(phi0(1)) sin(phi0(1)) 0]));
M(:,1) = T*M0;

% Prepare variables for recovery of z-magnetization
xidx=1:3:N;
yidx=2:3:N;
zidx=3:3:N;
zf = (1-E1);

if diffusion_calc
    %%% Build a k-space filter function
    kk = -floor((Niso)/2):floor((Niso-1)/2); 
    %kk = kk-1;
    %filtk1 = exp(-d*kk.^2);
    %filtk2 = exp(-d*(kk.^2+kk+1/3));
    
    %%% 
    gmT = 42.58e6 * 1e-3 * 2*pi; % sec^-1 mT^-1
    tau = diff.tau(:)*1e-3; % convert to sec
    dur = sum(tau); % total duration of dephasing period to consider
    G = diff.G(:)';
    dk = gmT*G*tau; %<--- dk is total dephasing between two EPG states
 
   
    % longitudinal is easy to calculate
    bLong = @(n)((n*dk).^2*dur);
    
    % transverse from the helper function
    bTrans = @(n)(bfactors(n));
    
    filtk1=zeros(size(kk));filtk2=zeros(size(kk));
    kk=abs(kk);
    for ii=1:length(kk)
        filtk1(ii) = exp(-diff.D*bLong(kk(ii)));
        filtk2(ii) = exp(-diff.D*bTrans(kk(ii)));
    end
    
    filtk1=ifftshift(filtk1);
    filtk2=ifftshift(filtk2);
end

% We are looking at F0 (FID) here
% Loop over pulses: have already computed first FID. First dephase then
% flip
for jj=2:Npulse
    
    % Get matrix ready for flip
    build_T_matrix_sub_implicit(rotmat(theta*[cos(phi0(jj)) sin(phi0(jj)) 0]));

    % Dephase and T1 recovery
    M(:,jj) = Reg*M(:,jj-1);
    M(zidx,jj) = M(zidx,jj) + zf;
    
    if diffusion_calc

        % Now apply diffusion filtering
        mxy = M(xidx,jj)+1i*M(yidx,jj);
        mz = M(zidx,jj);
        
        % Apply
        mxyf = ifft(fft(mxy(:)).*filtk2(:));
        mzf = ifft(fft(mz(:)).*filtk1(:));
        
        % Put these back in the array
        M(xidx,jj)=real(mxyf);
        M(yidx,jj)=imag(mxyf);
        M(zidx,jj)=real(mzf);
    end
    
    % Flip
    M(:,jj) = T*M(:,jj);
    
end

% Now generate signal and demodulate it
mxy = M(1:3:end,:) + 1i*M(2:3:end,:);
% Get signal from mean
s = 1i*mean(mxy,1);
% demodulate this
s = s .* exp(-1i*phi0(1:Npulse));

%%% Also return mz if needed
if nargout==3
    mz=M(3:3:end,:);
end

    function build_T_matrix_sub_implicit(AA)
        %%% This function operates on the existing T matrix, rather than
        %%% re-declare it each time. This is much faster 
        ix = 1:(N/3);
        T(3*N*(ix-1)+3*(ix-1)+1)=AA(1);
        T(3*N*(ix-1)+3*(ix-1)+2)=AA(2);
        T(3*N*(ix-1)+3*(ix-1)+3)=AA(3);
        T(3*N*ix-2*N+3*(ix-1)+1)=AA(4);
        T(3*N*ix-2*N+3*(ix-1)+2)=AA(5);
        T(3*N*ix-2*N+3*(ix-1)+3)=AA(6);
        T(3*N*ix-N+3*(ix-1)+1)=AA(7);
        T(3*N*ix-N+3*(ix-1)+2)=AA(8);
        T(3*N*ix-N+3*(ix-1)+3)=AA(9);
    end

    % Rotation matrix function
    function R = rotmat(u)
        
        % check zero input
        if any(u)
            th=norm(u);
            u=u/th;
            ct=cos(th);
            st=sin(th);
            R = [[ct + u(1)^2*(1-ct) u(1)*u(2)*(1-ct)-u(3)*st u(1)*u(3)*(1-ct)+u(2)*st];...
                [u(2)*u(1)*(1-ct)+u(3)*st ct+u(2)^2*(1-ct) u(2)*u(3)*(1-ct)-u(1)*st];
                [u(3)*u(1)*(1-ct)-u(2)*st u(2)*u(3)*(1-ct)+u(1)*st] ct+u(3)^2*(1-ct)];
            
        else
            R=eye(3);
        end
        
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
        for ll=1:length(G)
            kf = ki+gmT*G(ll)*tau(ll);
            bT = bT + (tau(ll)/3)*(ki.^2+kf.^2+ki.*kf);
            
            % now update ki for next iteration
            ki=kf;
        end
    end

end
