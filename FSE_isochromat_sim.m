function [s,mxy,mz] = FSE_isochromat_sim(theta,Niso,varargin)

% [s,mxy,mz] = FSE_isochromat_sim(theta, varargin)
%    ***** Required Arguments ******
%       theta = flip/rad (set of all flip angles including excitation and all
%               refocusing pulses)
%       Niso = number of isochromats to incluse
%    ***** Optional Arguments ******
%       To specify these use text argument 'diff' or 'kmax' followed by
%       the structure/value, as next argument.
%
%       T1 =    T1/ms
%       T2 =    T2/ms
%       ESP=    echo spacing / ms
%       psi=    isochromat gradient dephasing angles
%
%
% Shaihan Malik Oct.2015
% Note: Diffusion is not implemented in this version of the function.

%%% Arguments
ESP = 7; % echo spacing...
T1 = inf;
T2 = inf;
%%% Set up dephasing gradients
psi = linspace(-0.5,0.5,Niso)*2*pi;

for ii=1:length(varargin)
    
    if strcmpi(varargin{ii},'T1')
        T1 = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T2')
        T2 = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'ESP')
        ESP = varargin{ii+1};
    end
    
    % Range of dephasing angles
    if strcmpi(varargin{ii},'psi')
        psi = varargin{ii+1};
        if length(psi)~=Niso
            disp('error')
            return
        end
    end
end


Necho = length(theta)-1;

%%% CPMG condition
alpha = abs(theta);phi=angle(theta);
phi(2:end) = phi(2:end) + pi/2;



%%% Number of states
N = 3*Niso;

%%% Gradient dephasing matrix
rg={};
for jj=1:Niso
    rg{jj} = rotmat([0 0 psi(jj)]);
end
Rg = blkdiag(rg{:});
    
    
%%% Relaxation matrices - half dephasing
E1=exp(-0.5*ESP/T1);
E2=exp(-0.5*ESP/T2);
E=eye(N);
ii = 1:(N/3);
E(3*N*(ii-1)+3*(ii-1)+1)=E2;
E(3*N*ii-2*N+3*(ii-1)+2)=E2;
E(3*N*ii-N+3*(ii-1)+3)=E1;

%%% composite z-rotate and relaxation matrix
Reg = E*Rg;
Reg=sparse(Reg);

%%% variables for recovery of z-magnetization
Z0 = zeros([N 1]);
Z0(3:3:N) = (1-E1);


%%% Now run the sim
M = zeros([N Necho]);

% Initialize
M0 = zeros([N 1]);
M0(3:3:end)=1;

%%% Initialize RF rotation matrix, which is modified but not re-declared
%%% each time
T=sparse(zeros([N N]));
build_T_matrix_sub_implicit(rotmat(alpha(1)*[cos(phi(1)) sin(phi(1)) 0]));

%%% Excitation pulse, apply RF then dephase/relax
M0 = Reg*T*M0 + Z0;

%%% First refocus pulse
build_T_matrix_sub_implicit(rotmat(alpha(2)*[cos(phi(2)) sin(phi(2)) 0]));
M(:,1) = Reg*T*M0 + Z0;

% Loop over refocus pulses
for ii=2:Necho
    % Update T matrix
    build_T_matrix_sub_implicit(rotmat(alpha(ii+1)*[cos(phi(ii+1)) sin(phi(ii+1)) 0]));

    % Apply Rz/relax, flip, Rz/relax 
    M(:,ii) = Reg*T*(Reg*M(:,ii-1)+Z0)+Z0;
end

%%% Extract value
mxy = M(1:3:N,:) + 1i*M(2:3:N,:);
%%% Get signal from mean
s = 1i*mean(mxy,1);

%%% Also return mz if needed
if nargout==3
    mz=M(3:3:N,:);
end

    function build_T_matrix_sub_implicit(AA)
        %%% This function operates on the existing T matrix, rather than
        %%% re-declare it each time. This is much faster and perhaps could
        %%% also be used for accelerating the EPG simulations
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

end




