function [F0,Fn,Zn,F] = FSE_EPG_sim(theta,varargin)

% [F0,Fn,Zn,F] = FSE_EPG_sim(theta, varargin)
%    ***** Required Arguments ******
%       theta = flip/rad (set of all flip angles including excitation and all
%               refocusing pulses)
%
%    ***** Optional Arguments ******
%       To specify these use text argument 'diff' or 'kmax' followed by
%       the structure/value, as next argument.
%
%       T1 =    T1/ms
%       T2 =    T2/ms
%       ESP=    echo spacing / ms
%       Npathway = Number of pathways to include
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
% Note: Diffusion is not implemented in this version of the function. See
% the SPGR version for an idea of how to do this


np = length(theta);
kmax = 2*np - 1; % up to just before the next RF pulse
N=3*(2*(np-1)+1); % number of states in total

% starting state
E = zeros([3 1]);
E(3)=1;

% split magnitude and phase
alpha = abs(theta);phi=angle(theta);

% add CPMG phase
phi(2:end) = phi(2:end) + pi/2;

%% get variables
T1 = inf;
T2 = inf;
ESP=7;
Npathway = inf;
klimit=false;


for ii=1:length(varargin)
    
    if strcmp(varargin{ii},'T1')
        T1 = varargin{ii+1};
    end
    if strcmp(varargin{ii},'T2')
        T2 = varargin{ii+1};
    end
    if strcmp(varargin{ii},'ESP')||strcmp(varargin{ii},'TE')
        ESP = varargin{ii+1};
    end
    % # of coherence pathways to consider (default is inf)
    if strcmp(varargin{ii},'Npathway')||strcmp(varargin{ii},'npath')
        Npathway = varargin{ii+1};
    end
    % more drastic version of above (see appendix to Malik et al, doi:10.1002/mrm.24153)
    if strcmp(varargin{ii},'klimit')
        klimit=true;
    end
end

% enforce pathway limit
if (N>Npathway)&&~klimit
    N=Npathway;
end

if klimit
    nr = length(theta)-1; % number of refocus pulses
    kmax = 2*nr;
    
    %%% half this
    kmax=fix(kmax/2);
    if mod(nr,2) 
        % odd
        KMAX = [1:2:kmax (kmax-2):-2:1];
    else
        %even
        KMAX = [1:2:kmax kmax:-2:1];
    end
    NMAX = 3*(KMAX+1);
else
    % previous implementation
    NMAX = 6:6:6*(np-1);
    NMAX(NMAX>N)=(N-mod(N,3));
end
    
    
%% ==== build Shift matrix, S with indices 
S = zeros([N N]);
%%% F(k>1) look @ states just BEFORE np+1 pulse
kidx = 4:3:N; % miss out F1+
sidx = kidx-3;
idx = kidx + N*(sidx-1);
S(idx)=1;

%%% F(k<1) <--- start at F-1 (not F0* which is ignored @#2)
kidx = 5:3:N;
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+3;
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% Z states
kidx = 3:3:N;
ix = kidx + N*(kidx-1);
S(ix)=1;

%%% finally F0+ - relates to F-1
S(1,5)=1;

%%% also need F0*, also relate this to F-1
S(2,5)=1;

%% Relaxation =====
E1=exp(-0.5*ESP/T1);
E2=exp(-0.5*ESP/T2);
R=eye(N);
ii = 1:(N/3);
R(3*N*(ii-1)+3*(ii-1)+1)=E2;
R(3*N*ii-2*N+3*(ii-1)+2)=E2;
R(3*N*ii-N+3*(ii-1)+3)=E1;

%%% regrowth
b = zeros([N 1]);
b(3) = 1-E1;

%%% Composite rotate-shift
RS=R*S;

%% F matrix (many elements zero, not efficient)
F = zeros([N np-1]); %% This will now record state *at each echo* 19-2-15

%% Excitation pulse: Uncommented on 24-1-12
A = Trot(alpha(1),phi(1));
F_ex = A*E; %<---- state straight after excitation [F0 F0* Z0]

% state just after this excitation
F1 = zeros(N,1);
F1(1:3) = F_ex;

%%% relax/shift
F1 = RS*F1+b;

%% First refocusing pulse
A = Trot(alpha(2),phi(2));
T = build_T_matrix_sub(A,6);

%%% apply RF pulse
F(1:6,1) = T*F1(1:6);

%%% relax/shift - this is now the state for the first echo
F(:,1) = RS*F(:,1)+b;
% Deal with complex conjugate after shift
F(1,1)=conj(F(1,1)); %F0 comes from F-1 so conjugate

%% Simulate next refocusing pulses

for jj=2:np-1 
    A = Trot(alpha(jj+1),phi(jj+1));
    
    kidx=1:NMAX(jj);
    
    T = build_T_matrix_sub(A,NMAX(jj));

    % First evolve half ESP 
    FF = RS*F(:,jj-1)+b;
    
    % Deal with complex conjugate after shift
    FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate
    
    % Now flip
    FF(kidx) = T*FF(kidx);
    
    % Now evolve half of next ESP to get the echo
    FF = RS*FF+b;
    FF(1)=conj(FF(1));
    
    % This is now the echo, so store it
    F(:,jj)=FF;
end
F0 = F(1,:)*1i;
%%% Construct Fn and Zn
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)]; 
kvals = -2*(np-1):2*(np-1);

%%% Remove the lowest two negative states since these are never populated
%%% at echo time
idx(1:2)=[];
kvals(1:2)=[];

%%% Now reorder
Fn = F(idx,:);
%%% Conjugate
Fn(kvals<0,:)=conj(Fn(kvals<0,:));

%%% Similar for Zn
Zn = F(3:3:end,:);


    % Transition matrix. Operate on the same matrix in memory rather than
    % redefine each time.
    function T =  build_T_matrix_sub(AA,nn)
        T=zeros([nn nn]);
        ii = 1:(nn/3);
        T(3*nn*(ii-1)+3*(ii-1)+1)=AA(1);
        T(3*nn*(ii-1)+3*(ii-1)+2)=AA(2);
        T(3*nn*(ii-1)+3*(ii-1)+3)=AA(3);
        T(3*nn*ii-2*nn+3*(ii-1)+1)=AA(4);
        T(3*nn*ii-2*nn+3*(ii-1)+2)=AA(5);
        T(3*nn*ii-2*nn+3*(ii-1)+3)=AA(6);
        T(3*nn*ii-nn+3*(ii-1)+1)=AA(7);
        T(3*nn*ii-nn+3*(ii-1)+2)=AA(8);
        T(3*nn*ii-nn+3*(ii-1)+3)=AA(9);
    end

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
end
