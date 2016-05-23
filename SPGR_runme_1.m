%%% SPGR, NO diffusion. This generates figures 2&3 from abstract
% Shaihan Malik, October 2015

TR=5;
alpha= 10;
phi_0= 120;
T1=1500;
T2=500;
npulse = 100;
nmax=npulse-1;


%%% helper functions and quantities
psi = @(n)(2*pi*(0:fix(n)-1)/fix(n));
ft = @(m)(fftshift(fft(m,[],1),1)/size(m,1)); %<- shifted form of DFT with 1/N normalisation
ift = @(f)(ifft(ifftshift(f,1),[],1)*size(f,1));
n_indices = @(Niso)(-floor((Niso)/2):floor((Niso-1)/2));
d2r = @(x)(x*pi/180);
r2d = @(x)(x*180/pi);

%%% 1. Simulate with sufficient number of isochromats
[s,Fn,Zn] = SPGR_EPG_sim(d2r(alpha), d2r(phi_0),TR, T1, T2,npulse,'kmax',inf);

Fs = ift(Fn);

[ss0,mxy0,mz0] = SPGR_isochromat_sim(d2r(alpha),d2r(phi_0),TR, T1, T2,npulse,'psi',psi(2*nmax+1));


%%% Get colormap
colormap_fade
%
figure(1);
clf
nr=3;nc=2;
subplot(nr,nc,1:2)
pp=plot(abs([s(:) ss0(:) 10*(s(:)-ss0(:))]));
legend('EPG','Isochromat summation','Difference x10')
grid on
title('Predicted signal (|F_0|)')
xlabel('Pulse number')
ylabel('|F_0| / M_0')
set(pp(1),'marker','*','markersize',3)

subplot(nr,nc,3)
imagesc(1:npulse,-nmax:nmax,abs(Fn),[0 0.1])
axis xy
title('EPG: F_n')
xlabel('Pulse number')
ylabel('n')
colormap(jetfade)


subplot(nr,nc,4)
imagesc(1:npulse,1:2*nmax+1,abs(Fs))
title('ifft(F_n)')
xlabel('Pulse number')
ylabel('m')
colormap(jet)


subplot(nr,nc,5)
imagesc(1:npulse,psi(2*nmax+1),abs(mxy0))
xlabel('Pulse number')
ylabel('\psi (rad)')
title('Isochromats: M_+')
colormap(jet)


subplot(nr,nc,6)
imagesc(1:npulse,-nmax:nmax,abs(ft(mxy0)),[0 0.1])
axis xy
xlabel('Pulse number')
ylabel('n')
title('fft(M_+)')
colormap(jetfade)


set(gcf,'position',[100 100 500 600],'name','Abstract - Figure 2')

%% 2. More reduced versions

Niso = [10:5:npulse 2*npulse-1 2*npulse+50];

nrmse = @(x1,x2)(norm(x1(:)-x2(:))/norm(x2(:)));


[s,Fn,Zn] = SPGR_EPG_sim(d2r(alpha), d2r(phi_0),TR, T1, T2,npulse,'kmax',inf);

Fxy={};
s_iso=zeros([npulse length(Niso) ]);
err=zeros([length(Niso) 1]);
firsterr=zeros([length(Niso) 1]);
for ii=1:length(Niso)
    % Compute signal for reduced number of isochromats
    [sstmp,mxytmp] = SPGR_isochromat_sim(d2r(alpha),d2r(phi_0),TR, T1, T2,npulse,'psi',psi(Niso(ii)));

    s_iso(:,ii)=sstmp;
    err(ii) = nrmse(sstmp,s);
    Fxy{ii} = ft(mxytmp);
    
    %%% How many FIDs are correctly predicted?
    erridxtmp=find(abs(abs(sstmp(:))-abs(s(:)))>1e-10,1);%<- first non-zero error
    if ~isempty(erridxtmp)
        firsterr(ii)=erridxtmp;
    else
        firsterr(ii)=inf;
    end
end

%%% Display

figure(3);
clf
nr=3;nc=2;

%subplot(nr,nc,[1 3]);
subplot(nr,nc,(nc)+2*nc);
semilogy(Niso,err);
grid 
xlabel('Number of isochromats (N)')
ylabel('RMS difference in |F_0|')
title('RMS signal deviation from EPG')

%subplot(nr,nc,[7 9]);
subplot(nr,nc,(1)+2*nc);
ix = [3 7 12 17];
plot(abs(s_iso(:,ix)))
grid
hold
plot(abs(s),'k')
leg={};
for ii=1:length(ix),leg{ii}=sprintf('N = %d',Niso(ix(ii)));end
leg{end+1} = 'EPG';
ll=legend(leg);%,'location','eastoutside');
% ylim([0 0.35])
xlabel('RF pulse number')
ylabel('|F_0|/M_0')
title('Simulated signal (|F_0|) vs N')

%%% plots
ix = ([20 21 19 7]);
for ii=1:length(ix)
    %subplot(nr,nc,2*ii)
    subplot(nr,nc,ii)
    imagesc(1:npulse,n_indices(Niso(ix(ii))),abs(Fxy{ix(ii)}),[0 0.1])
    axis xy
    title(sprintf('FFT(M_+): N=%d',Niso(ix(ii))))
    ylabel('n')
    xlabel('pulse number')
end

colormap(jetfade)

set(gcf,'position',[100 100 865 650],'name','Abstract - Figure 3')
