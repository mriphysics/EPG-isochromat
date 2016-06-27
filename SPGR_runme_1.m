%%% SPGR simulations, performs Experiment 1 from the paper.
% Shaihan Malik 2016

TR=5;
alpha= 10;
phi_0= 120;
T1=1500;
T2=500;
npulse = 100;
nmax=npulse-1;


%%% helper functions and quantities
psi = @(n)(0.1*randn(1,n)+2*pi*(0:fix(n)-1)/fix(n));
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
ylabel('n','rotation',0)
colormap(jetfade)


subplot(nr,nc,4)
imagesc(1:npulse,1:2*nmax+1,abs(Fs))
title('iDFT(F_n)')
xlabel('Pulse number')
ylabel('m','rotation',0)
colormap(jet)


subplot(nr,nc,5)
imagesc(1:npulse,psi(2*nmax+1),abs(mxy0))
xlabel('Pulse number')
ylabel('\psi_m (rad)')
title('Isochromats: M_+')
colormap(jet)


subplot(nr,nc,6)
imagesc(1:npulse,-nmax:nmax,abs(ft(mxy0)),[0 0.1])
axis xy
xlabel('Pulse number')
ylabel('n','rotation',0)
title('DFT(M_+)')
colormap(jetfade)


set(gcf,'position',[100 100 500 600],'paperpositionmode','auto')

%%% Annotate
gg=get(gcf,'Children');
axes(gg(6))
text(-10,-0.01,'(a)','fontsize',18,'fontweight','bold')
axes(gg(4))
text(-20,-110,'(b)','fontsize',18,'fontweight','bold')
axes(gg(3))
text(-20,220,'(c)','fontsize',18,'fontweight','bold')
axes(gg(2))
text(-20,2*pi+0.2,'(d)','fontsize',18,'fontweight','bold')
axes(gg(1))
text(-20,-110,'(e)','fontsize',18,'fontweight','bold')

print -dpng -r300 Figure2.png








%% 2. More reduced versions

Niso = [10:5:npulse 2*npulse-1 2*npulse+50];

nrmse = @(x1,x2)(norm(x1(:)-x2(:))/norm(x2(:)));


[s,Fn,Zn] = SPGR_EPG_sim(d2r(alpha), d2r(phi_0),TR, T1, T2,npulse,'kmax',inf);

Fxy={};
mxy={};
s_iso=zeros([npulse length(Niso) ]);
err=zeros([length(Niso) 1]);
firsterr=zeros([length(Niso) 1]);
for ii=1:length(Niso)
    % Compute signal for reduced number of isochromats
    [sstmp,mxytmp] = SPGR_isochromat_sim(d2r(alpha),d2r(phi_0),TR, T1, T2,npulse,'psi',psi(Niso(ii)));

    s_iso(:,ii)=sstmp;
    err(ii) = nrmse(sstmp,s);
    Fxy{ii} = ft(mxytmp);
    mxy{ii}=mxytmp;
    
    %%% How many FIDs are correctly predicted?
    erridxtmp=find(abs(abs(sstmp(:))-abs(s(:)))>1e-10,1);%<- first non-zero error
    if ~isempty(erridxtmp)
        firsterr(ii)=erridxtmp;
    else
        firsterr(ii)=inf;
    end
end

%% Display

figure(3);
clf
nr=3;nc=4;

subplot(nr,nc,11:12);
semilogy(Niso,err);
grid 
xlabel('Number of isochromats (N)')
ylabel('S_{err}','rotation',0)
title('Error in predicted signal')

subplot(nr,nc,3:4);
ix = [ 7 19 21];
plot(abs(s_iso(:,ix)))
grid
hold
plot(abs(s),'k')
leg={};
for ii=1:length(ix),leg{ii}=sprintf('N = %d',Niso(ix(ii)));end
leg{end+1} = 'EPG';
ll=legend(leg);%,'location','eastoutside');
ylim([0 0.2])
xlabel('RF pulse number')
ylabel('|F_0|/M_0')
title('Simulated signal (|F_0|) vs N')

%%% plots

for ii=1:length(ix)
    
    subplot(nr,nc,nc*(ii-1)+1)
    imagesc(1:npulse,psi(Niso(ix(ii))),abs(mxy{ix(ii)}),[0 1])
    axis xy
    title(sprintf('M_+: N=%d',Niso(ix(ii))))
    ylabel('\psi_m (rad)')
    xlabel('pulse number')
    
    subplot(nr,nc,nc*(ii-1)+2)
    imagesc(1:npulse,n_indices(Niso(ix(ii))),abs(Fxy{ix(ii)}),[0 0.1])
    axis xy
    title(sprintf('DFT(M_+): N=%d',Niso(ix(ii))))
    ylabel('n','rotation',0)
    xlabel('pulse number')
end

colormap(jetfade)

set(gcf,'position',[100 100 900 500],'paperpositionmode','auto')

%%% Annotate
gg=get(gcf,'Children')

axes(gg(1))
cc=colorbar;
cc.Location = 'southoutside';
cc.Position = [0.38 0.04 0.12 0.022];

axes(gg(2))
cc=colorbar;
cc.Location = 'southoutside';
cc.Position = [0.11 0.04 0.12 0.022];


%
% ll = {'(f)','(e)','(d)','(c)','(b)','(a)'};
ll = {'(f)','(c)','(e)','(b)','(d)','(a)'};
for ii=[1 3 5]
    pos = gg(ii).Position;
    gg(ii).Position = [0.34 pos(2)+0.03 0.19 0.19];
    axes(gg(ii))
    text(-25,-Niso(ix(4-(ii+1)/2))*0.52,ll{ii},'fontsize',18,'fontweight','bold')
end

for ii=[2 4 6]
    pos = gg(ii).Position;
    gg(ii).Position = [0.07 pos(2)+0.03 0.19 0.19];
     axes(gg(ii))
    text(-25,-0.01,ll{ii},'fontsize',18,'fontweight','bold')
end

gg(8).Position = [0.62 0.6 0.33 0.32];
axes(gg(8))
text(-18,-0.01,'(g)','fontsize',18,'fontweight','bold')

gg(9).Position = [0.62 0.15 0.33 0.32];
axes(gg(9))
text(-50,5e-16,'(h)','fontsize',18,'fontweight','bold')

print -dpng -r300 Figure3.png
