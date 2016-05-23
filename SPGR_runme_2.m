%% SPGR diffusion script. Generates figure 5 from ISMRM abstract
% Shaihan Malik, October 2015

TR=5;
alpha= 10;
phi_0= 120;
T1=1500;
T2=500;
npulse = floor(5*T1/TR);
nmax=npulse-1;


%%% helper functions and quantities
psi = @(n)(2*pi*(0:fix(n)-1)/fix(n));
ft = @(m)(fftshift(fft(m,[],1),1)/size(m,1));
ift = @(f)(ifft(ifftshift(f,1),[],1)*size(f,1));
n_indices = @(Niso)(-floor((Niso)/2):floor((Niso-1)/2));
d2r = @(x)(x*pi/180);
r2d = @(x)(x*180/pi);
colormap_fade


%% Diffusion: Specify gradient strengths and durations (note the structures 
% below can take periods of zero gradient value, not included here)
% Both EPG and Isochromat codes take the same structure

G = [-2.8 6.1 3.2]; % mT/m
tau = [1 2 3.3]; %ms
D = 3e-9; % m^2/s

diff=struct;
diff.G = G; % mT/m ms
diff.tau = tau; %ms
diff.D = D;




%% Perform calculation for multiple RF spoiling phase angles. Perform isochromat
% calculation with Niso=100, and do EPG and Isochromats both with and without diffusion.

phi_arr = 1:0.5:180;
nphi=length(phi_arr);
Niso=100;
npulse = floor(5*T1/TR);

if 0 %<-- this loop takes ages to run, hence I've stored results. But feel free to run it again
    %%% Store signals
    Sig = zeros(nphi,4);
    figure(1)
    clf
    for ii=1:nphi
        [tmp0,Fn] = SPGR_EPG_sim(d2r(alpha), d2r(phi_arr(ii)),TR, T1, T2,npulse);
        [tmp1,mxy] = SPGR_isochromat_sim(d2r(alpha),d2r(phi_arr(ii)),TR,T1,T2,npulse,'Niso',Niso);
        [tmp2,Fn2] = SPGR_EPG_sim(d2r(alpha), d2r(phi_arr(ii)),TR, T1, T2,npulse,'diff',diff);
        [tmp3,mxy2] = SPGR_isochromat_sim(d2r(alpha),d2r(phi_arr(ii)),TR,T1,T2,npulse,'Niso',Niso,'diff',diff);
        
        Sig(ii,:) = abs([tmp0(end) tmp1(end) tmp2(end) tmp3(end)]);
        disp([ii nphi])
        %save Sig Sig
        plot(Sig,'-')
        drawnow
        pause(0.0001)
        save Sig Sig
    end
else 
    load Sig
end

%% Generate data for selected phase angle
phisel=108;
[tmp0,Fn] = SPGR_EPG_sim(d2r(alpha), d2r(phisel),TR, T1, T2,npulse,'kmax',inf);
[tmp1,mxy] = SPGR_isochromat_sim(d2r(alpha),d2r(phisel),TR,T1,T2,npulse,'Niso',Niso);
[tmp2,Fn2] = SPGR_EPG_sim(d2r(alpha), d2r(phisel),TR, T1, T2,npulse,'diff',diff,'kmax',inf);
[tmp3,mxy2] = SPGR_isochromat_sim(d2r(alpha),d2r(phisel),TR,T1,T2,npulse,'Niso',Niso,'diff',diff);


%%
Sideal = sind(alpha).*(1-exp(-TR./T1))./(1-exp(-TR./T1).*cosd(alpha));

figure(1);clf
nr = 2;nc=3;
fs=15;fs2=13;
subplot(nr,nc,1:2);
plot(phi_arr,Sig(:,1:2)/Sideal);
grid
hold
plot(phi_arr,abs(Sig(:,1)-Sig(:,2))/Sideal);

legend('EPG','100 isochromats','difference')
xlabel('\Phi / deg')
ylabel('|F_0| / S_{ideal}')
title('No diffusion')
set(gca,'fontsize',fs)

subplot(nr,nc,4:5);
plot(phi_arr,Sig(:,3:4)/Sideal);
grid
hold
plot(phi_arr,abs(Sig(:,3)-Sig(:,4))/Sideal);
legend('EPG','100 isochromats','difference')
xlabel('\Phi / deg')
ylabel('|F_0| / S_{ideal}')
title('With diffusion')
set(gca,'fontsize',fs)


%%% images
yl=[-150 150];
winl=[0 0.1];
winlog=[-5 0];
subplot(nr,nc,3)
imagesc(1:npulse,n_indices(npulse),log10(abs(Fn)/Sideal),winlog);ylim(yl)
axis xy
colormap(jetfade)
hold
patch([1 1 npulse npulse],[-50 yl(1) yl(1) -50],[0 0 0],'facealpha',0.3,'edgealpha',0.)
patch([1 1 npulse npulse],[50 yl(2) yl(2) 50],[0 0 0],'facealpha',0.3,'edgealpha',0.)
title('EPG, \Phi=108^\circ No diffusion')
xlabel('Pulse number')
ylabel('n','rotation',0,'fontweight','bold')
set(gca,'fontsize',fs)
set(gca,'fontsize',fs2)

subplot(nr,nc,6)
imagesc(1:npulse,n_indices(npulse),log10(abs(Fn2)/Sideal),winlog);ylim(yl)

axis xy
colormap(jetfade)
hold
patch([1 1 npulse npulse],[-50 yl(1) yl(1) -50],[0 0 0],'facealpha',0.3,'edgealpha',0.)
patch([1 1 npulse npulse],[50 yl(2) yl(2) 50],[0 0 0],'facealpha',0.3,'edgealpha',0.)
title('EPG, \Phi=108^\circ With diffusion')
xlabel('Pulse number')
ylabel('n','rotation',0,'fontweight','bold')
set(gca,'fontsize',fs2)

cc = colorbar;
set(cc,'position',[0.955 0.34 0.018 0.3],'fontsize',fs2)

ha = annotation('arrow',[0.4000 0.4203],[0.8745 0.8381]);
set(gcf,'position',[300 300 900 550])

%% Make a bunch of different EPGs for different phase offsets

theta_vec = [10 10 10 45 45 45];
phi_vec = [5 60 117 5 60 117];

Niso = 40;
npulse=100;
nrmse = @(x1,x2)(norm(x1(:)-x2(:))/norm(x2(:)));

Fn={};
err=[];
for ii=1:6
    [s,Fn{ii}] = SPGR_EPG_sim(d2r(theta_vec(ii)),d2r(phi_vec(ii)),TR, T1, ...
        T2,npulse,'kmax',inf);  
    [sstmp,mxytmp] = SPGR_isochromat_sim(d2r(theta_vec(ii)),d2r(phi_vec(ii)),...
        TR, T1, T2,npulse,'psi',psi(Niso));
   
     err(ii) = nrmse(sstmp,s);
end

figfp(4)
nr=2;nc=3;
for ii=1:6
    subplot(nr,nc,ii)
    %imagesc(1:npulse,n_indices(2*nmax+1),abs(Fn{ii}),[0 0.05])
    imagesc(1:npulse,n_indices(2*npulse-1),log10(abs(Fn{ii})),[-5 -1])
    axis xy
    title(sprintf('F_n     \\theta = %d^\\circ \\Phi = %1.1f ^\\circ',theta_vec(ii),phi_vec(ii)))
    xlabel('Pulse number')
    ylabel('n')
    
    text(1,-80,sprintf('nRMS diff = %1.3f',err(ii)),'color',[0 0 0],'fontweight','bold','fontsize',13)
end

colormap(jetfade)

set(gcf,'position',[100 200 800 450])


