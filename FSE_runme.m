%%% FSE simulations, performs Experiment 2 from the paper.
% Shaihan Malik 2016

%%% helper functions and quantities
psi = @(n)(2*pi*(0:n-1)/n);
ft= @(m)(circshift(fft(m,[],1),[floor(size(m,1)/2-1) 0])/size(m,1)); %<-- correct shift once lower negative F-states are dropped
ift = @(f)(ifft(circshift(f,[-floor(size(f,1)/2-1) 0]))*size(f,1));
n_indices = @(N)(-(2*N-2):2*N);
d2r = @(x)(x*pi/180);
r2d = @(x)(x*180/pi);

%%% Define sequence and relaxation properties
ESP=5;
T1=1500;
T2=500;
a0 = d2r([90 50*ones(1,100)]); %<-- Excitation & refocusing flips. Use 50 degrees to populate higher order configuration states

Npulse = length(a0);
Necho = Npulse-1;


%% Fully sampled Isochromats, vs EPG

%%% EPG
[ss, Fn,Zn,F] = FSE_EPG_sim(a0,'ESP',ESP,'T1',T1,'T2',T2);

%%% Isochromat summation
Niso = 4*Necho-1;
[ss0,mxy] = FSE_isochromat_sim(a0,Niso,'ESP',ESP,'T1',T1,'T2',T2,'psi',psi(Niso));

% Perform FTs
mf = ft(mxy);
Fs = ift(Fn);

figure(1);
clf
colormap_fade

nr=3;nc=2;
subplot(nr,nc,1:2)
pp=plot(abs([ss(:) ss0(:) 10*(ss(:)-ss0(:))]));
legend('EPG','Isochromat summation','Difference x10','location','southeast')
grid on
title('Predicted signal (|F_0|)')
xlabel('Echo number')
ylabel('|F_0| / M_0')

set(pp(1),'marker','*','markersize',3)

subplot(nr,nc,3)
imagesc(1:Necho,n_indices(Necho),abs(Fn),[0 0.1])
axis xy
title('EPG: F_n')
xlabel('Echo number')
ylabel('n')
colormap(jetfade)

subplot(nr,nc,4)
imagesc(1:Necho,[],abs(Fs))
title('ifft(F_n)')
xlabel('Echo number')
ylabel('m')
colormap(jet)

subplot(nr,nc,5)
imagesc(1:Necho,psi(Niso),abs(mxy))
xlabel('Echo number')
ylabel('\psi (rad)')
title('Isochromats: M_+')
colormap(jet)


subplot(nr,nc,6)
imagesc(1:Necho,n_indices(Necho),abs(ft(mxy)),[0 0.1])
axis xy
xlabel('Echo number')
ylabel('n')
title('fft(M_+)')
colormap(jetfade)

set(gcf,'position',[100 100 500 600])



%% Same simulation with reduced numbers of isochromats
Niso = 10:4*Necho-1;
nrmse = @(x1,x2)(norm(x1(:)-x2(:))/norm(x2(:)));

if 1
err = zeros(length(Niso));
mxy={};
for ii = 1:length(Niso)
    [sf,mxy{ii}] = FSE_isochromat_sim(a0,Niso(ii),'ESP',ESP,'T1',T1,'T2',T2,'psi',psi(Niso(ii)));
    err(ii)=nrmse(abs(sf(:)),abs(ss(:)));
end
end

%% Generate figure 4
figure(1)
clf
nr=3;nc=3;
subplot(nr,nc,5);
odix = find(mod(Niso,2)==1);
evix = find(mod(Niso,2)==0);
semilogy(Niso(odix),err(odix),'.-','markersize',10)
hold
semilogy(Niso(evix),err(evix),'.-','markersize',10)
legend('Odd number of isochromats','Even number of isochromats')
grid 
title('Error in predicted signal')
hx=xlabel('Number of isochromats, N');
ylabel('\epsilon','rotation',0,'fontsize',20)
hxpos=get(hx,'position');hxpos(1)=120;set(hx,'position',hxpos)
set(gca,'fontsize',14)

%%% Now generate separate frequency spectra diagrams
NisoR = [399 202 121 120 87];%<--- these are the plots to use
spidx = [6 8 7 2 1];%<- which subplot each of the above will appear in
for ix=1:length(spidx)
    subplot(nr,nc,spidx(ix))
    yt = -floor((NisoR(ix)-1)/2-2):floor((NisoR(ix)-1)/2);
    mf = ft(mxy{find(Niso==NisoR(ix))});
    imagesc(1:Necho,yt,abs(mf),[0 0.1])
    axis xy
    xlabel('Echo number')
    ylabel('n','rotation',0)
    title(sprintf('DFT(M_+): N=%d',NisoR(ix)))
    colormap(jetfade)
    grid on
    set(gca,'fontsize',11.5)
end

pp=get(gcf,'children');
set(gcf,'position',[100 100 1100 600],'paperpositionMode','auto')
cc = colorbar('southoutside');

set(pp(7),'position',[0.15 0.35 0.54 0.3])
% frequency subplot axes
ww=0.2;hh=0.2;
set(pp(1),'position',[0.07 0.75 ww hh])
set(pp(2),'position',[0.4 0.75 ww hh])
set(pp(3),'position',[0.07 0.05 ww hh])
set(pp(4),'position',[0.4 0.05 ww hh])
set(pp(5),'position',[1-ww-0.02 (1-hh*2)/2 ww hh*2])
set(cc,'position',[0.8 0.1 0.1 0.02],'fontsize',12)
aa=[];
aa{1} = annotation('arrow',[0.2136 0.3136],[0.2933 0.4900]);
aa{2} = annotation('arrow',[0.2418 0.2600],[0.7083 0.6400]);
aa{3} = annotation('arrow',[0.3900 0.3164],[0.7350 0.6400]);
aa{4} = annotation('arrow',[0.4591 0.4273],[0.2850 0.3517]);
aa{5} = annotation('arrow',[0.7455 0.6945],[0.4167 0.3633]);
for ii=1:5
    aa{ii}.HeadStyle = 'plain';
    aa{ii}.Color = [1 0 0];
end
text(377,-338,'Amplitude / M_0','fontsize',12)

print -dpng -r300 Figure4.png
