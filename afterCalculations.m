% for Ludwin
figure(1), close 1
fpom; % sampl. freq
lSyg; % signal length after zero-padding
segment(1).kat;

v; %surowe
Widma; 
wyglWidma;
CentrWidm;
dists_chebyM;

CC;
CCE;

wybrJ; % selectet ID in v matrix with author recorded muscles

%%
% eg. plots
figure(1),nexttile,
plot(v(1).dataB), hold on, plot(v(1).dataR), hold off; legend
nexttile, 
L = lSyg; % signal length
f = fpom/L*(0:(L/2-1));
wl = length(wyglWidma(1,1).Af);
plot(f,Widma(1,1).Ayf, "DisplayName","Widma"); hold on,
plot(f(1:wl),wyglWidma(1,1).Af, "DisplayName","wyglWidma"); 
plot(f(1:wl),CentrWidm(1,2).AfM,"DisplayName","CentrWidm"); hold off; 
legend; axis tight
xlim([0 150])

%  OR from fig Library
% x = v(wybrJ(1)).dataB;
% figFFTplot(x, fpom);