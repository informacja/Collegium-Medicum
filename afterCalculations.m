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

Sb;

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

%%
figure
ww=[]; vNum = wybrJ(1);
for(i = 1:10)%size(wyglWidma,2))
    ww(i,:) = wyglWidma(vNum,i).Af/wyglWidma(vNum,i).maxAf;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight
% hold on; plot(f(1:length(ww)),CentrWidm(vNum).Af2M); 
subtitle("a)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250])

ww=[]; vNum = wybrJ(2);
for(i = 1:10)%size(wyglWidma,2))
ww(i,:) = wyglWidma(vNum,i).Af/wyglWidma(vNum,i).maxAf;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("b)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250])

ww=[]; vNum = wybrJ(1);
for(i = 11:20)%size(wyglWidma,2))
ww(i,:) = wyglWidma(vNum,i).Af/sum(wyglWidma(vNum,i).Af);
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("c)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250])

ww=[]; vNum = wybrJ(2);
for(i = 11:20)%size(wyglWidma,2))
ww(i,:) = wyglWidma(vNum,i).Af/sum(wyglWidma(vNum,i).Af);
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("d)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250]);
hold on
 % plot(f(1:length(ww)),CentrWidm(vNum).AfE)
% nexttile
% ww(1,:) = wyglWidma(2,:).Af

% stdshade(ww,0.2,'g')

%  OR from fig Library
% x = v(wybrJ(1)).dataB;
% figFFTplot(x, fpom);