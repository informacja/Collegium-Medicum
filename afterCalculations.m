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

%% Figure 5
figure(), txStd = "std"; txAS = "ave. spect.";
ww=[]; vNum = wybrJ(1);
for(i = 1:10)%size(wyglWidma,2))
    ww(i,:) = wyglWidma(vNum,i).Af2/wyglWidma(vNum,i).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("a)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);
% 

ww=[]; vNum = wybrJ(2);
for(i = 1:10)
    ww(i,:) = wyglWidma(vNum,i).Af2/wyglWidma(vNum,i).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("b)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);


ww=[]; vNum = wybrJ(1); ip = 11; ik=20;
for(ii = ip:ik)
    i = ii-ip+1;
    ww(i,:) = wyglWidma(vNum,ii).Af2/wyglWidma(vNum,ii).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("c)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);


ww=[]; vNum = wybrJ(2);
for(ii = ip:ik)
    i = ii-ip+1;
    ww(i,:) = wyglWidma(vNum,ii).Af2/wyglWidma(vNum,ii).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("d)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);


% second row norm by energy
ww=[]; vNum = wybrJ(1);
for(i = 1:10)%size(wyglWidma,2))
    ww(i,:) = wyglWidma(vNum,i).Af2/sum(wyglWidma(vNum,i).Af2);
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("e)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);

ww=[]; vNum = wybrJ(2);
for(i = 1:10)
    ww(i,:) = wyglWidma(vNum,i).Af2/sum(wyglWidma(vNum,i).Af2);
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("f)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);


ww=[]; vNum = wybrJ(1);
for(ii = ip:ik)
    i = ii-ip+1;
    ww(i,:) = wyglWidma(vNum,ii).Af2/sum(wyglWidma(vNum,ii).Af2);
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("g)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);

ww=[]; vNum = wybrJ(2);
for(ii = ip:ik)
    i = ii-ip+1;
    ww(i,:) = wyglWidma(vNum,ii).Af2/sum(wyglWidma(vNum,ii).Af2);
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("h)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);
hold on



 width = 13; % inches
 hight =  6.5000;
 set(gcf,'Units','inches');                        % jednostka wymiarowania okna
 set(gcf,'Position', [0 0 width*2 hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
           

% check by plots
figure(), txStd = "std"; txAS = "ave. spect."; zakres = 11:20;
ww=[]; vNum = wybrJ(1);
for(i = zakres)%size(wyglWidma,2))
    ww(i,:) = Widma(vNum,i).Ayf;
end
nexttile, plot(f(1:length(ww)),smoothdata( ww)), axis tight; subtitle("a)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); 
legend(string(zakres));
% 





ww=[]; vNum = wybrJ(2);
for(i = 1:10)
    ww(i,:) = wyglWidma(vNum,i).Af2/wyglWidma(vNum,i).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("b)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);


ww=[]; vNum = wybrJ(1);
for(i = 11:20)
    ww(i,:) = wyglWidma(vNum,i).Af2/wyglWidma(vNum,i).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("c)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);


ww=[]; vNum = wybrJ(2);
for(i = 11:20)
    ww(i,:) = wyglWidma(vNum,i).Af2/wyglWidma(vNum,i).maxAf2;
end
nexttile, stdshade(ww,0.2,'g',f(1:length(ww))), axis tight; subtitle("d)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]"); xlim([0 250]); legend(txStd, txAS);

%% Figure 6
figure(6)
ww=[];
for(i = 1:size(CentrWidm,1))
    ww(i,:) = CentrWidm(i,1).Af2M;
end
nexttile, stdshade(ww(:,:),0.2,'b',f(1:length(ww))),
axis tight; subtitle("a)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250]);legend(txStd, txAS);

ww=[];
for(i = 1:size(CentrWidm,1))
    ww(i,:) = CentrWidm(i,2).Af2M;
end
nexttile, stdshade(ww(:,:),0.2,'b',f(1:length(ww))),
axis tight; subtitle("b)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250]);legend(txStd, txAS);

ww=[];
for(i = 1:size(CentrWidm,1))
    ww(i,:) = CentrWidm(i,1).Af2E;
end
nexttile, stdshade(ww(:,:),0.2,'b',f(1:length(ww))),
axis tight; subtitle("c)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250]);legend(txStd, txAS);

ww=[];
for(i = 1:size(CentrWidm,1))
    ww(i,:) = CentrWidm(i,2).Af2E;
end
nexttile, stdshade(ww(:,:),0.2,'b',f(1:length(ww))),
axis tight; subtitle("d)"); xlabel("Frequency [Hz]"); ylabel("Power [a. u.]");xlim([0 250]);legend(txStd, txAS);
hold on

%  OR from fig Library
% x = v(wybrJ(1)).dataB;
% figFFTplot(x, fpom);