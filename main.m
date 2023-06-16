% po odklejeniu elektrody, nie można porównywać amplitud
%todo
%centroidy
%liczba widm
% centroidy dla jednego ćw od wszystkich osób
% autokorelacja
% sprawdicz czygrupy się nie mies`ają
% roznicde dla posredniego BR-BR
% dists_cheby OR dists_chebyM
DEBUG = 1;
if(DEBUG)
    clear all;
    close all;
%     delete segments.mat % 10s
%     delete signals.mat  % hann window 10s
%     delete spectrums.mat%           
    delete centroids.mat
    DEBUG = 1;
end
allElapsedTime = tic;
windowing = 1;
printCentroids = 0; % a lot of console tables
plotAllFigures = 0;
deprecated = 0; % unused code
dirname = 'Archiwum';
dirname = '23.05.23';
dirname = '10';
clear v; % Clear files data
files = dir(fullfile(dirname,'**','*.mat'));
datafiles = fullfile({files.folder},{files.name});

txPr = "Pośredni"; txPc = "Podchwyt";
txBR = "Brachioradialis"; txBB = "Biceps brachii";
minActions = 10; skipedTrainingReppetinons = 0; % for segmentation selecting
k = 0;
for f = datafiles    
    k = k + 1;
    df = load(string(f));
    fieldname = fieldnames(df);
    vars = fieldnames(df);
    for i = 1:length(vars)
        assignin('base', "tmp", df.(vars{i}));        
        tmpData = tmp.movements.sources.signals.signal_1.data;
        tmpName = tmp.movements.sources.signals.signal_1.name;
        if ( tmp.movements.sources.signals.signal_1.name == "Ultium EMG-BRACHIORAD. RT")
            v(k).dataR = tmpData; % Radialis
            v(k).infoRecord = tmp.info.record_name;
            v(k).infoRName = tmpName;
            v(k).infoRDisp = txBR; 
        end
        if ( tmp.movements.sources.signals.signal_1.name == "Ultium EMG-BICEPS BR. RT")
            v(k).dataB = tmpData; % Biceps
            v(k).infoRecord = tmp.info.record_name;
            v(k).infoBName = tmpName;
            v(k).infoBDisp = txBB; 
        end
        
        tmpData = tmp.movements.sources.signals.signal_2.data;
        tmpName = tmp.movements.sources.signals.signal_2.name;
        if ( tmp.movements.sources.signals.signal_2.name == "Ultium EMG-BRACHIORAD. RT")
            v(k).dataR = tmpData; % Radialis
            v(k).infoRecord = tmp.info.record_name;
            v(k).infoRName = tmpName;
            v(k).infoRDisp = txBR; 
        end
        if ( tmp.movements.sources.signals.signal_2.name == "Ultium EMG-BICEPS BR. RT")
            v(k).dataB = tmpData; % Biceps
            v(k).infoRecord = tmp.info.record_name;
            v(k).infoBName = tmpName;
            v(k).infoBDisp = txBB; 
        end
    end
end
fprintf(1,"Liczba pacjentów: %d\n", length(v)/2); toc(allElapsedTime); % / liczba mięśni mierzonych równocześnie lub wykonywanego ćwiczenia
fprintf(1,"Macierz plików v: %d (ilość pomiarów * liczba ćwiczeń)\n", length(v));

% get from last sample
Yunits = tmp.movements.sources.signals.signal_1.units;
fpom = tmp.movements.sources.signals.signal_1.frequency; dtpom=1/fpom; %  Hz

if(exist("segments.mat"))
    load segments.mat
else
    nrF = 20; segmentActions;
end

if (skipedTrainingReppetinons) fprintf(1,"Rozkład: %d (pośrednich) %d (podchwytów) Sumarycznie pominiętych ćwiczeń: %d \n", countPosredni, countPodchwyt, skipedTrainingReppetinons); end;
fprintf(1,"Liczba segmentów: %d (liczba plików * liczba mięśni * ilość powtórzeń)\nSzacowanych segmentów: (%d)\n", length(segment), length(v)*2*10); % liczba plików * liczba mięśni * ilość powtórzeń
% figure;
% for i = 1:length(segment)
%     plot(segment(i).data)
% end

%--------------------------------------------------------------------------
nrFw = 201; % nr fig widma
lfrow= 4; lc=2; 
Tsyg=ceil(ceil(m*dtpom)/2)*2; % [sek] Maximal time of action duration
lSyg=round(Tsyg/dtpom);
clear Syg;
ksyg=0;     kol='kbrm'; kf=0;%figure(j)

if(exist("signals.mat"))
    load signals.mat
else
    
    tic
    for( i = 1:length(segment)) % uzupełnianie zerami segmentów
        SygRawLen(i) = length(segment(i).data');
        if (windowing)
            win = hann(SygRawLen(i));
            segment(i).data = segment(i).data.*win;
        end
        Syg(i,1:lSyg) = [segment(i).data' zeros(1, lSyg-length(segment(i).data))];
    
        SygKat(i) = segMio(i)+(segTraining(i)-1)*2; %plikSegMio(fileSegNr(nrs),nrB).i=2;  n+v(j).infoTraining-1*2; % training
        % 1 - Pośred BR, 2 Poś BB, 3 - Podchwyt BRadialis, 4 Podch BBiceps          
    end
    fprintf(1, "Okienkowanie Hanna ")
    toc;
    save signals.mat Syg SygKat SygRawLen
end

MTF(1).Tu = []; MTF(2).Tu = []; MTF(3).Tu = [];

if(exist("spectrums.mat"))
    load spectrums.mat
else
    nrF = 300; fprintf(1, "Liczenie widm... "); spectrumTrend;
end
fprintf(1,"Rozmiar widm: %dx%d\n", size(Widma)); 

if(exist("centroids.mat"))
    load centroids.mat
else
    nrF = 400; fprintf(1, "Liczenie centroidów... ");centroid;
end
fprintf(1,"Rozmiar centroidów: %dx%d\n", size(CentrWidm)); 

dCentr;

nrF = 800; bf = 0; disppolt; cc4; bf = 4; disppolt;
toc(allElapsedTime)
return

mnoznik = 24;
BR = v(11).dataR; BB = v(11).dataB; % radialis, biceps
L = resample(BR,mnoznik,1); R = resample(BB,mnoznik,1);
L = L/max(abs(L)); R = R/max(abs(R));
stereo_mtx = [L, R];%, v(11).dataB/max(abs(v(11).dataB))];
audiowrite('stereo sound normalized.wav', stereo_mtx, fpom*mnoznik);


% 11 to ja, 22 też
 X = []; 
for(i=1:length(v))
    if (v(i).infoRecord == """11 pośredni""")
        figure(i), subplot(211); plot(v(i).dataR); hold on; plot(v(i).dataB); axis('tight');
        Lf = length(v(i).dataB)/2; dtpom = 1/fpom*25;
        subplot(212); X = abs(fft(v(i).dataR)); plot([0:Lf-1]*dtpom, X(1:Lf)/Lf); hold on
        X = abs(fft(v(i).dataB)); plot([0:Lf-1]*dtpom, X(1:Lf)/Lf); hold off; axis('tight')
    end
    if (v(i).infoRecord == """22 pośredni NORAXON ELEKTRODY """)
        figure(i), subplot(211); plot(v(i).dataR); hold on; plot(v(i).dataB); axis('tight');
        Lf = length(v(i).dataB)/2; dtpom = 1/fpom*25;
        subplot(212); X = abs(fft(v(i).dataR)); plot([0:Lf-1]*dtpom, X(1:Lf)/Lf); hold on
        X = abs(fft(v(i).dataB)); plot([0:Lf-1]*dtpom, X(1:Lf)/Lf); hold off; axis('tight')
    end
    if (v(i).infoRecord == """22 podchwyt NORAXON ELEKTRODY""")
        figure(i), subplot(211); plot(v(i).dataR); hold on; plot(v(i).dataB); axis('tight');
        Lf = length(v(i).dataB)/2; dtpom = 1/fpom*25;
        subplot(212); X = abs(fft(v(i).dataR)); plot([0:Lf-1]*dtpom, X(1:Lf)/Lf); hold on
        X = abs(fft(v(i).dataB)); plot([0:Lf-1]*dtpom, X(1:Lf)/Lf); hold off; axis('tight')
    end
end


nrF = 900;
figureNubers = [134 136];
for( i = 200:100:nrF) figureNubers = [figureNubers i+16 i+17]; end
tic;
figPSW
toc;
% % 
% % figure(nrFw), subplot(1,2,1); hold off; subplot(1,2,2); hold off;
% % for (j = 1:length(v))
% %     dW=Widma(j,1).Ayf-Widma(j,2).Ayf; S(j,1)=sqrt(dW'*dW);
% %     dWwygl=wyglWidma(j,1).Af-wyglWidma(j,2).Af; S(j,2)=sqrt(dWwygl*dWwygl');
% %     Esr(j)=mean(Esyg(j,:)); % cecha
% %     lk = size(Widma,2); % liczba powtórzeń gestów j-tej kategrorii 
% %     Wsr(j).Wsr = zeros(1, length(Widma(j,1).Ayf'));
% %     for (k = 1:lk)
% %         Wsr(j).Wsr=Wsr(j).Wsr+Widma(j,k).Ayf';
% %     end
% %     Wsr(j).Wsr=Wsr(j).Wsr/lk;    %sredne widma % wzorzec
% % %     Wsr(j).Wsr=mean(Widma(j,:).Ayf')
% % %     Wsrwygl(j).Wsr=(wyglWidma(j,1).Af+wyglWidma(j,2).Af)/2; %S(j,2)=sqrt(dWwygl*dWwygl');
% %     minLenWidma = min(length(wyglWidma(j,1).Af),length(wyglWidma(j,1).Af2))
% % 
% %     Wsrwygl(j).Wsr = zeros(1, length(wyglWidma(j,1).Af));
% %     Wsrwygl(j).Wsr2 = zeros(1, length(wyglWidma(j,1).Af2));
% %     for (k = 1:lk)
% %         Wsrwygl(j).Wsr=Wsrwygl(j).Wsr+wyglWidma(j,k).Af;
% %         Wsrwygl(j).Wsr2=Wsrwygl(j).Wsr2+wyglWidma(j,k).Af2; % mocy
% % %         figure, plot(wyglWidma(j,k).Af2)
% %     end
% %     Wsrwygl(j).Wsr=Wsrwygl(j).Wsr/lk;    %sredne widma % wzorzec
% %     Wsrwygl(j).Wsr2=Wsrwygl(j).Wsr2/lk;    %sredne widma % wzorzec
% %     % D2
% %     dW2=Widma(j,1).Ayf2-Widma(j,2).Ayf2; S2(j,1)=sqrt(dW2'*dW2);
% %     dW2wygl=wyglWidma(j,1).Af2-wyglWidma(j,2).Af2; S2(j,2)=sqrt(dW2wygl*dW2wygl');
% % end
% % Ssr = [];
% % DWsr=Wsr(1).Wsr-Wsr(2).Wsr; Ssr(1)=sqrt(DWsr*DWsr');
% % DWsr=Wsrwygl(1).Wsr-Wsrwygl(2).Wsr; Ssr(2)=sqrt(DWsr*DWsr');
% % 
% %  figure(8), % dla widm nie wygładzonych
% %  for(j=1:2)
% %  subplot(2,2,j)
% %  plot([1:length(Wsr(j).Wsr)]/Tsyg, Wsr(j).Wsr,'r', ...
% %       [1:length(Widma(j,1).Ayf)]/Tsyg, Widma(j,1).Ayf,'c', ...
% %       [1:length(Widma(j,2).Ayf)]/Tsyg, Widma(j,2).Ayf,'k') % czerwone srednie widmo
% %  subplot(2,2,j+2)
% %  plot([1:length(Wsrwygl(j).Wsr)]/Tsyg, Wsrwygl(j).Wsr,'r', ...
% %       [1:length(wyglWidma(j,1).Af)]/Tsyg, wyglWidma(j,1).Af,'c', ...
% %       [1:length(wyglWidma(j,2).Af)]/Tsyg, wyglWidma(j,2).Af,'k') % czerwone srednie wygładzone widmo
% %  axis("tight"); xlabel("Hz")
% %  end
% %  dCG=0; dEG=0;
% % dists_cheby =[];
% nf=2;  
% for (j = 1:length(v) )
%     %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
%     for(k=1:lk) d=Wsrwygl(j).Wsr-wyglWidma(j,k).Af; 
%         dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d))/100; dists_cheby(j,k) = max(abs(d),[],2);  %uwaga przesunięcie przecinka
%         d2=Wsrwygl(j).Wsr2-wyglWidma(j,k).Af2; 
%         dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2))/100; dists_cheby2(j,k) = max(abs(d2),[],2);  %2-mocy
%         dEsyg(j,k)=abs(Esyg(j,k)-Esr(j));
%     end % odległość w grupie
%     figure(9), subplot(2,2,nf), plot(abs(d))
%     for(i=j+1:length(v))
%         dEsgr(j,i)=abs(Esr(i)-Esr(j))/2;
%         d=abs((Wsrwygl(j).Wsr-Wsrwygl(i).Wsr));  if(j==1) figure(9), subplot(1,2, 1), plot(d); end
%         dCG(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
%         dEG(j,i)=sqrt(sum(d.^2))/2; dists_chebyG(j,k) = max(abs(d),[],2)/2; % odległosv mięfzy grupowa
%         d=[];
%         d=abs((Wsrwygl(j).Wsr2-Wsrwygl(i).Wsr2));  if(j==1) figure(9), subplot(1,2, 1), plot(d); end
%         dCG2(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
%         dEG2(j,i)=sqrt(sum(d.^2))/2; dists_chebyG2(j,k) = max(abs(d),[],2)/2;
%         % dzielimy prrze 2 aby odl. m. grupowe były lepiej porównywanlne z
%         % liczonymi od centroidu
%         nf=4;
%     end
% end
% % dists_cheby
% dists_chebyG
% iloczyn wektorywy tylko w przestrzeni euclidesa
% miedzyGOrazodCentroidu = [[dEG; dCG ] dE dC ]
% miedzyGOrazodCentroiduGoraEuclidianBottomCity  = [[dEG(1,2);dCG(1,2)] dE(:,1) dC(:,1)]
% fprintf(1,"m.Group \t w.group")
% miedzyGOrazodCentroidu = [[dEG(1,2);dCG(1,2); dists_chebyG(1,2);dEsgr(1,2)] [dE(1,1);dC(1,1); dists_cheby(1,1);dEsyg(1,1)] [dE(2,1);dC(2,1);dists_cheby(2,1);;dEsyg(2,1)]]
% % miedzygrupowa
% % pomiędzygrpuami
% grupa 1, 2 , w gr. 1 2
% % stosunek odległośvi wewnątrz grupowej do mrfxygrupowej / druga jest spójniejsza/ w wierszach odległości E i C, w kolumnach grupy (osobh)
% % jakość klasyfikacji
% % można jeszcze na kwadratach spróbować 
% mg=[dEG(1,2);dCG(1,2);dists_chebyG(1,2);dEsgr(1,2)]; miedzyGOrazodCentroidu = [[dE(1,1);dC(1,1);dists_cheby(1,1);dEsyg(1,1)]./mg [dE(2,1);dC(2,1);dists_cheby(2,1);dEsyg(2,1)]./mg]
% mocy wygł
% mg2=[dEG2(1,2);dCG2(1,2);dists_chebyG2(1,2);dEsgr(1,2)]; miedzyGOrazodCentroidu2 = [[dE2(1,1);dC2(1,1);dists_cheby2(1,1);dEsyg(1,1)]./mg2 [dE2(2,1);dC2(2,1);dists_cheby2(2,1);dEsyg(2,1)]./mg2]
% Widmo = S;
% widmoMocy = sqrt(S2);
% porownywalnoscIrozroznialnosc=[S(1,:)./S(2,:);S(:,1)'./S(:,2)']
% po wierszach duża różnica po kolumnach mała
% Ssr
% figPSW jedakowe długości, wzorzec musi mieć takie same
%konkluzja wyrównywać zerami widma 
return
kol='kr';
for(i=1:2)
    ys=segment(i).data'; l = length(ys); tsym=l*dtpom; tsym=1;
    W(i).Ays = abs(fft(ys))*2/l; LAs=l;
    ys0=[ys  zeros(1,lSyg*16-l)]; l0=length(ys0);  tsym0=l*dtpom; tsym0=1;
    W(i).Ays0=abs(fft(ys0))*2/l0;
    figure(55), subplot(1,2,1); plot([0:LAs-1]/tsym,W(i).Ays(1:LAs),[kol(i) '.-']), hold on;
    subplot(1,2,2),plot([0:16*LAs-1]/tsym0,W(i).Ays0(1:16*LAs),[kol(i) '.-']); hold on;
end
subplot(1,2,1); hold off; subplot(1,2,2); hold off; 

N=1000; xy=2*sin(2*pi/N*[1:5*N])-3*sin(4*pi/N*[1:5*N])+4*sin(6*pi/N*[1:5*N]); AaN=abs(fft(xy))/(2.5*N); A2aN=AaN.^2; Aa2N=abs(fft(xy.^2))/(2.5*N);
figure(33); plot([0:49]/(3),A2a(1:50),'k.-',[0:49]/(3),Aa2(1:50),'r.-',[0:49]/(2.5),A2aN(1:50),'k.--',[0:49]/(2.5),Aa2N(1:50),'r.--');

% 
% m.Group 	 w.group
% miedzyGOrazodCentroidu =
% 
%    33.8508   20.6523   18.7710
%    15.3530    7.7939    4.1931
%     1.4859    1.4481    2.2938
%     4.7549    4.9129    1.0254
% 
% 
% miedzyGOrazodCentroidu = raw
% 
%     0.6101    0.5545
%     0.5076    0.2731
%     0.9746    1.5437
%     1.0332    0.2156
% 

% miedzyGOrazodCentroidu = lA
% 
%     0.6436    0.3943
%     0.5371    0.2292
%     0.9943    0.8487
%     1.0332    0.2156
% % 
% miedzyGOrazodCentroidu2 =
% 
%     0.7173    1.8997
%     0.7032    0.8733
%     0.8243    2.7688
%     1.0332    0.2156
% 
% 
% porownywalnoscIrozroznialnosc =
% 
%     0.9826    1.1002
%     2.6156    2.9287
% 
% 
% Ssr =
% 
%    92.8360   67.7016
% 
%    lA
%    


