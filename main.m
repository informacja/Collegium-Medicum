DEBUG = 0;
dirname = 'Archiwum';
clear v
v.data = [];
v.info = [];
files = dir(fullfile(dirname,'**','*.mat'));
datafiles = fullfile({files.folder},{files.name});
k = 0;
for f = datafiles
    k = k + 1;
    df = load(string(f));
    fieldname = fieldnames(df);
    vars = fieldnames(df);
    for i = 1:length(vars)
        assignin('base', "tmp", df.(vars{i}));
        v(k).data = tmp.movements.sources.signals.signal_1.data;
        v(k).info = tmp.info.record_name;
    end
end
% get from last sample
Yunits = tmp.movements.sources.signals.signal_1.units;
fpom = tmp.movements.sources.signals.signal_1.frequency; dtpom=1/fpom; %  Hz
% posredni = record_2023_04_20_15_06_BR_po_redni.movements.sources.signals.signal_1.data;
% podchwyt = record_2023_04_20_15_07_BR_podchwyt.movements.sources.signals.signal_1.data;

% 
%     a(1).d = [v(1).data; v(2).data];
%     a(2).d = [v(2).data; v(1).data];
% v(1).data = a(1).d;
% v(2).data = a(2).d;

if(exist("segment.mat"))
    load segment.mat
else
    segmentActions
end
%--------------------------------------------------------------------------
lfrow= 4; lc=2; 
Tsyg=ceil(ceil(m*dtpom)/2)*2; % [sek]
lSyg=round(Tsyg/dtpom);
clear Syg;
ksyg=0;     kol='kbrm'; kf=0;%figure(j)
for( i = 1:length(segment))
    SygRawLen(i) = length(segment(i).data');
    Syg(i,1:lSyg) = [segment(i).data' zeros(1, lSyg-length(segment(i).data))];
%     figure(i), plot(Syg(i,:))
end

for(j = 1:length(v)) % grupa
    figure(j), nxf = 0; % nie drukuj figur
    fileSegNr; %TODO
    for (i = 1:length(find(fileSegNr==j))) % akcje w. grupy
        i = mod(i,2);
        figure(j)
        clear y;
        ksyg=ksyg+1;
        y = Syg(ksyg,:)'; %v(j).data(n1:Nbf);
        X = (y.^2); Esyg(j,i+1)=sum(X)*dtpom/lSyg; 
        %segment(nrs).data = y;
        Nf=length(y); %todo
        nx=[0:Nf-1];
        subplot(lfrow,lc,1+i),  plot(nx*dtpom, y); xlabel(sprintf("Ruch pośredni %d: y(t) t[sek]", i));
        if( i == 0) title("                                                                                                                                   Dziedzina czasu"); end
        ylabel(['Amplituda [' Yunits ']'])
        subplot(lfrow,lc,1+i+lc),  plot(nx*dtpom, X);
        % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Twygl=0.25; nTu = Tsyg/Twygl;
        Tu=Twygl/dtpom;
        run("../MTF/filtrWidma.m");
        xf = [0:LwAm-1];
        hold on; plot(xf/Tsyg,Ayf,'c',[0:Ldf]/Tsyg,Af,'k'); axis('tight');  hold off;
        xlabel(sprintf("Kwadrat y i wygł. y(t)^2 %d (Tu=%.1fms): t[sek]", i,Tu*dtpom*1000));
        A = fft(y); lA = length(A);
        Afw = abs(A(1:round(lA/2)))/SygRawLen(ksyg); %lA; 
        Podzial=4; if(j==2) Podzial=10; end
        nf=round(Nf/2); 
        X = Afw(1:nf);
        Twygl=0.05; nTu = Tsyg/Twygl; % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Tu=Twygl/dtpom;
        run("../MTF/filtrWidma.m");
        nk=round(Nf/Podzial);
        subplot(lfrow,lc,1+i+2*lc),   plot([0:nk-1]/Tsyg,Ayf(1:nk),'c',[0:Ldf]/Tsyg,Af,'k');
        Widma(j,i+1).Ayf=Ayf; wyglWidma(j,i+1).Af=Af;% i*2+j
        figure(nrFw), subplot(1,2,1); hold on; kf=mod(kf,4)+1; plot([0:Ldf],Af,kol(kf)); %plot(wyglWidma(j,i+1).Af); hold off; 
%         figPW("png")
        figure(j)
        if( 1+i+2*lc == 5 ) title("                                                                                                                                 Dziedzina częstotliwości"); end
        xlabel(sprintf("Widmo %d [Hz] Tu=%.1fms f_g=1/Tu=%.0fHz ",i,Tu*dtpom*1000,1/(Tu*dtpom)));
        Podzial=15; if(j==2) Podzial=30; end
        Twygl=0.025; nTu = Tsyg/Twygl;
%         nf=Nf;%round(Nf/Podzial);
%         Xx=fft(y.^2);Xx=abs(Xx(1:nf));  
        X=Afw(1:nf).^2;
%         figure(111), plot(X); hold on; plot(Xx)
        Tu=Twygl/dtpom;
        run("../MTF/filtrWidma.m");
        nf=round(Nf/Podzial);
        subplot(lfrow,lc,1+i+3*lc), plot([0:nf-1]/Tsyg,Ayf(1:nf),'c',[0:nf-1]/Tsyg,Af(1:nf),'k');
        xlabel(sprintf("Widmo mocy %d f_g=1/Tu Tu=%.1fms",i,Tu*dtpom*1000));
        Widma(j,i+1).Ayf2=Ayf; wyglWidma(j,i+1).Af2=Af;% i*2+j

        sgtitle( sprintf("%s %d",v(j).info, ksyg))
        figPW("png",5)
        figure(nrFw), subplot(1,2,2); hold on; plot([0:nf-1]/Tsyg,Af(1:nf),kol(kf)); %plot(wyglWidma(j,i+1).Af); hold off; 

        %nrs = nrs +1;
    end
end
figure(nrFw), subplot(1,2,1); hold off; subplot(1,2,2); hold off;
for (j = 1:length(v))
    dW=Widma(j,1).Ayf-Widma(j,2).Ayf; S(j,1)=sqrt(dW'*dW);
    dWwygl=wyglWidma(j,1).Af-wyglWidma(j,2).Af; S(j,2)=sqrt(dWwygl*dWwygl');
    Esr(j)=mean(Esyg(j,:));
    lk = size(Widma,2); % liczba powtórzeń gestów j-tej kategrorii 
    Wsr(j).Wsr = zeros(1, length(Widma(j,1).Ayf'));
    for (k = 1:lk)
        Wsr(j).Wsr=Wsr(j).Wsr+Widma(j,k).Ayf';
    end
    Wsr(j).Wsr=Wsr(j).Wsr/lk;    %sredne widma % wzorzec
%     Wsr(j).Wsr=mean(Widma(j,:).Ayf')
%     Wsrwygl(j).Wsr=(wyglWidma(j,1).Af+wyglWidma(j,2).Af)/2; %S(j,2)=sqrt(dWwygl*dWwygl');
    minLenWidma = min(length(wyglWidma(j,1).Af),length(wyglWidma(j,1).Af2))

    Wsrwygl(j).Wsr = zeros(1, length(wyglWidma(j,1).Af));
    Wsrwygl(j).Wsr2 = zeros(1, length(wyglWidma(j,1).Af2));
    for (k = 1:lk)
        Wsrwygl(j).Wsr=Wsrwygl(j).Wsr+wyglWidma(j,k).Af;
        Wsrwygl(j).Wsr2=Wsrwygl(j).Wsr2+wyglWidma(j,k).Af2; % mocy
%         figure, plot(wyglWidma(j,k).Af2)
    end
    Wsrwygl(j).Wsr=Wsrwygl(j).Wsr/lk;    %sredne widma % wzorzec
    Wsrwygl(j).Wsr2=Wsrwygl(j).Wsr2/lk;    %sredne widma % wzorzec
    % D2
    dW2=Widma(j,1).Ayf2-Widma(j,2).Ayf2; S2(j,1)=sqrt(dW2'*dW2);
    dW2wygl=wyglWidma(j,1).Af2-wyglWidma(j,2).Af2; S2(j,2)=sqrt(dW2wygl*dW2wygl');
end
Ssr = [];
DWsr=Wsr(1).Wsr-Wsr(2).Wsr; Ssr(1)=sqrt(DWsr*DWsr');
DWsr=Wsrwygl(1).Wsr-Wsrwygl(2).Wsr; Ssr(2)=sqrt(DWsr*DWsr');

 figure(8), % dla widm nie wygładzonych
 for(j=1:2)
 subplot(2,2,j)
 plot([1:length(Wsr(j).Wsr)]/Tsyg, Wsr(j).Wsr,'r', ...
      [1:length(Widma(j,1).Ayf)]/Tsyg, Widma(j,1).Ayf,'c', ...
      [1:length(Widma(j,2).Ayf)]/Tsyg, Widma(j,2).Ayf,'k') % czerwone srednie widmo
 subplot(2,2,j+2)
 plot([1:length(Wsrwygl(j).Wsr)]/Tsyg, Wsrwygl(j).Wsr,'r', ...
      [1:length(wyglWidma(j,1).Af)]/Tsyg, wyglWidma(j,1).Af,'c', ...
      [1:length(wyglWidma(j,2).Af)]/Tsyg, wyglWidma(j,2).Af,'k') % czerwone srednie wygładzone widmo
 axis("tight"); xlabel("Hz")
 end
%  dCG=0; dEG=0;
% dists_cheby =[];
nf=2;  
for (j = 1:length(v) )
    %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
    for(k=1:lk) d=Wsrwygl(j).Wsr-wyglWidma(j,k).Af; 
        dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d))/100; dists_cheby(j,k) = max(abs(d),[],2);  %uwaga przesunięcie przecinka
        d2=Wsrwygl(j).Wsr2-wyglWidma(j,k).Af2; 
        dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2))/100; dists_cheby2(j,k) = max(abs(d2),[],2);  %2-mocy
        dEsyg(j,k)=abs(Esyg(j,k)-Esr(j));
    end % odległość w grupie
    figure(9), subplot(2,2,nf), plot(abs(d))
    for(i=j+1:length(v))
        dEsgr(j,i)=abs(Esr(i)-Esr(j))/2;
        d=abs((Wsrwygl(j).Wsr-Wsrwygl(i).Wsr));  if(j==1) figure(9), subplot(1,2, 1), plot(d); end
        dCG(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
        dEG(j,i)=sqrt(sum(d.^2))/2; dists_chebyG(j,k) = max(abs(d),[],2)/2; % odległosv mięfzy grupowa
        d=[];
        d=abs((Wsrwygl(j).Wsr2-Wsrwygl(i).Wsr2));  if(j==1) figure(9), subplot(1,2, 1), plot(d); end
        dCG2(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
        dEG2(j,i)=sqrt(sum(d.^2))/2; dists_chebyG2(j,k) = max(abs(d),[],2)/2;
        % dzielimy prrze 2 aby odl. m. grupowe były lepiej porównywanlne z
        % liczonymi od centroidu
        nf=4;
    end
end
% dists_cheby
% dists_chebyG
% iloczyn wektorywy tylko w przestrzeni euclidesa
miedzyGOrazodCentroidu = [[dEG; dCG ] dE dC ]
miedzyGOrazodCentroiduGoraEuclidianBottomCity  = [[dEG(1,2);dCG(1,2)] dE(:,1) dC(:,1)]
fprintf(1,"m.Group \t w.group")
miedzyGOrazodCentroidu = [[dEG(1,2);dCG(1,2); dists_chebyG(1,2);dEsgr(1,2)] [dE(1,1);dC(1,1); dists_cheby(1,1);dEsyg(1,1)] [dE(2,1);dC(2,1);dists_cheby(2,1);;dEsyg(2,1)]]
% miedzygrupowa
% pomiędzygrpuami
grupa 1, 2 , w gr. 1 2
% stosunek odległośvi wewnątrz grupowej do mrfxygrupowej / druga jest spójniejsza/ w wierszach odległości E i C, w kolumnach grupy (osobh)
% jakość klasyfikacji
% można jeszcze na kwadratach spróbować 
mg=[dEG(1,2);dCG(1,2);dists_chebyG(1,2);dEsgr(1,2)]; miedzyGOrazodCentroidu = [[dE(1,1);dC(1,1);dists_cheby(1,1);dEsyg(1,1)]./mg [dE(2,1);dC(2,1);dists_cheby(2,1);dEsyg(2,1)]./mg]
mocy wygł
mg2=[dEG2(1,2);dCG2(1,2);dists_chebyG2(1,2);dEsgr(1,2)]; miedzyGOrazodCentroidu2 = [[dE2(1,1);dC2(1,1);dists_cheby2(1,1);dEsyg(1,1)]./mg2 [dE2(2,1);dC2(2,1);dists_cheby2(2,1);dEsyg(2,1)]./mg2]
Widmo = S;
widmoMocy = sqrt(S2);
porownywalnoscIrozroznialnosc=[S(1,:)./S(2,:);S(:,1)'./S(:,2)']
% po wierszach duża różnica po kolumnach mała
Ssr
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


