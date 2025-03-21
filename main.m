% po odklejeniu elektrody, nie można porównywać amplitud
%todo
% centroidy dla jednego ćw od wszystkich osób
% roznicde dla posredniego BR-BR
% dists_cheby OR dists_chebyM

DEBUG = 1;
forArticle = 1; % save figures
% Parseval = 1; % setted from generate....m
% -1 = ujemny relatywne Twygl względem length(segment)
%  0 = zero nie przeliczamy Twygl i segmenty bez zer
%  1 = dodaj zera
PL = 1;
EN = 2;
lang = EN;
if(isdeployed)
    DEBUG = 0;
    % Parseval = -1; % 1  dodaj zera, ujemy relatywne Twygl wzgledem lentth(segment), zero nie przeliczamy Twygl i segmenty bez zer
end
if(forArticle)
    delete spectrums.mat
    delete centroids.mat
    delete figBase
end
if(DEBUG)
    % clear all;
    % close all;
    % delete segments.mat % 10s
    % delete signals.mat  % hann window 10s
    % delete spectrums.mat           
    % delete centroids.mat
    DEBUG = 1;
    if(~exist("Parseval","var"))
        Parseval = -1;
    end
    % Parseval = 0; % 0 - nie dodawaj zer
end

fWyswieltCentroidow = 555; % [Hz]

global cntCol;

allElapsedTime = tic;
windowing = 1;
printCentroids = 0; % a lot of console tables
plotAllFigures = 0;
compareExampleData = 0;
deprecated = 0; % unused code
% dirname = 'Archiwum';
% dirname = '23.05.23';
dirname = '33unique'; % 2 times author was recorded (34 patients number in sum)
% dirname = '/Users/puler/Documents/MATLAB/Ortheo3D/data'; % data via Qualisys from Delsys hardware
if(compareExampleData)
    dirname = '"../2przypadki/dane/Miopatia zapalna/Emg_Qemf.003/MVA_000 (1).wav"';
end
v = []; % Clear files data
files = dir(fullfile(dirname,'**','*.mat'));
datafiles = fullfile({files.folder},{files.name});

txPr = "Pośredni"; txPc = "Podchwyt";
if(lang == EN)
    txPr = "Intermediate grip"; txPc = "Supinated grip";
end
txBR = "Brachioradialis"; txBB = "Biceps brachii";
minActions = 10; skipedTrainingReppetinons = 0; % for segmentation selecting

k = 0;
for f = datafiles    
    k = k + 1;
    df = load(string(f));
%     fieldname = fieldnames(df);
    vars = fieldnames(df);
    for i = 1:length(vars)
        assignin('base', "tmp", df.(vars{i}));
        if(isfield(tmp,"Analog")) % Qualisys
            Qualisys = 1;
            tmpData = tmp.Analog.Data;
            tmpData = tmpData-mean(tmpData,2);
            tmpName = tmp.Analog.BoardName;
            if( tmpName == 'Delsys Trigno API')
                for(c = tmp.Analog.ChannelNumbers)
                    switch(tmp.Analog.Labels{c})
                        case 'L_Rectus Femoris' 
                            v(k).infoRDisp = txBR; 
                            v(k).dataR = tmpData(c,:); % Radialis
                            v(k).infoRecord = tmp.File(end-4);
                            v(k).infoRName = tmp.Analog.Labels{c};                                   
                        case 'L_Vastus Lateralis'
                            v(k).dataB = tmpData(c,:); % Biceps
                            v(k).infoRecord = tmp.File(end-4);              %to improve
                            v(k).infoBName = tmp.Analog.Labels{c};                 
                            v(k).infoBDisp = txBB; 
                    end           
                end
            else
                error("Unimplemented BoardName for Qualisys");
            end
        else % Noraxon Ultium
            Qualisys = 0;
            tmpData = tmp.movements.sources.signals.signal_1.data;
            tmpName = tmp.movements.sources.signals.signal_1.name;
            if ( tmpName == "Ultium EMG-BRACHIORAD. RT")
                v(k).dataR = tmpData; % Radialis
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoRName = tmpName;
                v(k).infoRDisp = txBR; 
            end
            if ( tmpName == "Ultium EMG-BICEPS BR. RT")
                v(k).dataB = tmpData; % Biceps
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoBName = tmpName;
                v(k).infoBDisp = txBB; 
            end
            
            tmpData = tmp.movements.sources.signals.signal_2.data;
            tmpName = tmp.movements.sources.signals.signal_2.name;
            if ( tmpName == "Ultium EMG-BRACHIORAD. RT")
                v(k).dataR = tmpData; % Radialis
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoRName = tmpName;
                v(k).infoRDisp = txBR; 
            end
            if ( tmpName == "Ultium EMG-BICEPS BR. RT")
                v(k).dataB = tmpData; % Biceps
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoBName = tmpName;
                v(k).infoBDisp = txBB; 
            end
        end % Qualisys or Ultium 
    end % vars
end % datafiles
fprintf(1,"Liczba pacjentów: %d\n", length(v)/2); toc(allElapsedTime); % / liczba mięśni mierzonych równocześnie lub wykonywanego ćwiczenia
fprintf(1,"Macierz plików v: %d (ilość pomiarów * liczba ćwiczeń)\n", length(v));

% get from last sample
if(Qualisys)
    Yunits = 'uV';
    fpom =  tmp.Analog.Frequency; 
    minActions = 2; 
else
    Yunits = tmp.movements.sources.signals.signal_1.units;
    fpom = tmp.movements.sources.signals.signal_1.frequency; 
end
dtpom=1/fpom; %  Hz

%------SEGMENTATION--------------------------------------------------------
if(exist("segments.mat"))
    load segments.mat
else
    nrF = 20; segmentActions;
end

if (skipedTrainingReppetinons) fprintf(1,"Rozkład: %d (pośrednich) %d (podchwytów) Sumarycznie pominiętych ćwiczeń: %d \n", countPosredni, countPodchwyt, skipedTrainingReppetinons); end;
fprintf(1,"Liczba segmentów: %d (liczba plików * liczba mięśni * ilość powtórzeń)\nSzacowanych segmentów: (%d)\n", length(segment), length(v)*2*10); % liczba plików * liczba mięśni * ilość powtórzeń
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
    if(1)%0<Parseval)
        sLmin = 1e10; % minimal segment length 
        for( i = 1:length(segment))
            if(sLmin > length(segment(i).data))
                sLmin = length(segment(i).data);
            end
        end
        sLmax = 0; % maximal segment length 
        for( i = 1:length(segment))
            if(sLmax < length(segment(i).data))
                sLmax = length(segment(i).data);
            end
        end
    end
    sLmax = lSyg;
    if (windowing) fprintf(1, "Okienkowanie Hanna "); end;
    for( i = 1:length(segment)) % uzupełnianie zerami segmentów
        SygRawLen(i) = length(segment(i).data');
        % if(Parseval)
        %     % from noParseval branch
        %     f = find(segment(i).data == max(segment(i).data));
        %     begIndx = f(1)-(l/2);
        %     if(begIndx<1) begIndx = 1; end
        %     endIndx = begIndx + l-1;
        %     if(endIndx > SygRawLen(i)) endIndx = SygRawLen(i); begIndx = endIndx-l+1; end
        %     tmp = segment(i).data(begIndx:endIndx);
        %     win = hann(l);
        %     if (windowing)           
        %         tmp = tmp.*win;
        %     end
        %     Syg(i,1:l) = [tmp];
        % else % noParseval            
            if (windowing)
                win = hann(SygRawLen(i));
                if(iscolumn(segment(i).data))
                    segment(i).data = segment(i).data.*win;
                else
                    segment(i).data = segment(i).data'.*win;
                end
            end           
        % end
        Syg(i,1:lSyg) = [segment(i).data' zeros(1, lSyg-length(segment(i).data))];
        SygKat(i) = segMio(i)+(segTraining(i)-1)*2; %plikSegMio(fileSegNr(nrs),nrB).i=2;  n+v(j).infoTraining-1*2; % training
        sygKat(i) = segment(i).miesien+(segment(i).gest-1)*2;
        segment(i).kat = sygKat(i);
        segment(i).len = SygRawLen(i); 
        % 1 - Pośred BR, 2 Poś BB, 3 - Podchwyt BRadialis, 4 Podch BBiceps          
    end    
    toc;
    save signals.mat Syg SygKat SygRawLen segment sLmax sLmin
end
%------SPECTRUMS-----------------------------------------------------------
nrF = 500; selectTraining; % wybór indeksów "j"
MTF(1).Tu = []; MTF(2).Tu = []; MTF(3).Tu = [];

if(exist("spectrums.mat"))
    load spectrums.mat
else
    nrF = 1000; fprintf(1, "Liczenie widm... "); spectrumTrend;
end
fprintf(1,"Rozmiar widm: %dx%d\n", size(Widma));

%------CENTROIDS-----------------------------------------------------------
if(exist("centroids.mat"))
    load centroids.mat
else
    nrF = 2000; fprintf(1, "Liczenie centroidów... ");centroid;
end
fprintf(1,"Rozmiar centroidów: %dx%d\n", size(CentrWidm)); 
fprintf(1,"Rozmiar centroidów uśrednionych: %dx%d\n", size(lpacj)); 

dCentr;
%------MINKOWSKI-DISTANCE--CentrWidm-wyglWidma/maxAf-..-CCE-wyglWidma------
nrF = 4000; 
for(jakieDist = 1:4)
    if(mod(jakieDist,2) == 1 ) nrF = nrF+1; bf = 0; else bf = 4; end
    if(jakieDist>2) flagaMaxima = 0; else flagaMaxima = 1; end 
    nrFig = jakieDist*2+i+500; 
    minkowskiDist;    
    figure(nrF);
    disppolt;
end

%------MINKOWSKI-DISTANCE--CC-CentrWidm--CCE-CentrWidm---------------------
for(jakieDist = 5:6)
    if(mod(jakieDist,2) == 1 ) nrF = nrF+1; bf = 0; figure(nrF); else bf = 4; end
    if(jakieDist>5) flagaMaxima = 0; else flagaMaxima = 1; end
    nrFig = jakieDist*2+500; 
    minkowskiDist;
    figure(nrF); 
    ddplot; 
end
fprintf(1, "main = "); toc(allElapsedTime);

nrF = 5e3; 
clear s;
  PL = 1; EN = 2;
  s.lang = EN;
  % s.figNrList = [ 1 22 33 34 ];
  % s.figPath = strcat("figBase/",string(k),"_",date,nameWithoutExt,"/");
  % s.exportPath = "article/fig/";
  % s.prefix = string(k);

  % h = strfind(folder,"/data/");
  % d = folder(h+6:end);
  % g = strfind(d,"/");
  % d = d(1:g-1);
  % s.filename = strcat("_", d, date, "_");
  % s.skipSavedFig = 1;
 s.figNrList = [    
     251
        253
        % 501
        % 502 % done by shringking
        503
        509
        1013
        1023
        2005
        2059
        2069
        2078
        4001
        4002
        4003 ];

if(forArticle)
    tic; save2folder(s); toc; 
    close all force
    if(forArticle) s.exportPath = "article/fi/"; s.inch3dot25 = -2; end
    save4article(s);
    fprintf(1, "main + forArticle(1) = "); toc(allElapsedTime);
end
return

%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mnoznik = 24;
BR = v(11).dataR; BB = v(11).dataB; % radialis, biceps
L = resample(BR,mnoznik,1); R = resample(BB,mnoznik,1);
L = L/max(abs(L)); R = R/max(abs(R));
% stereo_mtx = [L, R];%, v(11).dataB/max(abs(v(11).dataB))];
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

nrF = 4000;
figureNubers = [134 136];
for( i = 200:100:nrF) figureNubers = [figureNubers i+16 i+17]; end
tic;
figPSW
toc;

% iloczyn wektorywy tylko w przestrzeni euclidesa

return
kol='kr';
for(i=1:2)
    ys=segment(i).data'; sLmin = length(ys); tsym=sLmin*dtpom; tsym=1;
    W(i).Ays = abs(fft(ys))*2/sLmin; LAs=sLmin;
    ys0=[ys  zeros(1,lSyg*16-sLmin)]; l0=length(ys0);  tsym0=sLmin*dtpom; tsym0=1;
    W(i).Ays0=abs(fft(ys0))*2/l0;
    figure(55), subplot(1,2,1); plot([0:LAs-1]/tsym,W(i).Ays(1:LAs),[kol(i) '.-']), hold on;
    subplot(1,2,2),plot([0:16*LAs-1]/tsym0,W(i).Ays0(1:16*LAs),[kol(i) '.-']); hold on;
end
subplot(1,2,1); hold off; subplot(1,2,2); hold off; 

N=1000; xy=2*sin(2*pi/N*[1:5*N])-3*sin(4*pi/N*[1:5*N])+4*sin(6*pi/N*[1:5*N]); AaN=abs(fft(xy))/(2.5*N); A2aN=AaN.^2; Aa2N=abs(fft(xy.^2))/(2.5*N);
figure(33); plot([0:49]/(3),A2a(1:50),'k.-',[0:49]/(3),Aa2(1:50),'r.-',[0:49]/(2.5),A2aN(1:50),'k.--',[0:49]/(2.5),Aa2N(1:50),'r.--');
