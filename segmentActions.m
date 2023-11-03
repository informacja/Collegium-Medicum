tic; nrs = 1; %nr segmentu | licznik globalny
m = 0; % max length of segment 
countPodchwyt = 0;  countPosredni = 0;
names = [];
%--------------------------------------------------------------------------
% Segmentacja plików
%--------------------------------------------------------------------------
% wybrane = [3 4 19 20];
selectTraining
j = 1; 
while (j <= length(v))
    nrR = 0; nrB = 0;
    n1 = 1;  

    y = v(j).dataR; %a(j).d;%v(j).data;
    tmpData = v(j).dataR + v(j).dataB; % założenie o byciu w fazie
    y = tmpData;
%     Nokno = 12;
%     lY = length(y);
%     if( 0 && j == 8 ) 
%         s = zeros(1,lY); s(1)=y(1);
%         for (n = 2:lY) %Nokno:lY)
%             s(n)=s(n-1)+y(n); %s(n-Nokno+Nokno/2) = mean(y(n-Nokno+1:n));
%         end
%         ys=y; y=s;
%         
%         figure(100), plot(ys,'g'); hold on; plot(y,'r');
%         y = s;
%     end
    if ~DEBUG
        if isempty(y) continue; end;
    end
    name = [ v(j).infoRecord v(j).infoRDisp ]; 
    %         run("../MTF/filtrWidma.m");
    %         figure, plot(Af);
    md = 3.9; %[ seconds] 3.8 was
%     md = length(y)/(numOfActions*1.9);
    minAmplitude = 500;%peak2peak(y)/10;
    mp = minAmplitude; %[uV]
    peaksOfSignal = findpeaks(y,fpom,'MinPeakDistance',md,'MinPeakProminence',mp,'Threshold',15);
%     findpeaks(y,fpom,'Annotate','extents','WidthReference','halfheight')
%     MinPeakWidth
% title('Signal Peak Widths')
    cutoffA = max(peaksOfSignal)/5;
    for( i = 1:length(peaksOfSignal) )
        if(peaksOfSignal(i) <= cutoffA)
            peaksOfSignal(i) = 0;
        else
            peaksOfSignal(i) = 1;
        end
    end
%     TF = islocalmax(peaksOfSignal,'SamplePoints',y); 
% TF = islocalmax(A,'MinSeparation',minutes(45),'SamplePoints',t);
    numOfActions = sum(peaksOfSignal);
    if numOfActions < minActions %||  numOfActions > 12 %|| j > 2
        figure(100+(j)*2), findpeaks(y,fpom,'MinPeakProminence',mp,'MinPeakDistance',md,'Threshold',15,"Annotate","peaks"); title(name);  xlabel("Czas [s]"); ylabel("Amplituda [uV]");     
        v(j) = []; % usuń record z macierzy wejściowej
        skipedTrainingReppetinons = skipedTrainingReppetinons+1; 
        continue; 
    end % Skip not well segmented training
%     names = [names; string(v(j).infoRecord)];
    if (length(find(v(j).infoRecord=='r'))) countPosredni = countPosredni+1; v(j).infoTraining = 1; end;
    if (length(find(v(j).infoRecord=='c'))) countPodchwyt = countPodchwyt+1; v(j).infoTraining = 2; end;
%     figure(50+j), plot(peaksOfSignal); hold on;plot(TF*cutoffA*10); hold off; 
%     p2p = peak2peak(y);
    %mp = p2p/3; not recomended
%     md = length(y)/(numOfActions*1.9);

% if (v(j).infoRecord == """11 pośredni""") || (v(j).infoRecord == """11 podchwyt""") 
%     figure(100+(j+1)*2), findpeaks(y,fpom, 'MinPeakProminence',mp,'MinPeakDistance',md,"Annotate","peaks",'Threshold',15,'Annotate','extents'); title(name);
% figPW('nomargin');
% end
    %     clear pks locs,width,prominence
    [pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',md*fpom,"Annotate","peaks",'Threshold',15,'Annotate','extents');
    
%     figure(50+j), plot(y); hold on;plot(locs, pks); hold off; 

    locs(length(locs)+1) = length(y);
    for(k = 1:numOfActions)
        v(j).segLen = numOfActions;
        e = envelope(y(locs(k):locs(k+1)),1111,'rms'); %[n1,locs(k),locs(k+1), locs(k)+find(min(e)==e)]
        
        if wybrJ(1) == j && DEBUG
            figure("Name", "RMS Envelope")
            envelope(y(locs(k):locs(k+1)),1111,'rms'); hold on; axis('tight'); xlabel("Nr próbki (względny)"); ylabel("Amplituda")
            n = find(min(e)==e);
            ax = axis; plot([n n], ax(3:4), 'k--');
%             ax = axis; plot([Nbf Nbf], ax(3:4), 'k--');
        end

%         figure(200+nrs), plot(e);
        Nbf = (locs(k)+find(min(e)==e)); % minimum punkt podziału
        if Nbf > length(v(j).dataB) 
            Nbf = length(v(j).dataB); paraSegmentowOgraniczonych = [nrs nrs+1]
%             continue;
        end
%         Nbf = Nbf(1);
        for(n = 1:2) % muscules signals const TODO 1
            if (n==1) s = v(j).dataR(n1:Nbf); fileSegMio(nrs) = txBR; segMio(nrs) = 1; nrR=nrR+1; segTraining(nrs) = v(j).infoTraining; plikSegMio(j,nrR).i=1;
%             plikSegMio(j,nrR). = nrs dla segmentu
            else   s = v(j).dataB(n1:Nbf); fileSegMio(nrs) = txBB; segMio(nrs) = 2; nrB=nrB+1; segTraining(nrs) = v(j).infoTraining; plikSegMio(j,nrB).i=2; end
    
    %         figure(100+nrs), plot (s);
            fileSegNr(nrs) = j; % w przypadku różnych ilości segmentów w plikach
            segment(nrs).data = s;
            segment(nrs).miesien = n;
            segment(nrs).gest = v(j).infoTraining;
%             v(j,k).klasyfikacjaSegmentów = nrs;
            m = max(m, length(s));
            %n1 = Nbf+1; %Nbf = length(v(j).dataR); [n1 Nbf]
            nrs=nrs+1; 
        end    
        n1 = Nbf+1;
    end
    j = j + 1;
end

% if liczba segmentó mniejsza niż 9 skip

save segments.mat segment fileSegNr fileSegMio m skipedTrainingReppetinons v countPodchwyt countPosredni segTraining segMio

if(DEBUG) 
%      figure,histogram(fileSegNr, length(v), 'BinWidth',1); title("Rozkład segmentów w plikach"); ylabel("Segmenty per ćwiczenie"); xlabel("Liczba ćwiczeń");  % ilość segmentów per plik
    sLens = []; for i = 1:length(segment) sLens(i) = length(segment(i).data);end
    figure,
    subplot(211),mhistogr(fileSegNr, length(v),1); axis("tight"); title("Rozkład segmentów w plikach"); ylabel("Segmenty per ćwiczenie"); xlabel("Liczba ćwiczeń");
    subplot(212),plot(sLens/fpom,'.'); hold on; axis("tight"); title("Rozkład długości segmentów"); ylabel("Długość segmentu [s]"); xlabel("Nr segmentu");
    axis auto 
%     figPW("nomargin") TODO?
    if(deprecated) 
        for i = 554:555%length(segment) % color anomaly
            sLens(i) = length(segment(i).data);
            e = envelope(segment(i).data, 9999,'rms');
            length(e)/length(segment(i).data)
            TF = islocalmin(e,'MinProminence', peak2peak(e)/1000);
    %         sum(TF)
            if(sum(TF) > 1)
    %             figure, plot(e); hold on; plot(segment(i).data); plot(e(TF), segment(i).data(TF)/fpom,'r.')
    %             figure, plot(segment(i).data); hold on; plot(e); 
            end
        end
    %     hold off; axis('tight')
    end
end
  toc;