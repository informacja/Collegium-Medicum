nrs = 1; %nr segmentu | licznik globalny
m = 0; % max length of segment 
countPodchwyt = 0;  countPosredni = 0;
names = [];
%--------------------------------------------------------------------------
% Segmentacja plików
%--------------------------------------------------------------------------
wybrane = [3 4 19 20];
j = 1;
while (j <= length(v))
    n1 = 1;  

    y = v(j).dataR; %a(j).d;%v(j).data;
    tmpData = v(j).dataR + v(j).dataB;
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
    name = [ v(j).infoRrecord v(j).infoRDisp ];
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
        figure(100+(j+1)*2), findpeaks(y,fpom,'MinPeakProminence',mp,'MinPeakDistance',md,'Threshold',15,"Annotate","peaks"); title(name);        
        v(j) = [];
        skipedTrainingReppetinons = skipedTrainingReppetinons+1;continue; 
    end % Skip not well segmented training
%     names = [names; string(v(j).infoRecord)];
    if (v(j).infoRecord(8) == 'c') countPodchwyt = countPodchwyt+1; end;
    if (v(j).infoRecord(8) == 'r') countPosredni = countPosredni+1; end;
%     figure(50+j), plot(peaksOfSignal); hold on;plot(TF*cutoffA*10); hold off; 
%     p2p = peak2peak(y);
    %mp = p2p/3; not recomended
%     md = length(y)/(numOfActions*1.9);

    figure(100+(j+1)*2), findpeaks(y,fpom, 'MinPeakProminence',mp,'MinPeakDistance',md,"Annotate","peaks",'Threshold',15); title(name);
%     clear pks locs,width,prominence
    [pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',md*fpom,"Annotate","peaks",'Threshold',15,'Annotate','extents');
    
%     figure(50+j), plot(y); hold on;plot(locs, pks); hold off; 

    locs(length(locs)+1) = length(y);
    for(k = 1:numOfActions)
        e = envelope(y(locs(k):locs(k+1)),1111,'rms'); %[n1,locs(k),locs(k+1), locs(k)+find(min(e)==e)]
%         figure(200+nrs), plot(e);
        Nbf = (locs(k)+find(min(e)==e)); % minimum punkt podziału
%         Nbf = Nbf(1);
        for(n = 0:1) % muscules signals const TODO 1
            if (n) s = v(j).dataR(n1:Nbf); fileSegMio(nrs) = txBR;
            else   s = v(j).dataB(n1:Nbf); fileSegMio(nrs) = txBB; end

    %         figure(100+nrs), plot (s);
            fileSegNr(nrs) = j; % w przypadku różnych ilości segmentów w plikach
            segment(nrs).data = s;
            m = max(m,length(s));
            %n1 = Nbf+1; %Nbf = length(v(j).dataR); [n1 Nbf]
            nrs=nrs+1; 
        end    
        n1 = Nbf+1;
    end
    j = j + 1;
end

% if liczba segmentó mniejsza niż 9 skip

save segments.mat segment fileSegNr fileSegMio m skipedTrainingReppetinons v countPodchwyt countPosredni

if(DEBUG) figure(100+(j+1)*2+1); hist(fileSegNr); title("Rozkład segmentów w plikach"); ylabel("Segmenty per ćwiczenie"); xlabel("Liczba ćwiczeń"); end; % ilość segmentów per plik