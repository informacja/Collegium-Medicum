nrs = 1; %nr segmentu | licznik globalny
m = 0; % max length of segment 
%--------------------------------------------------------------------------
% Segmentacja plików
%--------------------------------------------------------------------------
for (j = 1:length(v))
    n1 = 1;  

    y = v(j).dataR; %a(j).d;%v(j).data;
    tmpData = v(j).dataR + v(j).dataB;
    y = tmpData;
    if ~DEBUG
        if isempty(y) continue; end;
    end
    name = [ v(j).infoRrecord v(j).infoRDisp ];
    %         run("../MTF/filtrWidma.m");
    %         figure, plot(Af);
    md = fpom; %[1 second]
    mp = 200; %[uV]
    peaksOfSignal = findpeaks(y,'MinPeakDistance',md,'MinPeakProminence',mp);
    cutoffA = max(peaksOfSignal)/10;
    for( i = 1:length(peaksOfSignal) )
        if(peaksOfSignal(i) <= cutoffA)
            peaksOfSignal(i) = 0;
        end
    end
    TF = islocalmax(peaksOfSignal);
    numOfActions = sum(TF);
%     figure(50+j), plot(peaksOfSignal); hold on;plot(TF*cutoffA*10); hold off; 
    p2p = peak2peak(y);
    %mp = p2p/3; not recomended
    md = length(y)/(numOfActions*1.4);
    figure((j+1)*2), findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',md,"Annotate","peaks"); title(name);
%     clear pks locs,width,prominence
    [pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',md,"Annotate","peaks",'Annotate','extents');
    
%     figure(50+j), plot(y); hold on;plot(locs, pks); hold off; 

    locs(length(locs)+1) = length(tmpData);
    for(k = 1:numOfActions)
        e = envelope(y(locs(k):locs(k+1)),1111,'rms'); %[n1,locs(k),locs(k+1), locs(k)+find(min(e)==e)]
%         figure(200+nrs), plot(e);
        Nbf = (locs(k)+find(min(e)==e)); % minimum punkt podziału
        for(n = 0:1) % muscules signals const
            if (n) s = v(j).dataR(n1:Nbf); 
            else   s = v(j).dataB(n1:Nbf); end

    %         figure(100+nrs), plot (s);
            fileSegNr(nrs) = j; % w przypadku różnych ilości segmentów w plikach
            segment(nrs).data = s;
            m = max(m,length(s));
            n1 = Nbf+1; %Nbf = length(v(j).data);
            nrs=nrs+1; 
        end             
    end
end

save segment.mat segment fileSegNr m

if(DEBUG) hist(fileSegNr); title("Ilość segmentów w plikach"); ylabel("Segmenty"); xlabel("Nr pliku"); end; % ilość segmentów per plik