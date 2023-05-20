%biurko, matlab ten skrypt
dirname = 'Archiwum';
clear v
v.data = [];
v.info = [];
files = dir(fullfile(dirname,'**','*.mat'));
datafiles = fullfile({files.folder},{files.name});
fpom=2000; dtpom=1/fpom; % Hz
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

% posredni = record_2023_04_20_15_06_BR_po_redni.movements.sources.signals.signal_1.data;
% podchwyt = record_2023_04_20_15_07_BR_podchwyt.movements.sources.signals.signal_1.data;
lfrow= 4; lc=2; n1 = 18000; nrs = 1; %nr segmentu
m = 0;
for (j = 1:length(v))
    n1 = 1;
    y = v(j).data;
    %         run("../MTF/filtrWidma.m");
    %         figure, plot(Af);
    p2p = peak2peak(y);
    mp = p2p/2;
    md = length(y)/3;
    %     figure((j+1)*2), findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',md,"Annotate","peaks")
    [pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',md,"Annotate","peaks",'Annotate','extents');

    e = envelope(y(locs(1):locs(2)),1111,'rms');
    %    figure(55), plot(e);
    Nbf = (locs(1)+find(min(e)==e));

    for( i = 0:1)
        y = v(j).data(n1:Nbf); X = y.^2;
        segment(nrs).data = y;
        m = max(m,length(y));
        n1 = Nbf+1; Nbf = length(v(j).data);
        nrs=nrs+1;
    end
end
Tsyg=ceil(ceil(m*dtpom)/2)*2; % [sek]
lSyg=round(Tsyg/dtpom);
clear Syg;
ksyg=0;
for(j = 1:length(v))
    figure(j), nxf = 0; % nie drukuj figur
    for( i = 1:length(segment))
        Syg(i,1:lSyg) = [segment(i).data' zeros(1,lSyg-length(segment(i).data))];
    end
    figure(j)
    for (i = 0:1)
        clear y;
        ksyg=ksyg+1;
        y = Syg(ksyg,:)'; %v(j).data(n1:Nbf);
        X = y.^2;
        %segment(nrs).data = y;
        Nf=length(y); %todo
        nx=[0:Nf-1];
        title("Dziedzina czasu");
        subplot(lfrow,lc,1+i),  plot(nx*dtpom, y); xlabel(sprintf("Ruch pośredni %d: y(t) t[sek]", i));
        subplot(lfrow,lc,1+i+lc),  plot(nx*dtpom, X);
        % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Twygl=0.25; nTu = Tsyg/Twygl;
        Tu=Twygl/dtpom;
        run("../MTF/filtrWidma.m");
        xf = [0:LwAm-1];
        hold on; plot(xf/Tsyg,Ayf,'c',[0:Ldf]/Tsyg,Af,'k'); axis('tight');  hold off;
        xlabel(sprintf("Kwadrat y i wygł.y(t)^2 %d (Tu=%.1fms):     t[sek]", i,Tu*dtpom*1000));
        A = fft(y); lA = length(A);
        Afw = abs(A(1:round(lA/2)))/lA;
        Podzial=4; if(j==2) Podzial=10; end
        nf=round(Nf/2); 
        X = Afw(1:nf);
        Twygl=0.05; nTu = Tsyg/Twygl; % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Tu=Twygl/dtpom;
        run("../MTF/filtrWidma.m");
        nk=round(Nf/Podzial);
        subplot(lfrow,lc,1+i+2*lc),   plot([0:nk-1]/Tsyg,Ayf(1:nk),'c',[0:Ldf]/Tsyg,Af,'k');
        Widma(j,i+1).Ayf=Ayf; wyglWidma(j,i+1).Af=Af;% i*2+j
        title("Dziedzina częstotliwości");
        xlabel(sprintf("Widmo %d [Hz] Tu=%.1fms f_g=1/Tu=%.0fHz ",i,Tu*dtpom*1000,1/(Tu*dtpom)));
        Podzial=15; if(j==2) Podzial=30; end
        Twygl=0.025; nTu = Tsyg/Twygl;
        nf=round(Nf/Podzial);
        X=Afw(1:nf).^2;
        Tu=Twygl/dtpom;
        run("../MTF/filtrWidma.m");
        subplot(lfrow,lc,1+i+3*lc), plot([0:LwAm-1]/Tsyg,Ayf,'c',[0:Ldf]/Tsyg,Af,'k');
        xlabel(sprintf("Widmo mocy %d f_g=1/Tu Tu=%.1fms",i,Tu*dtpom*1000));
        Widma(j,i+1).Ayf2=Ayf; wyglWidma(j,i+1).Af2=Af;% i*2+j

        sgtitle(v(j).info)
        %nrs = nrs +1;
    end
end
for (j = 1:length(v))
    dW=Widma(j,1).Ayf-Widma(j,2).Ayf; S(j,1)=sqrt(dW'*dW);
    dWwygl=wyglWidma(j,1).Af-wyglWidma(j,2).Af; S(j,2)=sqrt(dWwygl*dWwygl');
    %sredne widma
    Wsr(j).Wsr=(Widma(j,1).Ayf'+Widma(j,2).Ayf')/2; %S(j,1)=sqrt(dW'*dW); % wzorzec
    Wsrwygl(j).Wsr=(wyglWidma(j,1).Af+wyglWidma(j,2).Af)/2; %S(j,2)=sqrt(dWwygl*dWwygl');

    % D2
    dW2=Widma(j,1).Ayf2-Widma(j,2).Ayf2; S2(j,1)=sqrt(dW2'*dW2);
    dW2wygl=wyglWidma(j,1).Af2-wyglWidma(j,2).Af2; S2(j,2)=sqrt(dW2wygl*dW2wygl');
end
return
DWsr=Wsr(1).Wsr-Wsr(2,:); 
Ssr(1)=sqrt(DWsr'*DWsr);
DWsr=Wsrwygl(1,:)-Wsrwygl(2,:); Ssr(2)=sqrt(DWsr'*DWsr);

 figure, plot(1:length(Wsr), Wsr,'r', 1:length(Widma(1,1).Ayf),Widma(1,1).Ayf,'c', 1:length(Widma(1,2).Ayf),Widma(1,2).Ayf,'k')% czerwone srednie widmo

Widmo = S
widmoMocy = sqrt(S2)
porownywalnoscIrozroznialnosc=[S(1,:)./S(2,:);S(:,1)'./S(:,2)']
% po wierszach duża różnica po kolumnach mała
Ssr
% figPSW jedakowe długości, wzorzec musi mieć takie same