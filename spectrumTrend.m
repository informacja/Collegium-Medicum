tic; 
ksyg = 0;
lAfMin = 10e1000;
lAfMax = 0;
mnoznik = 2;

lingua = dictionary2(lang, Yunits);

if(Parseval>0)    mnoznik = .5; % z zerami
else              mnoznik = .15; % bez zer
end

for(j = 1:length(v)) % grupa
%     figure(nrF+j),
    nxf = 0; % nie drukuj figur
    fileSegNr; %TODO
   
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg))
        ifig = segMio(nseg(i))-1; % ifigure SygKat
%         ifig = SygKat(nseg(i))-1; % ifigure SygKat
        nj = find(wybrJ==j);
        ifLastSeg = (~isempty(nj) && (i == length(nseg) || i == length(nseg)-1));
        if( ifLastSeg || plotAllFigures) figure(j+nrF); end
       
        clear y;
        ksyg=nseg(i);
        if(Parseval>0)
            y = Syg(ksyg,:)'; %v(jdata(n1:Nbf);
        else
            y = Syg(ksyg, 1:segment(ksyg).len)';
        end
        X = (y.^2); Esyg(j,i)=sum(X)*dtpom;%/lSyg;         
        %segment(nrs).data = y;
        Nf=length(y); %todo
        nx=[0:Nf-1];
       if( ifLastSeg || plotAllFigures )
            subplot(lfrow,lc,1+ifig), plot(nx*dtpom, y); axis('tight');
            if( ifig == 0) title(lingua.t1); end
            subtitle(sprintf(lingua.subt1, fileSegMio(ksyg), ksyg));
            xlabel(sprintf(lingua.xl1, i));
            ylabel(lingua.yl1)
            subplot(lfrow,lc,1+ifig+lc),  plot(nx*dtpom, X);axis('tight');
        end 
        % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Twygl=0.25*mnoznik; nTu = Tsyg/Twygl;
        Tu=Twygl/dtpom; nfw = 1;
        run("../MTF/filtrWidma.m");
        xf = [0:LwAm-1];
        if( ifLastSeg || plotAllFigures )
            hold on; plot(xf*dtpom,Ayf,'Color','#0072BD'); % ciemy niebieski - domyślny
            plot([0:Ldf]*dtpom,Af,'k'); axis('tight');  hold off; axis('tight');
            xlabel(sprintf(lingua.xl2, i,Tu*dtpom*1000));
        end
        
        s.y = y; s.l = SygRawLen(ksyg); s.Nf = Nf; s.mnoznik = mnoznik; s.Parseval = Parseval; s.Tsyg = Tsyg; s.sL = Nf; s.dtpom = dtpom; s.MTF = MTF; s.nxf = nxf; s.sLmax = sLmax; 
        s.doKwadratu = 0;
        [Ayf, Af, ~, nk, dx, ix, x] = st(s);

        Widma(j,i).Ayf=Ayf;     % raw spectrum
        wyglWidma(j,i).Af=Af;   % interpolowane
        
        % figure(1),   plot([0:length(Ayf)-1],Ayf,'c',[0:dx:(length(Af)-1)*dx],Af,'k.-',[0:length(Afor)-1],Afor,'r*'); axis('tight');length(nseg) % interpolata
        % nf=round(Nf/Podzial);
        % kf=SygKat(nseg(i));

        Ldf = length(Af)-1; % by PSW 27.02.24
        if( ifLastSeg || plotAllFigures )
            subplot(lfrow,lc,1+ifig+2*lc),   plot([0:nk-1]/Tsyg,Ayf(1:nk),'c',x,Af(1:ix),'k'); axis('tight'); %
            % figure(nrFw), subplot(1,2,1); hold on; 
             %kf=mod(kf,4)+1; 
            % plot([0:dx:(length(Af)-1)*dx],Af,kol(kf)); axis('tight'); %plot(wyglWidma(j,i).Af); hold off; 
            figure(nrF+j)
            if( 1+ifig+2*lc == 5 ) title(lingua.t2); end
            xlabel(sprintf(lingua.xl3,i,Tu*dtpom*1000,1/(Tu*dtpom)));
        end
        
        s.doKwadratu = 1;
        [Ayf, Af, Afw, nf, dx, ix, x] = st(s);

        Widma(j,i).Ayf2=Ayf; wyglWidma(j,i).Af2=Af;
        % nf=round(Nf/Podzial);

        if( ifLastSeg || plotAllFigures )
%             figure(nrFw+j),
            subplot(lfrow,lc,1+ifig+3*lc), plot([0:nf-1]/Tsyg,Ayf(1:nf),'c',x,Af(1:ix),'k');
            xlabel(sprintf(lingua.xl4,i,Tu*dtpom*1000));
            % sgtitle( sprintf("%s",v(j).infoBDisp)); 
            axis('tight');
            % figure(nrFw), subplot(1,2,2); hold on; plot([0:nf-1]/Tsyg,Af(1:nf),kol(kf)); %plot(wyglWidma(j,i).Af); hold off; 
            maxAf(j)=max(Af); 
            sgtitle(lingua.sgt); axis('tight');
        end
        lAfMin = min(length(Af),lAfMin);
        lAfMax = max(length(Af),lAfMax);
    end
end
hold off;

toc;
save spectrums.mat Widma wyglWidma Esyg lAfMin lAfMax

function [Ayf Af Afw nk dx ix x] = st(s) % spectrum&trend
    A = fft(s.y); lA = length(A);
    Afw = abs(A(1:round(lA/2)))/s.l; %lA; 
    Podzial=4; if(j==2) Podzial=10; end
    nf=round(s.Nf/2); 
    X = Afw(1:nf); podstawa = 0.05;
    if(s.doKwadratu) X=Afw(1:nf).^2;  podstawa = 0.025; end
    Podzial=15; if(j==2) Podzial=30; end
    Twygl=podstawa*s.mnoznik; % jednakowe dla wszystkich pacjentów 
    
    dx=s.sL/s.sLmax;
    if(s.Parseval<0) Twygl=Twygl*dx;  end % wzgledne
   
    nTu = s.Tsyg/Twygl; % LSyg / nTu jest liczbą próbek w oknie wygładzania
    Tu=Twygl/s.dtpom; nfw = 2;
    MTF = s.MTF;
    nxf = s.nxf;
    run("../MTF/filtrWidma.m");
    nk=round(s.sL/Podzial); 
    snL=length(Af);
    
    nf=round(s.Nf/Podzial);
    if(s.Parseval<0) % interpolacja
        x=1:dx:snL; length(x)-s.sLmax;
        Afor=Af; % wtgł nie interpolowane
        Af=csapi(1:snL,Afor,x);
    
        ix = round(nf*(1/dx));% interpolowana skala x
        if(snL<ix) ix = snL; end
        x = [0:dx:(ix-1)*dx]/s.Tsyg;
    else
        ix = nf;
        x = [0:(ix-1)]/s.Tsyg;
    end
end

function [d] = dictionary2(lang,Yunits)
    PL = 1;
    EN = 2;
    dict(PL).t1 = "                                                                                               Dziedzina czasu";
    dict(PL).subt1 = "Mięsień: %s, segment nr: %d ";
    dict(PL).xl1 = "Chwyt pośredni nr (względny) %d : y(t) t[sek]";
    dict(PL).yl1 = ['Amplituda [' Yunits ']'];
    dict(PL).xl2 = "Kwadrat y i wygł. y(t)^2 %d (Tu=%.1fms): t[sek]";
    dict(PL).t2 = "                                                                                            Dziedzina częstotliwości";
    dict(PL).xl3 = "Widmo %d [Hz] Tu=%.1fms f_g=1/Tu=%.0fHz ";
    dict(PL).xl4 = "Widmo mocy %d f_g=1/Tu Tu=%.1fms";
    dict(PL).sgt = "Widma sygnałów wewn. grupy dla raw i mocy";

    dict(EN).t1 = "                                                                                                    Time domain";
    dict(EN).subt1 = "Muscle: %s, segment nr: %d ";
    dict(EN).xl1 = "Intermediate grip nr (relative) %d : y(t) t[s]";
    dict(EN).yl1 = ['Amplitude [' Yunits ']'];
    dict(EN).xl2 = "Square y and smooth y(t)^2 %d (Tu=%.1fms): t[s]";
    dict(EN).t2 = "                                                                                            Frequency domain";
    dict(EN).xl3 = "Spectrum %d [Hz] Tu=%.1fms f_g=1/Tu=%.0fHz ";
    dict(EN).xl4 = "Spectrum of power %d f_g=1/Tu Tu=%.1fms";
    dict(EN).sgt = "Spectrum for intra-group raw & power signals";

    d = dict(lang);
end

% disp("Normowanie widma")

% normowanie widma - depracated
% for(j = 1:length(v)) % grupa
%     lAyf=length(Widma(j,1).Ayf);
%     CentrWidm(j,1).Ayf=zeros(lAyf,1); nAyf(1)=0; nAyf(2)=0; 
%     CentrWidm(j,1).Ayf2=zeros(lAyf,1);
%     CentrWidm(j,2).Ayf=zeros(lAyf,1); nAyf(1)=0; nAyf(2)=0; 
%     CentrWidm(j,2).Ayf2=zeros(lAyf,1);
%     nseg=find(fileSegNr==j);
%     for (i = 1:length(nseg)) 
%         Widma(j,i).maxAyf = max(Widma(j,i).Ayf);
%         Widma(j,i).maxAyf2 = max(Widma(j,i).Ayf2);
%         kat = segment(nseg(i)).miesien;
%         nAyf(kat)=nAyf(kat)+1;
%         CentrWidm(j, kat).Ayf=CentrWidm(j, kat).Ayf+Widma(j,i).Ayf/Widma(j,i).maxAyf; 
%         CentrWidm(j, kat).Ayf2=CentrWidm(j, kat).Ayf2+Widma(j,i).Ayf2/Widma(j,i).maxAyf2; 
%     end
%     CentrWidm(j, 1).Ayf=CentrWidm(j,1).Ayf/nAyf(1); 
%     CentrWidm(j, 2).Ayf=CentrWidm(j,2).Ayf/nAyf(2);
%     CentrWidm(j, 1).Ayf2=CentrWidm(j,1).Ayf2/nAyf(1); 
%     CentrWidm(j, 2).Ayf2=CentrWidm(j,2).Ayf2/nAyf(2); 
% end
