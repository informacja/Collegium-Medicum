tic; 
ksyg = 0;
for(j = 1:length(v)) % grupa
%     figure(nrF+j),
    nxf = 0; % nie drukuj figur
    fileSegNr; %TODO
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg))
        ifig = segMio(nseg(i))-1; % ifigure SygKat
%         ifig = SygKat(nseg(i))-1; % ifigure SygKat
       
        if(j == length(v))
            figure(nrF+j)
        end
        clear y;
        ksyg=ksyg+1;
        y = Syg(ksyg,:)'; %v(jdata(n1:Nbf);
        X = (y.^2); Esyg(j,i)=sum(X)*dtpom;%/lSyg;         
        %segment(nrs).data = y;
        Nf=length(y); %todo
        nx=[0:Nf-1];
        if j == length(v)
            subplot(lfrow,lc,1+ifig),  plot(nx*dtpom, y); xlabel(sprintf("Ruch pośredni %d: y(t) t[sek]", i));
            if( ifig == 0) title("                                                                                                                                   Dziedzina czasu"); end
            ylabel(['Amplituda [' Yunits ']'])
            subplot(lfrow,lc,1+ifig+lc),  plot(nx*dtpom, X);
        end 
        % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Twygl=0.25; nTu = Tsyg/Twygl;
        Tu=Twygl/dtpom; nfw = 1;
        run("../MTF/filtrWidma.m");
        xf = [0:LwAm-1];
        if j == length(v)
            hold on; plot(xf/Tsyg,Ayf,'c',[0:Ldf]/Tsyg,Af,'k'); axis('tight');  hold off;
            xlabel(sprintf("Kwadrat y i wygł. y(t)^2 %d (Tu=%.1fms): t[sek]", i,Tu*dtpom*1000));
        end
        A = fft(y); lA = length(A);
        Afw = abs(A(1:round(lA/2)))/SygRawLen(ksyg); %lA; 
        Podzial=4; if(j==2) Podzial=10; end
        nf=round(Nf/2); 
        X = Afw(1:nf);
        Twygl=0.05; nTu = Tsyg/Twygl; % LSyg / nTu jest liczbą próbek w oknie wygładzania
        Tu=Twygl/dtpom; nfw = 2;
        run("../MTF/filtrWidma.m");
        nk=round(Nf/Podzial);
        Widma(j,i).Ayf=Ayf; wyglWidma(j,i).Af=Af;% i*2+j
        v(j).kat = n+v(j).infoTraining-1*2; % training
        if j == length(v)
            subplot(lfrow,lc,1+ifig+2*lc),   plot([0:nk-1]/Tsyg,Ayf(1:nk),'c',[0:Ldf]/Tsyg,Af,'k');
            figure(nrFw), subplot(1,2,1); hold on; kf=mod(kf,4)+1; plot([0:Ldf],Af,kol(kf)); %plot(wyglWidma(j,i).Af); hold off; 
%         figPW("png")
            figure(nrF+j)
            if( 1+ifig+2*lc == 5 ) title("                                                                                                                                 Dziedzina częstotliwości"); end
            xlabel(sprintf("Widmo %d [Hz] Tu=%.1fms f_g=1/Tu=%.0fHz ",i,Tu*dtpom*1000,1/(Tu*dtpom)));
        end
        Podzial=15; if(j==2) Podzial=30; end
        Twygl=0.025; nTu = Tsyg/Twygl;
%         nf=Nf;%round(Nf/Podzial);
%         Xx=fft(y.^2);Xx=abs(Xx(1:nf));  
        X=Afw(1:nf).^2;
%         figure(111), plot(X); hold on; plot(Xx)
        Tu=Twygl/dtpom; nfw = 3;
        run("../MTF/filtrWidma.m");
        nf=round(Nf/Podzial);
        Widma(j,i).Ayf2=Ayf; wyglWidma(j,i).Af2=Af;% i*2+j

%         Widma(j,i).kat = v(j).kat = n+v(j).infoTraining-1*2; % training
%         to erase
%         figPW("png",5)
        if j == length(v)
            subplot(lfrow,lc,1+ifig+3*lc), plot([0:nf-1]/Tsyg,Ayf(1:nf),'c',[0:nf-1]/Tsyg,Af(1:nf),'k');
            xlabel(sprintf("Widmo mocy %d f_g=1/Tu Tu=%.1fms",i,Tu*dtpom*1000));
            sgtitle( sprintf("%s %d",v(j).infoBDisp, ksyg))
            figure(nrFw), subplot(1,2,2); hold on; plot([0:nf-1]/Tsyg,Af(1:nf),kol(kf)); %plot(wyglWidma(j,i).Af); hold off; 
            maxAf(j)=max(Af); 
            sgtitle("Widma sygnałów wewn. grupy dla raw i mocy")
        end
        %nrs = nrs +1;
    end
 
end

% normowanie widma
for(j = 1:length(v)) % grupa
    lAyf=length(Widma(j,1).Ayf);
    CentrWidm(j,1).Ayf=zeros(lAyf,1); nAyf(1)=0; nAyf(2)=0; 
    CentrWidm(j,1).Ayf2=zeros(lAyf,1);
    CentrWidm(j,2).Ayf=zeros(lAyf,1); nAyf(1)=0; nAyf(2)=0; 
    CentrWidm(j,2).Ayf2=zeros(lAyf,1);
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg)) 
        Widma(j,i).maxAyf = max(Widma(j,i).Ayf);
        Widma(j,i).maxAyf2 = max(Widma(j,i).Ayf2);
        if (segMio(nseg(i)) == 1 ) kat = 1; else kat = 2; end
        nAyf(kat)=nAyf(kat)+1;
        CentrWidm(j, kat).Ayf=CentrWidm(j, kat).Ayf+Widma(j,i).Ayf/Widma(j,i).maxAyf; 
        CentrWidm(j, kat).Ayf2=CentrWidm(j, kat).Ayf2+Widma(j,i).Ayf2/Widma(j,i).maxAyf2; 
    end
    CentrWidm(j, 1).Ayf=CentrWidm(j,1).Ayf/nAyf(1); 
    CentrWidm(j, 2).Ayf=CentrWidm(j,2).Ayf/nAyf(2);
    CentrWidm(j, 1).Ayf2=CentrWidm(j,1).Ayf2/nAyf(1); 
    CentrWidm(j, 2).Ayf2=CentrWidm(j,2).Ayf2/nAyf(2); 
end
toc;
% figure(nrFw+1); subplot(1,2,1) 
save spectrums.mat Widma wyglWidma Esyg CentrWidm
