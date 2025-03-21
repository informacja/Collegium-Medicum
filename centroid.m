tic;
% if(size(CentrWidm)==[2, length(Widma)]) depracated
%     CentrWidm = CentrWidm';
% end
lingua = dictionary2(lang);

figure(nrF);
tmpTxt = "";
for(j = 1:2) % grupa
    tmpTxt = [tmpTxt;v(j).infoRecord];
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg))
        nrs = nseg(i);    nsp=SygKat(nrs);
        %         if(fileSegMio(nrs) == txBR)
        L = length(segment(nrs).data);
        subplot(2,2,nsp), plot([1:L]/fpom, segment(nrs).data); hold on; title(nsp); xlabel(lingua.xl1)

        title(sprintf(lingua.t1, v(j).infoRecord))
        %     end
        %         if(fileSegMio(nrs) == txBB) subplot(2,2,ns), plot(segment(nrs).data); hold on; title(txBB);end
        %         if(fileSegMio(nrs) == txBR) subplot(2,2,3), plot(segment(nrs).data); hold on; end
        %         if(fileSegMio(nrs) == txBB) subplot(2,2,4), plot(segment(nrs).data); hold on; end
    end
    subplot(2,2,1), hold off;
    subplot(2,2,2), hold off;
    subplot(2,2,3), hold off;
    subplot(2,2,4), hold off;
end
% sgtitle(tmpTxt);

% f = [1:8500+1]; f = [1:6200+1]; f = 1:(length(wyglWidma(j,1).Af)); % max
% f = 1:155*max(Tsyg)+1; % [Hz]*[samples]
f = 1:fWyswieltCentroidow*max(Tsyg)+1; % [Hz]*[samples]
xf = (f-1)/Tsyg;     % [Hz]

clear dCentrM dCentrE;
nf=2;
nrs = 0;
nrF = nrF + length(v);
% normowanie widma
for(j = 1:length(v)) % grupa
    %     figure(nrFw+j);close(nrFw+j);
    %     figure(nrF+j);
    % normowanie widma
    % lAf=length(wyglWidma(j,1).Af);
    lAf = lAfMax;
    %     CentrWidm(j, 1).AfM=zeros(lAf,1); nAf(1)=0;
    CentrWidm(j, 1).AfM=zeros(lAf,1); nAf(1)=0;
    CentrWidm(j, 1).Af2M=zeros(lAf,1);
    CentrWidm(j, 1).AfE=zeros(lAf,1);
    CentrWidm(j, 1).Af2E=zeros(lAf,1);
    CentrWidm(j, 2).AfM=zeros(lAf,1); nAf(2)=0;
    CentrWidm(j, 2).Af2M=zeros(lAf,1);
    CentrWidm(j, 2).AfE=zeros(lAf,1);
    CentrWidm(j, 2).Af2E=zeros(lAf,1);
    %     Centr(j, 1).AfM = zeros(lAf,1);
    %     Centr(j, 1).AfE = zeros(lAf,1);
    %     Centr(j1 2).AfM = zeros(lAf,1);
    %     Centr(j, 2).AfE = zeros(lAf,1);
    nj = find(wybrJ==j);
    if( length(nj) > 0 || plotAllFigures) figure(j+nrF); end
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg)) 
        nrs = nseg(i);

        %         if(v(j).plikeSegMio(i) == txBR) kat = 1; end

        if(fileSegMio(nrs) == txBR) nrskat(1) = nrs; kat = 1; end
        if(fileSegMio(nrs) == txBB) nrskat(2) = nrs; kat = 2; end
        wyglWidma(j,i).maxAf = max(wyglWidma(j,i).Af);
        wyglWidma(j,i).maxAf2 = max(wyglWidma(j,i).Af2);
        %         if (wyglWidma() todo
        nAf(kat)=nAf(kat)+1;

        Psyg(j,i)=Esyg(j,i)/SygRawLen(nrs);
        Psr(j) = mean(Psyg(j,:)); % Power
        Esr(j) = mean(Esyg(j,:)); % cecha

        % Ujednolicanie długości widm zerami
        wyglWidma(j,i).Af(1:lAf) = [wyglWidma(j,i).Af zeros(1, lAf-length(wyglWidma(j,i).Af ))];
        wyglWidma(j,i).Af2(1:lAf) = [wyglWidma(j,i).Af2 zeros(1, lAf-length(wyglWidma(j,i).Af2 ))];

        CentrWidm(j, kat).AfM= CentrWidm(j, kat).AfM +wyglWidma(j,i).Af'/wyglWidma(j,i).maxAf;
        CentrWidm(j, kat).Af2M=CentrWidm(j, kat).Af2M+wyglWidma(j,i).Af2'/wyglWidma(j,i).maxAf2; %mocy

        Ps=sum(wyglWidma(j,i).Af);
        CentrWidm(j, kat).AfE= CentrWidm(j, kat).AfE +wyglWidma(j,i).Af'/Ps;
        Ps=sum(wyglWidma(j,i).Af2);
        CentrWidm(j, kat).Af2E=CentrWidm(j, kat).Af2E+wyglWidma(j,i).Af2'/Ps;
        if isnan(CentrWidm(j, kat).AfM)
            e_ind = [j,i,kat]
        end
        %         Centr(j, kat).AfM = Centr(j, kat).AfM + CentrWidm(j,kat).AfM;
        %         Centr(j, kat).AfE = Centr(j, kat).AfE + CentrWidm(j,kat).AfE;

        if( length(nj) > 0 || plotAllFigures )

            subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
            title(fileSegMio(nrs)); xlabel(lingua.xl21) %end;
            subplot(2,2,kat+2);  hold on; plot(xf, wyglWidma(j,i).Af(f)/sum(wyglWidma(j,i).Af));
            %         subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
            %         subplot(2,2,2+kat);  hold on; plot(xf, wyglWidma(j,i).Af(f)/(Psyg(j,i);
            title(segment(nrs).miesien); xlabel(lingua.xl22)
        end
    end
    sgtitle(sprintf(lingua.sgt));%, v(fileSegNr(nrs)).infoRecord))
    for(kat = 1:2)
        CentrWidm(j, kat).AfM=CentrWidm(j, kat).AfM/nAf(kat); CentrWidm(j, kat).Af2M=CentrWidm(j,kat).Af2M/nAf(kat);
        CentrWidm(j, kat).AfE=CentrWidm(j, kat).AfE/nAf(kat); CentrWidm(j, kat).Af2E=CentrWidm(j,kat).Af2E/nAf(kat);

        if( length(nj) > 0 || plotAllFigures )
            subplot(2,2,kat); hold on; plot(xf, CentrWidm(j, kat).AfM(f),'k--'); hold off; axis('tight'); 
            title(fileSegMio(nrskat(kat)));
            subplot(2,2,kat+2); hold on; plot(xf, CentrWidm(j, kat).AfE(f),'k--'); hold off; axis('tight');
            title(fileSegMio(nrskat(kat))); xlabel(lingua.xl3)
        end
    end
    % TODO
    %     dCentrM(j,:)=abs(CentrWidm(j, 1).AfM'-CentrWidm(j, 2).AfM')/2;
    %     dCentrE(j,:)=abs(CentrWidm(j, 1).AfE'-CentrWidm(j, 2).AfE')/2;
    %City
    %     dCentrM(j,:)=abs(CentrWidm(j, 1).AfM'-CentrWidm(j, 2).AfM')/2;
    %     dCentrE(j,:)=abs(CentrWidm(j, 1).AfE'-CentrWidm(j, 2).AfE')/2;
    %Euclid
    %Cheby

    %         dCentr dla miejsciej eucli
    % TODO        subplot(4,2,kat+4);  hold on; plot(xf, CentrWidm(j,kat).AfM(f)'./dCentrM(j,f),'k--'); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]')
    %         subplot(4,2,kat+6);  hold on; plot(xf, dC(j,kat).AfE(f)'./dCentrE(j,f),'k--');  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]')

end

%średnia z centoidów
cwrr = []; cwrb = []; cwcr = []; cwcb = [];
cwrrE = []; cwrbE = []; cwcrE = []; cwcbE = [];
lpacj = 0;
for( i = 1:length(v))
    if(find(i==wybrJ)) continue; end
    for mkat = 1:2
        if v(i).infoTraining == 1 %&& v(i).mięsien
            if(mkat == 1)
                cwrr = [cwrr; CentrWidm(i,mkat).AfM']; % CentroidWidmposRedni
                cwrrE = [cwrrE; CentrWidm(i,mkat).AfE'];
            else
                cwrb = [cwrb; CentrWidm(i,mkat).AfM']; % CentroidWidmposRedniBiceps
                cwrbE = [cwrbE; CentrWidm(i,mkat).AfE'];
            end
        end
        if v(i).infoTraining == 2
            if(mkat == 1)
                %     size(CentrWidm(i,v(i).infoTraining).AfM)
                cwcr = [cwcr; CentrWidm(i,mkat).AfM']; % CentroidWidmpodChwyt
                cwcrE = [cwcrE; CentrWidm(i,mkat).AfE'];
            else
                cwcb = [cwcb; CentrWidm(i,mkat).AfM']; % CentroidWidmpodChwyt
                cwcbE = [cwcbE; CentrWidm(i,mkat).AfE'];
            end
        end
    end    
    lpacj = lpacj + 1;
end

CC =[ (mean(cwrr,1)); (mean(cwrb,1)); (mean(cwcr,1)); (mean(cwcb,1)); ]; % Centroidy Centroidów
CCE =[ (mean(cwrrE,1)); (mean(cwrbE,1)); (mean(cwcrE,1)); (mean(cwcbE,1)); ]; % Centroidy Centroidów

figure(nrF+1), subplot(2,1,1), plot(xf,CC(1,[1:length(xf)])); hold on; plot(xf,CC(2,[1:length(xf)])); title(lingua.ps); legend(txBR,txBB);hold off; axis tight; xlabel(lingua.fq); ylabel(lingua.amp);
figure(nrF+1), subplot(2,1,2), plot(xf,CC(3,[1:length(xf)])); hold on; plot(xf,CC(4,[1:length(xf)])); title(lingua.pc); legend(txBR,txBB);hold off; axis tight; xlabel(lingua.fq); ylabel(lingua.amp);
sgtitle(lingua.sgt41)

figure(nrF+2), hold on; plot(xf,CC(1,[1:length(xf)])); plot(xf,CC(2,[1:length(xf)])); plot(xf,CC(3,[1:length(xf)])); plot(xf,CC(4,[1:length(xf)])); title(lingua.t4);
legend(lingua.l4);xlabel (lingua.fq); ylabel(lingua.amp); axis tight; hold off;
%--------------------------------------------------------------------------
% dCentr posrRad, posrBiceps, podChwytRad, podChwytBiceps

figure(nrF+3), subplot(2,1,1), plot(xf, CCE(1,[1:length(xf)])); hold on; plot(xf, CCE(2,[1:length(xf)])); title(lingua.ps); legend(txBR,txBB);hold off; axis tight; xlabel(lingua.fq); ylabel(lingua.amp);
figure(nrF+3), subplot(2,1,2), plot(xf, CCE(3,[1:length(xf)])); hold on; plot(xf, CCE(4,[1:length(xf)])); title(lingua.pc); legend(txBR,txBB);hold off; axis tight; xlabel(lingua.fq); ylabel(lingua.amp);
sgtitle(lingua.sgt42)

figure(nrF+4), hold on; plot(xf,CCE(1,[1:length(xf)])); plot(xf,CCE(2,[1:length(xf)])); plot(xf,CCE(3,[1:length(xf)])); plot(xf,CCE(4,[1:length(xf)])); title(lingua.t6); hold off;
legend(lingua.l4); xlabel (lingua.fq); ylabel(lingua.amp); axis tight; hold off;
%--------------------------------------------------------------------------

%  rysowanie
for(j = 1:length(v))
    nj = find(wybrJ==j);
    for(kat=1:2)
        %Centr(j, kat).AfE = Centr(j,kat).AfE
        figure(nrF-50)
        if(v(j).infoTraining == 1) % pośredni
            subplot(4,2,kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); plot(xf, CC(kat, f),'r--','LineWidth',1); hold off; axis('tight'); xlabel(lingua.xl51); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
            subplot(4,2,kat+2);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--'); plot(xf, CCE(kat, f),'r--','LineWidth',1);  hold off; axis('tight'); xlabel(lingua.xl52); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
        else
            %         if(v(j).infoTraining == 2)
            subplot(4,2,4+kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); plot(xf, CC(kat+2, f),'r--','LineWidth',1); hold off; axis('tight'); xlabel(lingua.xl51); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
            subplot(4,2,6+kat);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--'); plot(xf, CCE(kat+2, f),'r--','LineWidth',1); hold off; axis('tight'); xlabel(lingua.xl52); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
        end
    end
    sgtitle(sprintf(lingua.sgt5, length(CentrWidm)));

    if( length(nj) > 0 || plotAllFigures )
        figure(j+nrF);
        for(kat = 1:2)
            subplot(2,2,kat);    hold on;
            plot(xf, CC(kat,f),'k--','LineWidth',1);
            plot(xf, CentrWidm(j,kat).AfM(f),'r--','LineWidth',1);
            subtitle(lingua.subt);
            hold off;
            %         title(fileSegMio(nrs)); xlabel("Widmo unormowane") %end;
            subplot(2,2,kat+2);  hold on;
            plot(xf, CCE(kat+2,f),'k--','LineWidth',1);
            plot(xf, CentrWidm(j,kat).AfE(f),'r--','LineWidth',1);
            hold off;

            %         subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
            %         subplot(2,2,2+kat);  hold on; plot(xf, wyglWidma(j,i).Af(f)/(Psyg(j,i);
            %         title(fileSegMio(nrs)); xlabel("Widmo mocy")
        end
    end
end

% if(v(j).infoTraining == 1) % pośredni
%             subplot(4,2,kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--','LineWidth',2); plot(xf, CC(kat, f),'r--','LineWidth',2); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
%             subplot(4,2,kat+2);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--','LineWidth',2); plot(xf, CC(kat, f),'r--','LineWidth',2);  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
%         else
% %         if(v(j).infoTraining == 2)
%             subplot(4,2,4+kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); plot(xf, CC(kat+2, f),'r--','LineWidth',2); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
%             subplot(4,2,6+kat);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--'); plot(xf, CC(kat+2, f),'r--','LineWidth',2); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
%         end
nrs = 0;  nf=2; nrF = nrF+length(v);
% odleglosci

% global Psyg, dEM, dCM, dists_chebyM;
for(j = 1:length(v)) % grupa trainig
    %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg))
        nrs = nseg(i);
        k = i;
        kat = segment(nrs).miesien;
        d=CentrWidm(j, kat).AfM-wyglWidma(j,k).Af'/wyglWidma(j,k).maxAf;
        dEM(j,k)=sqrt(sum(d.^2));
        dCM(j,k)=sum(abs(d));
        dists_chebyM(j,k) = max(abs(d));
        Psyg(j,k) = Esyg(j,k)/SygRawLen(nrs); % unormowany
        d=CentrWidm(j, kat).AfE-wyglWidma(j,k).Af'/Psyg(j,k);
        dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d)); dists_cheby(j,k) = max(abs(d));  %uwaga przesunięcie przecink
        d2=CentrWidm(j, kat).Af2M-wyglWidma(j,k).Af2';
        dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2)); dists_cheby2(j,k) = max(abs(d2));  %2-mocy
        dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
        %         figure(30), plot(j,dists_chebyM(j,k),'k.'); hold on;
    end % odległość w grupie

    figure(nrF-2), subplot(2,2,nf), plot(abs(d)); title("abs(d)"); hold on;
    for(i=j+1:length(v))
        %         abs(Esr(i)-Esr(j))/2
        %         dEsgr(j,i)=abs(Esr(i)-Esr(j))/2;
        %         for(kat = 1:2)
        %             % bezsensowne odejmowanie wzajemnych
        %             d=abs((CentrWidm(j,kat).AfM-CentrWidm(i, kat).AfM));  if(j==1) figure(nrF+j-3), subplot(1,2, 1), plot(d); title("distatans"); hold on; end
        %             dCG(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
        %             dEG(j,i)=sqrt(sum(d.^2))/2; dists_chebyG(j,k) = max(abs(d),[],1)/2; % odległosv mięfzy grupowa
        %             d=[];
        %             d=abs((CentrWidm(j, kat).Af2M-CentrWidm(i, kat).Af2M));  %if(j==1) figure(nrF+j), subplot(1,2, 1), plot(d); title("d. mocy"); hold on; end
        %             dCG2(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
        %             dEG2(j,i)=sqrt(sum(d.^2))/2; dists_chebyG2(j,k) = max(abs(d),[],1)/2;
        %         end
        % dzielimy prrze 2 aby odl. m. grupowe były lepiej porównywanlne z
        % liczonymi od centroidu
        nf=4;
    end
end
% figure(nrF-2), subplot(1,2, 1), plot(d); hold off;
% figure(nrF-2), subplot(1,2, 1), plot(d); hold off;
% disp("UWAGA CENrTroid")
% return;
nrs = 0; nf=2; %nrF = nrF+1;
% global Psyg, dEM, dCM, dists_chebyM; %ISTOTNE
if(1)
    for(j = 1:length(v)) % grupa training
        nseg=find(fileSegNr==j);
        for (i = 1:length(nseg))
            k = i;
            nrs = nseg(i);
            if(fileSegMio(nrs)==txBR)
                if v(j).infoTraining == 1
                    c=1;
                    %                 d=CC(c,:)-wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
                end
                if v(j).infoTraining == 2
                    c=2+1;
                    %                d=CC(c,:)-wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
                end
            end
            if(fileSegMio(nrs)==txBB)
                if v(j).infoTraining == 1
                    c=2;
                end
                if v(j).infoTraining == 2
                    c=2+2;
                end
            end
            d=CC(c,:)-wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
            kategoria = segment(nrs).miesien;
            dEM(j,k)=sqrt(sum(d.^2));
            dCM(j,k)=sum(abs(d));
            dists_chebyM(j,k) = max(abs(d));

            Ps = sum(wyglWidma(j,k).Af);
            Psyg(j,k) = Esyg(j,k)/SygRawLen(nrs); % unormowany
            d=CCE(c,:)-wyglWidma(j,k).Af/Ps;%Psyg(j,k);
            %         d=CentrWidm(j, kategoria).AfE-wyglWidma(j,k).Af'/Ps;%Psyg(j,k);
            dEE(j,k)=sqrt(sum(d.^2)); dCE(j,k)=sum(abs(d)); dists_chebyE(j,k) = max(abs(d));

            Ps2= sum(wyglWidma(j,k).Af2);
            d2=CentrWidm(j, kategoria).Af2E-wyglWidma(j,k).Af2'/Ps2; %TODO CCE2
            %         d2=CentrWidm(j, kategoria).Af2E-wyglWidma(j,k).Af2'/Ps2;
            dE2E(j,k)=sqrt(sum(d2.^2)); dC2E(j,k)=sum(abs(d2));
            dists_cheby2E(j,k) = max(abs(d2));  %2-mocy
            dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
        end % odległość w grupie
        figure(nrF-2), subplot(2,2,nf), plot(abs(d)); title("abs(d)"); hold on;
    end
end

% for i = length(wyglWidma) TEST(i) = wyglWidma(i,1).maxAf; end
% nag = ["między", "wew"];
% dists_cheby
% dists_chebyG
% iloczyn wektorywy tylko w przestrzeni euclidesa
if(printCentroids)
    % nag = ["między", "wew"]
    fprintf(1,'\n\teuc. Max\tCity\tCheby\tEnerg')
    for(j = 1:length(v)) % grupa
        fprintf(1,'\ngr.%-2d',j)
        %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
        if printCentroids
            for (k = 1:length(find(fileSegNr==j)))
                fprintf(1,';  %6.3f %.3f %.3f %8.3g',dEM(j,k),dCM(j,k),dists_chebyM(j,k),dEsyg(j,k)); %] [dE(2,1);dC(2,1);dists_cheby(2,1);;dEsyg(2,1)
            end
        end
        figure(nrF+90); plot(dEM(j,:),'k.')
    end
else
    disp("Pominęto wypisywanie odległości dla centroidów")
end

save centroids.mat CentrWidm wyglWidma dEM dCM dists_chebyM dEsyg Psyg dEM dCM dists_chebyM CC CCE Psr dC dE dists_cheby dEE dCE dists_chebyE dE2E dC2E dists_cheby2E dEsyg lpacj
toc;

function [d] = dictionary2(lang)
    PL = 1;
    EN = 2;
    dict(PL).t1 = "                                                                             Grupa dla ćwiczenia %s";
    dict(PL).xl1 = "Czas [s]";
    dict(PL).xl21 = "Widmo względem wartości maksymalnej";
    dict(PL).xl22 = "Widmo unormowane sumą amplitud";
    dict(PL).xl3 = "Widmo unormowane sumą amplitud";
    dict(PL).sgt = "Widma wewnątrzgrupowe";%, ćwiczenie %s";
    dict(PL).sgt41 = "Średnia z centroidów (Max)";
    dict(PL).ps = "Pośredni";
    dict(PL).pc = "Podchwyt";
    dict(PL).fq = "Częśtotliwość [Hz]";
    dict(PL).amp = "Amplituda";
    dict(PL).t4 = "Średnia z centroidów mięśni i ćwiczeń (Max)";
    dict(PL).l4 = ["Pośredni Radialis", "Pośredni Biceps", "Podchwyt Radialis", "Podchwyt Biceps"];
    dict(PL).sgt42 ="Średnia z centroidów (Energia)";
    dict(PL).sgt5 = "Centroidy z wszystkich ćwiczeń (%d) pośredni(1:4) i podchwyt(5:8)";
    dict(PL).xl51 = 'Widma wygładzone unorm.Max [Hz]';
    dict(PL).xl52 = 'Widma wygładzone unorm.Energia [Hz]';
    dict(PL).t6 = "Średnia z centroidów mięśni i ćwiczeń (Energia)";
    dict(PL).subt = "k--CC, r--CentrWidm";

    dict(EN).t1 = "                                                                             Group for training %s";
    dict(EN).xl1 = "Time [s]";
    dict(EN).xl21 = "Spectrum relative to the maximum value";
    dict(EN).xl22 = "Spectrum normalized by sum of amplitudes";
    dict(EN).xl3 = "Spectrum normalized by sum of amplitudes";
    dict(EN).sgt = "Intra-group spectra";%, training %s";
    dict(EN).sgt41 = "Mean of centroids (Max)";
    dict(EN).ps = "Intermediate";
    dict(EN).pc = "Supinated grip";
    dict(EN).fq = "Frequency [Hz]";
    dict(EN).amp = "Amplitude";
    dict(EN).t4 = "Muscle centroid and training average (Max)";
    dict(EN).l4 = [ strcat(dict(EN).ps, " Radialis"), strcat(dict(EN).ps, " Biceps"),
                    strcat(dict(EN).pc, " Radialis"), strcat(dict(EN).pc, " Biceps")];
    dict(EN).sgt42 ="Centroid average (Energy)";
    dict(EN).sgt5 = strcat("Centroids from all trainings (%d) ", dict(EN).ps,"(1:4) and ", dict(EN).pc,"(5:8)");
    dict(EN).xl51 = 'Smoothed normalized spectra by Max value [Hz]';
    dict(EN).xl52 = 'Smoothed normalized spectra by Energy [Hz]';
    dict(EN).t6 = "Muscle and training centroid average (Energy)";
    dict(EN).subt = "k--CC, r--CentrSpectra";
    
    d = dict(lang);
end
% TEST =[]; for i = 1:length(wyglWidma) TEST(i) = isempty(CentrWidm(i).AfM); end; max(TEST)