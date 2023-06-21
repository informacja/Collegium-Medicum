tic;
wybrJ = [];
for (j=1:length(v))
    if (v(j).infoRecord == """22 pośredni NORAXON ELEKTRODY """) wybrJ = [wybrJ j]; end
    if (v(j).infoRecord == """22 podchwyt NORAXON ELEKTRODY""") wybrJ = [wybrJ j]; end
    if (v(j).infoRecord == """11 pośredni""") wybrJ = [wybrJ j]; end
    if (v(j).infoRecord == """11 podchwyt""") wybrJ = [wybrJ j]; end 
end
if(size(CentrWidm)==[2, length(Widma)])
    CentrWidm = CentrWidm';
end
nrs = 0;
figure(nrF);
tmpTxt = "";
for(j = 1:2) % grupa
    tmpTxt = [tmpTxt;v(j).infoRecord];
     nseg=find(fileSegNr==j);
    for (i = 1:length(nseg))
        nrs = nseg(i);    nsp=SygKat(nrs);  
%         if(fileSegMio(nrs) == txBR) 
            subplot(2,2,nsp), plot(segment(nrs).data); hold on; title(nsp);
%     end
%         if(fileSegMio(nrs) == txBB) subplot(2,2,ns), plot(segment(nrs).data); hold on; title(txBB);end
%         if(fileSegMio(nrs) == txBR) subplot(2,2,3), plot(segment(nrs).data); hold on; end
%         if(fileSegMio(nrs) == txBB) subplot(2,2,4), plot(segment(nrs).data); hold on; end
    end
    subplot(2,2,1), hold off; 
    subplot(2,2,2), hold off; 
    subplot(2,2,3), hold off; 
    subplot(2,2,4), hold off; 
    sgtitle(sprintf("Grupa dla ćwiczenia %s", v(j).infoRecord))
end
sgtitle(tmpTxt);

f = [1:8500+1]; f = [1:6200+1]; xf = (f-1)/Tsyg;
clear dCentrM dCentrE;
nf=2;  
nrs = 0;
nrF = nrF + length(v);
% normowanie widma
for(j = 1:length(v)) % grupa
%     figure(nrFw+j);close(nrFw+j);
%     figure(nrF+j); 
       % normowanie widma
    lAf=length(wyglWidma(j,1).Af);
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
%     Centr(j, 2).AfM = zeros(lAf,1);
%     Centr(j, 2).AfE = zeros(lAf,1);
    nj = find(wybrJ==j);
    if( length(nj) > 0 || plotAllFigures) figure(j+nrF); end
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg)) 
        nrs = nseg(i);
        
%         if(v(j).plikeSegMio(i) == txBR) kat = 1; end

        if(fileSegMio(nrs) == txBR) kat = 1; end
        if(fileSegMio(nrs) == txBB) kat = 2; end
        wyglWidma(j,i).maxAf = max(wyglWidma(j,i).Af);
        wyglWidma(j,i).maxAf2 = max(wyglWidma(j,i).Af2);
%         if (wyglWidma() todo
        nAf(kat)=nAf(kat)+1;

        Psyg(j,i)=Esyg(j,i)/SygRawLen(nrs);
        Psr(j) = mean(Psyg(j,:)); % Power
        Esr(j) = mean(Esyg(j,:)); % cecha
        
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
            title(fileSegMio(nrs)); xlabel("Widmo unormowane") %end; 
            subplot(2,2,kat+2);  hold on; plot(xf, wyglWidma(j,i).Af(f)/sum(wyglWidma(j,i).Af)); 
    %         subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
    %         subplot(2,2,2+kat);  hold on; plot(xf, wyglWidma(j,i).Af(f)/(Psyg(j,i);
            title(segment(nrs).miesien); xlabel("Widmo mocy")
        end
    end
    sgtitle(sprintf("Widma wewnątrzgrupowe, ćwiczenie %s", v(fileSegNr(nrs)).infoRecord))
    for(kat = 1:2)
        CentrWidm(j, kat).AfM=CentrWidm(j, kat).AfM/nAf(kat); CentrWidm(j, kat).Af2M=CentrWidm(j,kat).Af2M/nAf(kat);
        CentrWidm(j, kat).AfE=CentrWidm(j, kat).AfE/nAf(kat); CentrWidm(j, kat).Af2E=CentrWidm(j,kat).Af2E/nAf(kat);
    
        if( length(nj) > 0 || plotAllFigures )
            subplot(2,2,kat); hold on; plot(xf, CentrWidm(j, kat).AfM(f),'k--'); hold off; axis('tight');
            subplot(2,2,kat+2); hold on; plot(xf, CentrWidm(j, kat).AfE(f),'k--'); hold off; axis('tight');
            title(fileSegMio(nrs)); xlabel("Widmo mocy")
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
for( i = 1:length(v))
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
end
CC =[ (mean(cwrr)); (mean(cwrb)); (mean(cwcr)); (mean(cwcb)); ]; % Centroidy Centroidów
CCE =[ (mean(cwrrE)); (mean(cwrbE)); (mean(cwcrE)); (mean(cwcbE)); ]; % Centroidy Centroidów

figure(nrF+1), subplot(2,1,1), plot(CC(1,:)); hold on; plot(CC(2,:)); title("Pośredni"); legend(txBR,txBB);hold off;
figure(nrF+1), subplot(2,1,2), plot(CC(3,:)); hold on; plot(CC(4,:)); title("Podchwyt"); legend(txBR,txBB);hold off;
sgtitle("Średnia z centroidów (Max)")

figure(nrF+2), hold on; plot(CC(1,:)); plot(CC(2,:)); plot(CC(3,:)); plot(CC(4,:)); title("Średnia z centroidów mięśni i ćwiczeń (Max)"); 
legend("Pośredni Radialis", "Pośredni Biceps", "Podchwyt Radialis", "Podchwyt Biceps");
%--------------------------------------------------------------------------
 % dCentr posrRad, posrBiceps, podChwytRad, podChwytBiceps

figure(nrF+3), subplot(2,1,1), plot(CCE(1,:)); hold on; plot(CCE(2,:)); title("Pośredni"); legend(txBR,txBB);hold off;
figure(nrF+3), subplot(2,1,2), plot(CCE(3,:)); hold on; plot(CCE(4,:)); title("Podchwyt"); legend(txBR,txBB);hold off;
sgtitle("Średnia z centroidów (Energia)")

figure(nrF+4), hold on; plot(CCE(1,:)); plot(CCE(2,:)); plot(CCE(3,:)); plot(CCE(4,:)); title("Średnia z centroidów mięśni i ćwiczeń (Energia)"); 
legend("Pośredni Radialis", "Pośredni Biceps", "Podchwyt Radialis", "Podchwyt Biceps");
%--------------------------------------------------------------------------

%  rysowanie
for(j = 1:length(v))
    nj = find(wybrJ==j);
        for(kat=1:2)
            %Centr(j, kat).AfE = Centr(j,kat).AfE 
            figure(nrF-50)
            if(v(j).infoTraining == 1) % pośredni
                subplot(4,2,kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); plot(xf, CC(kat, f),'r--','LineWidth',1); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
                subplot(4,2,kat+2);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--'); plot(xf, CCE(kat, f),'r--','LineWidth',1);  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
            else
    %         if(v(j).infoTraining == 2)
                subplot(4,2,4+kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); plot(xf, CC(kat+2, f),'r--','LineWidth',1); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
                subplot(4,2,6+kat);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--'); plot(xf, CCE(kat+2, f),'r--','LineWidth',1); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
            end
         end
    sgtitle(sprintf("Centroidy z wszystkich ćwiczeń (%d) pośredni(1:4) i podchwyt(4:8)", length(CentrWidm)));

    if( length(nj) > 0 || plotAllFigures )
        figure(j+nrF);
        for(kat = 1:2)
            subplot(2,2,kat);    hold on; 
            plot(xf, CC(kat,f),'k--','LineWidth',1); 
            plot(xf, CentrWidm(j,kat).AfM(f),'r--','LineWidth',1); 
            subtitle("k--CC, r--CentrWidm");
            hold off;
    %         title(fileSegMio(nrs)); xlabel("Widmo unormowane") %end;
            subplot(2,2,kat+2);  hold on; 
            plot(xf,CCE(kat+2,f),'k--','LineWidth',1);
            plot(xf, CentrWidm(j,kat).AfE(f),'r--','LineWidth',1);
            hold off;

            %         subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
            %         subplot(2,2,2+kat);  hold on; plot(xf, wyglWidma(j,i).Af(f)/(Psyg(j,i);
    %         title(fileSegMio(nrs)); xlabel("Widmo mocy")
        end
    end
end

% if(v(j).infoTraining == 1) % pośredni
%             subplot(4,2,kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--','LineWidth',2); plot(xf, CC(kat, f),'r--','LineWidth',2); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
%             subplot(4,2,kat+2);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--','LineWidth',2); plot(xf, CC(kat, f),'r--','LineWidth',2);  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
%         else
% %         if(v(j).infoTraining == 2)
%             subplot(4,2,4+kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); plot(xf, CC(kat+2, f),'r--','LineWidth',2); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
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

save centroids.mat CentrWidm wyglWidma dEM dCM dists_chebyM dEsyg Psyg dEM dCM dists_chebyM CC CCE Psr dC dE dists_cheby dEE dCE dists_chebyE dE2E dC2E dists_cheby2E dEsyg  
toc;

% TEST =[]; for i = 1:length(wyglWidma) TEST(i) = isempty(CentrWidm(i).AfM); end; max(TEST)