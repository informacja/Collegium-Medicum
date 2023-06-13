tic;
nrs = 0;
for(j = 1:2) % grupa
    figure(nrF-j);
    sgtitle(v(j).infoRecord);
    for (i = 1:length(find(fileSegNr==j)))
        nrs = nrs + 1;        
        if(fileSegMio(nrs) == txBR) subplot(2,2,1), plot(segment(nrs).data); hold on; title(txBR);end
        if(fileSegMio(nrs) == txBB) subplot(2,2,2), plot(segment(nrs).data); hold on; title(txBB);end
%         if(fileSegMio(nrs) == txBR) subplot(2,2,3), plot(segment(nrs).data); hold on; end
%         if(fileSegMio(nrs) == txBB) subplot(2,2,4), plot(segment(nrs).data); hold on; end
    end
    subplot(2,2,1), hold off; 
    subplot(2,2,2), hold off; 
    subplot(2,2,3), hold off; 
    subplot(2,2,4), hold off; 
    sgtitle(sprintf("Grupa dla ćwiczenia %s", v(j).infoRecord))
end

f = [1:8500+1]; xf = (f-1)/Tsyg;
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
    CentrWidm(j, 1).AfM=zeros(lAf,1); nAf(1)=0;
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
    
    if(j == length(v)) figure(j+nrF); end
    for (i = 1:length(find(fileSegNr==j))) 
        nrs = nrs + 1;
        
        if(fileSegMio(nrs) == txBR) kat = 1; end
        if(fileSegMio(nrs) == txBB) kat = 2; end
        wyglWidma(j,i).maxAf = max(wyglWidma(j,i).Af);
        wyglWidma(j,i).maxAf2 = max(wyglWidma(j,i).Af2);
        nAf(kat)=nAf(kat)+1;
        Psyg(j,i)=Esyg(j,i)/SygRawLen(nrs);

        CentrWidm(j, kat).AfM= CentrWidm(j, kat).AfM+wyglWidma(j,i).Af'/wyglWidma(j,i).maxAf; 
        CentrWidm(j, kat).Af2M=CentrWidm(j, kat).Af2M+wyglWidma(j,i).Af2'/wyglWidma(j,i).maxAf2; %mocy
        CentrWidm(j, kat).AfE= CentrWidm(j, kat).AfE+wyglWidma(j,i).Af'/Psyg(j,i); 
        CentrWidm(j, kat).Af2E=CentrWidm(j, kat).Af2E+wyglWidma(j,i).Af2'/Psyg(j,i); 

%         Centr(j, kat).AfM = Centr(j, kat).AfM + CentrWidm(j,kat).AfM;
%         Centr(j, kat).AfE = Centr(j, kat).AfE + CentrWidm(j,kat).AfE;

        if(j == length(v))
            subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
            title(fileSegMio(nrs)); xlabel("Widmo unormowane") %end; 
            subplot(2,2,kat+2);  hold on; plot(xf, wyglWidma(j,i).Af(f)/(Esyg(j,i)/SygRawLen(nrs))); 
    %         subplot(2,2,kat);    hold on; plot(xf, wyglWidma(j,i).Af(f)/wyglWidma(j,i).maxAf);
    %         subplot(2,2,2+kat);  hold on; plot(xf, wyglWidma(j,i).Af(f)/(Esyg(j,i)/SygRawLen(nrs)));
            title(fileSegMio(nrs)); xlabel("Widmo mocy")
        end
    end
    sgtitle(sprintf("Widma wewnątrzgrupowe, ćwiczenie %s", v(fileSegNr(nrs)).infoRecord))
    for(kat = 1:2)
        CentrWidm(j, kat).AfM=CentrWidm(j, kat).AfM/nAf(kat); CentrWidm(j, kat).Af2M=CentrWidm(j,kat).Af2M/nAf(kat);
        CentrWidm(j, kat).AfE=CentrWidm(j, kat).AfE/nAf(kat); CentrWidm(j, kat).Af2E=CentrWidm(j,kat).Af2E/nAf(kat);
    end
% TODO
%     dCentrM(j,:)=abs(CentrWidm(j, 1).AfM'-CentrWidm(j, 2).AfM')/2; 
%     dCentrE(j,:)=abs(CentrWidm(j, 1).AfE'-CentrWidm(j, 2).AfE')/2; 
    %City
%     dCentrM(j,:)=abs(CentrWidm(j, 1).AfM'-CentrWidm(j, 2).AfM')/2; 
%     dCentrE(j,:)=abs(CentrWidm(j, 1).AfE'-CentrWidm(j, 2).AfE')/2; 
    %Euclid
    %Cheby
    for(kat=1:2)
        %Centr(j, kat).AfE = Centr(j,kat).AfE 
        Psr(j) = mean(Psyg(j,:));
        Esr(j)=mean(Esyg(j,:)); % cecha
        figure(nrF+50)
        if(v(j).infoTraining == 1) % pośredni
            subplot(4,2,kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
            subplot(4,2,kat+2);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--');  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
        else
%         if(v(j).infoTraining == 2)
            subplot(4,2,4+kat);  hold on; plot(xf, CentrWidm(j,kat).AfM(f),'k--'); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
            subplot(4,2,6+kat);  hold on; plot(xf, CentrWidm(j,kat).AfE(f),'k--');  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]'); if kat == 1 title(v(1).infoRDisp); end; if kat == 2 title(v(1).infoBDisp); end
        end
        %         dCentr dla miejsciej eucli 
% TODO        subplot(4,2,kat+4);  hold on; plot(xf, CentrWidm(j,kat).AfM(f)'./dCentrM(j,f),'k--'); hold off; axis('tight'); xlabel('Widma wygładzone unorm.Max [hz]')
%         subplot(4,2,kat+6);  hold on; plot(xf, dC(j,kat).AfE(f)'./dCentrE(j,f),'k--');  hold off; axis('tight'); xlabel('Widma wygładzone unorm.Energia [Hz]')
    end
end
sgtitle(sprintf("Centroidy z wszystkich ćwiczeń (%d) pośredni i podchwyt)", length(CentrWidm)));

%średnia z centoidów
cwrr = []; cwrb = []; cwcr = []; cwcb = [];
for( i = 1:length(CentrWidm))
    for mkat = 1:2
        if v(i).infoTraining == 1 %&& v(i).mięsien
            if(mkat == 1)
                cwrr = [cwrr; CentrWidm(i,mkat).AfM']; % CentroidWidmposRedni
            else
                cwrb = [cwrb; CentrWidm(i,mkat).AfM']; % CentroidWidmposRedniBiceps
            end
        end
        if v(i).infoTraining == 2
            if(mkat == 1)
    %     size(CentrWidm(i,v(i).infoTraining).AfM)
                cwcr = [cwcr; CentrWidm(i,mkat).AfM']; % CentroidWidmpodChwyt
            else
                cwcb = [cwcb; CentrWidm(i,mkat).AfM']; % CentroidWidmpodChwyt
            end
        end
    end
end

figure(nrF+1), subplot(2,1,1), plot(mean(cwrr)); hold on; plot(mean(cwrb)); title("Posredni");
legend(txBR,txBB);hold off;
figure(nrF+1), subplot(2,1,2), plot(mean(cwcr(:,:))); hold on; plot(mean(cwcb(:,:))); title("Podchwyt");
legend(txBR,txBB);hold off;
sgtitle("Średnia z centroidów")

figure(nrF+2), hold on; plot(mean(cwrr(:,:))); plot(mean(cwrb(:,:)));
plot(mean(cwcr(:,:))); plot(mean(cwcb(:,:))); title("średnia z centroidów mięśni i ćwiczeń"); 
legend("Pośredni Radialis", "Pośredni Biceps", "Podchwyt Radialis", "Podchwyt Biceps");
CC =[ (mean(cwrr)); (mean(cwrb)); (mean(cwcr)); (mean(cwcb)); ]; % Centroidy Centroidów
%--------------------------------------------------------------------------
 % dCentr posrRad, posrBiceps, podChwytRad, podChwytBiceps

nrs = 0;  nf=2; nrF = nrF+100;
% odleglosci

% global Psyg, dEM, dCM, dists_chebyM;
for(j = 1:length(v)) % grupa trainig
    %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
    for (k = 1:length(find(fileSegNr==j))) %
        nrs = nrs + 1;
        d=CentrWidm(j).AfM-wyglWidma(j,k).Af'/wyglWidma(j,k).maxAf; 
        dEM(j,k)=sqrt(sum(d.^2)); 
        dCM(j,k)=sum(abs(d))/100; 
        dists_chebyM(j,k) = max(abs(d),[],1);  %uwaga przesunięcie przecinka
        Psyg(j,k) = Esyg(j,k)/SygRawLen(nrs); % unormowany
        d=CentrWidm(j).AfE-wyglWidma(j,k).Af'/Psyg(j,k); 
        dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d))/100; dists_cheby(j,k) = max(abs(d));  %uwaga przesunięcie przecink
        d2=CentrWidm(j).Af2M-wyglWidma(j,k).Af2'; 
        dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2))/100; dists_cheby2(j,k) = max(abs(d2));  %2-mocy
        dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
%         figure(30), plot(j,dists_chebyM(j,k),'k.'); hold on;
    end % odległość w grupie
    
    figure(nrF-2), subplot(2,2,nf), plot(abs(d)); title("abs(d)"); hold on;
    for(i=j+1:length(v))
%         abs(Esr(i)-Esr(j))/2
        dEsgr(j,i)=abs(Esr(i)-Esr(j))/2;
        d=abs((CentrWidm(j).AfM-CentrWidm(i).AfM));  if(j==1) figure(nrF+j-3), subplot(1,2, 1), plot(d); title("distatans"); hold on; end
        dCG(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
        dEG(j,i)=sqrt(sum(d.^2))/2; dists_chebyG(j,k) = max(abs(d),[],1)/2; % odległosv mięfzy grupowa
        d=[];
        d=abs((CentrWidm(j).Af2M-CentrWidm(i).Af2M));  %if(j==1) figure(nrF+j), subplot(1,2, 1), plot(d); title("d. mocy"); hold on; end
        dCG2(j,i)=sum(d)/200; % UWAGA przesunięcie przecinka o dwa miejsca w celu łatwiejszej interpretacji wyników
        dEG2(j,i)=sqrt(sum(d.^2))/2; dists_chebyG2(j,k) = max(abs(d),[],1)/2;
        % dzielimy prrze 2 aby odl. m. grupowe były lepiej porównywanlne z
        % liczonymi od centroidu
        nf=4;
    end
end
% figure(nrF-2), subplot(1,2, 1), plot(d); hold off;
% figure(nrF-2), subplot(1,2, 1), plot(d); hold off;

nrs = 0; nf=2; %nrF = nrF+1;
% global Psyg, dEM, dCM, dists_chebyM;
for(j = 1:length(v)) % grupa training
    for (k = 1:length(find(fileSegNr==j))) % training
        nrs = nrs + 1;
        if(fileSegMio(nrs)==txBR)
            if v(j).infoTraining == 1
                c=1;
                d=CC(c)-wyglWidma(j,k).Af'/wyglWidma(j,k).maxAf; 
            end
            if v(j).infoTraining == 2
               c=2+1;
               d=CC(c)-wyglWidma(j,k).Af'/wyglWidma(j,k).maxAf; 
            end
        end
        if(fileSegMio(nrs)==txBB)
            if v(j).infoTraining == 1
                c=2;
                d=CC(c)-wyglWidma(j,k).Af'/wyglWidma(j,k).maxAf; 
            end
            if v(j).infoTraining == 2
               c=2+2;
               d=CC(c)-wyglWidma(j,k).Af'/wyglWidma(j,k).maxAf; 
            end
        end
        dEM(j,k)=sqrt(sum(d.^2)); 
        dCM(j,k)=sum(abs(d))/100; 
        dists_chebyM(j,k) = max(abs(d),[],1);  %uwaga przesunięcie przecinka
        Psyg(j,k) = Esyg(j,k)/SygRawLen(nrs); % unormowany
        d=CentrWidm(j).AfE-wyglWidma(j,k).Af'/Psyg(j,k); 
        dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d))/100; dists_cheby(j,k) = max(abs(d),[],1);  %uwaga przesunięcie przecink
        d2=CentrWidm(j).Af2M-wyglWidma(j,k).Af2'; 
        dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2))/100; dists_cheby2(j,k) = max(abs(d2),[],1);  %2-mocy
        dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
    end % odległość w grupie
    figure(nrF-2), subplot(2,2,nf), plot(abs(d)); title("abs(d)"); hold on;
end
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
    disp("Pominęto wyposywanie odległości dla centroidów")
end

% cc4 %nadpisz odległości

% odległości centroidów
        % grupa 1, 2 , w gr. 1 2
 
save centroids.mat CentrWidm dEM dCM dists_chebyM dEsyg Psyg dEM dCM dists_chebyM
toc;