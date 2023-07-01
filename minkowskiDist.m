% 0 od. centoidów J i Max,
% 1 od centroidów J i sum(Af)
% 2 od CC J  i Max
% 3 od CC i cenroidów J
nrs = 0; nf=2; %nrF = nrF+1;

wybrJakieDist = (jakieDist == 1 || jakieDist == 3);
% global Psyg, dEM, dCM, dists_chebyM; %ISTOTNE
tic;
ddEM = []; ddCM = []; ddists_chebyM = []; dAll1 = []; dAll2 = []; dAll3 = []; dAll4=[];%zeros(1 length(SygKat)*length(CentrWidm(1,1).AfM));
if(1)
    jj= zeros(1,2); 
   
    for(j = 1:length(v)) % grupa training
        nseg=find(fileSegNr==j);
        for (i = 1:length(nseg))
            k = i;
            nrs = nseg(i);
            c = SygKat(nrs);  % index
            nrG = v(j).infoTraining; % Gestu
            kat = segment(nrs).miesien;
            jj(nrG) = jj(nrG)+1;
            switch (jakieDist)
                case 1 % CWłasnego
                    Afw = CentrWidm(j,kat).AfM';
                    Af = wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
                    Psyg(j,k) = Esyg(j,k)/SygRawLen(nrs);
                case 2 % EE
                    Afw = CC(c,:);
                    Af = wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
                case 3 % E
                    Ps = sum(wyglWidma(j,k).Af);
                    Afw = CentrWidm(j,kat).AfE';
                    Af = wyglWidma(j,i).Af/Ps;
                case 4
                    Ps = sum(wyglWidma(j,k).Af);
                    Afw = CCE(c,:);
                    Af = wyglWidma(j,i).Af/Ps; % powtórka
                case 5
                    Afw = CC(c,:);
                    Af = CentrWidm(j,kat).AfM';
                case 6
                    Afw = CCE(c,:);
                    Af = CentrWidm(j,kat).AfE'; %Chyba błąd
            end
            d=Afw-Af; % wzorcowe

            dE=sqrt(sum(d.^2));
            dC=sum(abs(d));
            dCZ = max(abs(d));
                  
            if(wybrJakieDist)
                switch(c) %jakieDist 1 i 3 oraz d E C CZ
                    case 1, dAll1 = [dAll1 dCZ];  
                    case 2, dAll2 = [dAll2 dCZ];  
                    case 3, dAll3 = [dAll3 dCZ];  
                    case 4, dAll4 = [dAll4 dCZ]; 
                end   
            end
                   
            switch (jakieDist)
                case 1
                    dCM(j,k)= dC;
                    dEM(j,k)= dE;
                    dists_chebyM(j,k) = dCZ;
                case 2
                    dCM(j,k)= dC;
                    dEM(j,k)= dE;
                    dists_chebyM(j,k) = dCZ;
                case 3
                    dCE(j,k)= dC;
                    dEE(j,k)= dE;
                    dists_chebyE(j,k) = dCZ;
                case 4
                    dCE(j,k)= dC;
                    dEE(j,k)= dE;                    
                    dists_chebyE(j,k) = dCZ;
                case 5
                    ddCM(nrG, jj(nrG), kat) = dC;
                    ddEM(nrG, jj(nrG), kat) = dE;
                    ddists_chebyM(nrG, jj(nrG), kat) = dCZ;
                case 6                    
                    ddCE(nrG, jj(nrG), kat)= dC;
                    ddEE(nrG, jj(nrG), kat)= dE;
                    ddists_chebyE(nrG, jj(nrG), kat) = dCZ;
            end

            % unormowany
            %         d=CCE(c,:)-wyglWidma(j,k).Af/Ps;%Psyg(j,k);
            % %         d=CentrWidm(j, kategoria).AfE-wyglWidma(j,k).Af'/Ps;%Psyg(j,k);
            %         dEE(j,k)=sqrt(sum(d.^2)); dCE(j,k)=sum(abs(d)); dists_chebyE(j,k) = max(abs(d));
            %
            %         Ps2= sum(wyglWidma(j,k).Af2);
            %         d2=CentrWidm(j, kategoria).Af2E-wyglWidma(j,k).Af2'/Ps2; %TODO CCE2
            % %         d2=CentrWidm(j, kategoria).Af2E-wyglWidma(j,k).Af2'/Ps2;
            %         dE2E(j,k)=sqrt(sum(d2.^2)); dC2E(j,k)=sum(abs(d2));
            %         dists_cheby2E(j,k) = max(abs(d2));  %2-mocy
            %         dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
        end % odległość w grupie
%         figure(nrF-2), subplot(2,2,nf), plot(abs(d)); title("abs(d)"); hold on;
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
        %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhattan)
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
txDist = "Niezdefinowano";
if (jakieDist == 1) txDist = "Moc (Czebyszew)"; end
if (jakieDist == 3) txDist = "Energia (Czebyszew)"; end

lbins = 30;
if(wybrJakieDist)

figure(nrFig-1);
    for( c = 1:4 )
        subplot(2,2,c);
        switch(c)
            case 1, nd = find(abs(dAll1)>0.1e-16); mhistf(dAll1(nd),lbins,1); nd = dAll1(nd); hold on; axis('tight'); subtitle(strcat("SygKat = ", string(c))); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; xlabel(sprintf("$\\bar{x}$ = %g $\\sigma$ = %g", m, s),'Interpreter','latex');
            case 2, nd = find(abs(dAll2)>0.1e-16); mhistf(dAll2(nd),lbins,1); nd = dAll2(nd); hold on; axis('tight'); subtitle(strcat("SygKat = ", string(c))); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; xlabel(sprintf("$\\bar{x}$ = %g $\\sigma$ = %g", m, s),'Interpreter','latex');
            case 3, nd = find(abs(dAll3)>0.1e-16); mhistf(dAll3(nd),lbins,1); nd = dAll3(nd); hold on; axis('tight'); subtitle(strcat("SygKat = ", string(c))); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; xlabel(sprintf("$\\bar{x}$ = %g $\\sigma$ = %g", m, s),'Interpreter','latex');
            case 4, nd = find(abs(dAll4)>0.1e-16); mhistf(dAll4(nd),lbins,1); nd = dAll4(nd); hold on; axis('tight'); subtitle(strcat("SygKat = ", string(c))); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; xlabel(sprintf("$\\bar{x}$ = %g $\\sigma$ = %g", m, s),'Interpreter','latex');
        end
        sgtitle(txDist)
    end
    figure(jakieDist+250), nd = [dAll1 dAll2 dAll3 dAll4]; mhistf(nd,lbins,1); hold on; axis('tight'); title(sprintf("Wykres zbiorczy, jakieDist = %d, %s", jakieDist, txDist));  
    m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; xlabel(sprintf("$\\bar{x}$ = %g $\\sigma$ = %g", m, s),'Interpreter','latex');
end
% [pvg,pvr,chi2e,chi2r,bins,Nemp]=mhistf(dAll1(nd),lbins,1); hold off;
% save centroids.mat CentrWidm dEM dCM dists_chebyM dEsyg Psyg dEM dCM dists_chebyM CC CCE Psr dC dE dists_cheby dEE dCE dists_chebyE dE2E dC2E dists_cheby2E dEsyg mx
toc;