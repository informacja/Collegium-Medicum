% 0 od. centoidów J i Max,
% 1 od centroidów J i sum(Af)
% 2 od CC J  i Max
% 3 od CC i cenroidów J
nrs = 0; nf=2; %nrF = nrF+1;

lingua = dictionary2(lang);

% wybrJakieDist = (jakieDist == 1 || jakieDist == 3);
% global Psyg, dEM, dCM, dists_chebyM; %ISTOTNE
tic;
for(nband=1:twoOrOne)
 Sb(nband).ddEM = []; Sb(nband).ddCM = []; Sb(nband).ddists_chebyM = []; %zeros(1 length(SygKat)*length(CentrWidm(1,1).AfM));
end; dAll1 = []; dAll2 = []; dAll3 = []; dAll4=[];
if(1)
    maxGestures = 2;
    jj= zeros(1,maxGestures); 
   
    for(j = 1:length(v)) % grupa training
        nseg=find(fileSegNr==j);
        for (i = 1:length(nseg))
            k = i;
            nrs = nseg(i);
            c = SygKat(nrs);  % index
            if( v(j).infoTraining > maxGestures) nrG = mod(v(j).infoTraining,maxGestures); % Gest
            else nrG = v(j).infoTraining; end
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
            
            % [dCM,dEM,dists_chebyM,dCE,dEE,dists_chebyE,ddCM,ddEM,ddists_chebyM,ddCE,ddEE,ddists_chebyE] = segDist(Afw(1:nflim),Af(1:nflim),wybrJakieDist,jakieDist,c,dAll1,dAll2,dAll3,dAll4,nrG,jj,kat);
             % podzielWidmo()
            % 
            for(nband=1:twoOrOne)
                if(nband-1)
                    d=Afw(nflim+1:end)-Af(nflim+1:end);
                else
                    d=Afw(1:nflim)-Af(1:nflim);
                end

                dC=sum(abs(d));
                dE=sqrt(sum(d.^2));
                dCZ = max(abs(d));

                if(wybrJakieDist)
                    switch(c) %jakieDist 1 i 3 oraz d E C CZ
                        case 1, dAll1 = [dAll1 dE];
                        case 2, dAll2 = [dAll2 dE];
                        case 3, dAll3 = [dAll3 dE];
                        case 4, dAll4 = [dAll4 dE];
                    end
                end

                switch (jakieDist)
                    case 1
                        Sb(nband).dCM(j,k)= dC;
                        Sb(nband).dEM(j,k)= dE;
                        Sb(nband).dists_chebyM(j,k) = dCZ;
                    case 2
                        Sb(nband).dCM(j,k)= dC;
                        Sb(nband).dEM(j,k)= dE;
                        Sb(nband).dists_chebyM(j,k) = dCZ;
                    case 3
                        Sb(nband).dCE(j,k)= dC;
                        Sb(nband).dEE(j,k)= dE;
                        Sb(nband).dists_chebyE(j,k) = dCZ;
                    case 4
                        Sb(nband).dCE(j,k)= dC;
                        Sb(nband).dEE(j,k)= dE;
                        Sb(nband).dists_chebyE(j,k) = dCZ;
                    case 5
                        Sb(nband).ddCM(nrG, jj(nrG), kat) = dC;
                        Sb(nband).ddEM(nrG, jj(nrG), kat) = dE;
                        Sb(nband).ddists_chebyM(nrG, jj(nrG), kat) = dCZ;
                    case 6
                        Sb(nband).ddCE(nrG, jj(nrG), kat)= dC;
                        Sb(nband).ddEE(nrG, jj(nrG), kat)= dE;
                        Sb(nband).ddists_chebyE(nrG, jj(nrG), kat) = dCZ;
                end
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
    
            % [dCM,dEM,dists_chebyM,dCE,dEE,dists_chebyE,ddCM,ddEM,ddists_chebyM,ddCE,ddEE,ddists_chebyE] = segDist(Afw(nflim+1:end),Af(nflim+1:end),wybrJakieDist,jakieDist,c,dAll1,dAll2,dAll3,dAll4,nrG,jj,kat);
            % nrF=nrF+1; figure(nrF); disppolt;
           
end

% for i = length(wyglWidma) TEST(i) = wyglWidma(i,1).maxAf; end
% nag = ["między", "wew"];
% dists_cheby
% dists_chebyG
% iloczyn wektorywy tylko w przestrzeni euclidesa
if(printCentroids)
    % nag = ["między", "wew"]
    fprintf(1,'\n\tEuc. Max\tCity\tCheby\tEnerg')
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
txDist = "Niezdefinowano"; txDist = "";
if (jakieDist == 1) txDist = lingua.me; end
if (jakieDist == 3) txDist = lingua.ee; end
txpi = ''; yTxt="Count"; xTxt="Distances"; 
lbins = 20;
if(wybrJakieDist)

    figure(nrFig-1+jakieDist); if(lang==PL) fprintf("nrFig = %d, jakieDist = %d flagaMaxima = %d\n", nrFig-1+jakieDist, jakieDist, flagaMaxima); end
    for( c = 1:4 )
        subplot(2,2,c); eps= 0.1e-16;
        switch(c)
            case 1, nd = find(abs(dAll1)>eps); [pvg,pve,pvM]=mhistMGE(dAll1(nd),lbins); nd = dAll1(nd); hold on; axis('tight'); subtitle( char(c+96)+")"); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; ylabel(yTxt); xlabel(xTxt); 
            % results = mhistMGEfromLMA(nd, lbins, false); 
            % pvgeMa = [ results.Normal.p_chi2, results.DoubleExponential.p_chi2 results.Maxwell.p_chi2]
            if(lang==PL) fprintf("$\\bar{x}$ = %g $\\sigma$ = %g $\\pi_g$=%.2f$\\%%$  $\\pi_e$=%.2f$\\%%$  $\\pi_M$=%.2f$\\%%$ \n", m, s, pvg*100, pve*100, pvM*100); end%,'Interpreter','latex'); end
            case 2, nd = find(abs(dAll2)>eps); [pvg,pve,pvM]=mhistMGE(dAll2(nd),lbins); nd = dAll2(nd); hold on; axis('tight'); subtitle( char(c+96)+")"); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; ylabel(yTxt); xlabel(xTxt);
            % results = mhistMGEfromLMA(nd, lbins, false); 
            % pvgeMb = [ results.Normal.p_chi2, results.DoubleExponential.p_chi2 results.Maxwell.p_chi2]  
            if(lang==PL) fprintf("$\\bar{x}$ = %g $\\sigma$ = %g $\\pi_g$=%.2f$\\%%$  $\\pi_e$=%.2f$\\%%$  $\\pi_M$=%.2f$\\%%$  \n", m, s, pvg*100, pve*100, pvM*100); end%,'Interpreter','latex'); end
            case 3, nd = find(abs(dAll3)>eps); [pvg,pve,pvM]=mhistMGE(dAll3(nd),lbins); nd = dAll3(nd); hold on; axis('tight'); subtitle( char(c+96)+")"); 
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; ylabel(yTxt); xlabel(xTxt); 
            % results = mhistMGEfromLMA(nd, lbins, false); 
            % pvgeMc = [ results.Normal.p_chi2, results.DoubleExponential.p_chi2 results.Maxwell.p_chi2]     
            if(lang==PL) fprintf("$\\bar{x}$ = %g $\\sigma$ = %g $\\pi_g$=%.2f$\\%%$  $\\pi_e$=%.2f$\\%%$  $\\pi_M$=%.2f$\\%%$ \n", m, s, pvg*100, pve*100, pvM*100); end %,'Interpreter','latex'); end
            case 4, nd = find(abs(dAll4)>eps); [pvg,pve,pvM]=mhistMGE(dAll4(nd),lbins); nd = dAll4(nd); hold on; axis('tight'); subtitle( char(c+96)+")");   
            m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; ylabel(yTxt); xlabel(xTxt); 
            % results = mhistMGEfromLMA(nd, lbins, true); 
            % pvgeMd = [ results.Normal.p_chi2, results.DoubleExponential.p_chi2 results.Maxwell.p_chi2]
      
            if(lang==PL) fprintf("$\\bar{x}$ = %g $\\sigma$ = %g $\\pi_g$=%.2f$\\%%$  $\\pi_e$=%.2f$\\%%$  $\\pi_M$=%.2f$\\%%$ \n", m, s, pvg*100, pve*100, pvM*100); end %,'Interpreter','latex');end
        end
        sgtitle(txDist)
    end    
    legend(["" "Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], "Position", [0.2965    0.3238    0.2228    0.0903] )
    figure(jakieDist+250), nd = [dAll1 dAll2 dAll3 dAll4]; [pvg,pve,pvM]=mhistMGE(nd,lbins); hold on; axis('tight'); title(sprintf(lingua.t1, jakieDist, txDist));  subtitle("a)"); ylabel(yTxt); xlabel(xTxt);
    m = mean(nd); s = std(nd); ax = axis; plot([s+m s+m], ax(3:4), 'k--'); plot([m-s m-s], ax(3:4), 'k--'); plot([m m], ax(3:4), "m--");  hold off; %xlabel(sprintf("$\\bar{x}$ = %g $\\sigma$ = %g $\\pi_g$=%.2f$\\%%$  $\\pi_e$=%.2f$\\%%$  $\\pi_M$=%.2f$\\%%$ ", m, s, pvg*100, pve*100, pvM*100),'Interpreter','latex');
end
% [pvg,pve,pvM,bins,Nemp]=mhistMGE(dAll1(nd),lbins); hold off;
% save centroids.mat CentrWidm dEM dCM dists_chebyM dEsyg Psyg dEM dCM dists_chebyM CC CCE Psr dC dE dists_cheby dEE dCE dists_chebyE dE2E dC2E dists_cheby2E dEsyg mx
toc;

% function [dCM,dEM,dists_chebyM,dCE,dEE,dists_chebyE,ddCM,ddEM,ddists_chebyM,ddCE,ddEE,ddists_chebyE] = segDist(Afw,Af,wybrJakieDist,jakieDist,c,dAll1,dAll2,dAll3,dAll4,nrG,jj,kat)
% % dCM=[]; dEM=[]; dists_chebyM=[]; dCE=[]; dEE=[]; dists_chebyE=[]; ddCM=[]; ddEM=[]; ddists_chebyM=[]; ddCE=[]; ddEE=[]; ddists_chebyE=[];
% d=Afw-Af; % wzorcowe
% 
% dC=sum(abs(d));
% dE=sqrt(sum(d.^2));
% dCZ = max(abs(d));
% 
% if(wybrJakieDist)
%     switch(c) %jakieDist 1 i 3 oraz d E C CZ
%         case 1, dAll1 = [dAll1 dE];  
%         case 2, dAll2 = [dAll2 dE];  
%         case 3, dAll3 = [dAll3 dE];  
%         case 4, dAll4 = [dAll4 dE]; 
%     end   
% end
% 
% switch (jakieDist)
%     case 1
%         dCM(j,k)= dC;
%         dEM(j,k)= dE;
%         dists_chebyM(j,k) = dCZ;
%     case 2
%         dCM(j,k)= dC;
%         dEM(j,k)= dE;
%         dists_chebyM(j,k) = dCZ;
%     case 3
%         dCE(j,k)= dC;
%         dEE(j,k)= dE;
%         dists_chebyE(j,k) = dCZ;
%     case 4
%         dCE(j,k)= dC;
%         dEE(j,k)= dE;                    
%         dists_chebyE(j,k) = dCZ;
%     case 5
%         ddCM(nrG, jj(nrG), kat) = dC;
%         ddEM(nrG, jj(nrG), kat) = dE;
%         ddists_chebyM(nrG, jj(nrG), kat) = dCZ;
%     case 6                    
%         ddCE(nrG, jj(nrG), kat)= dC;
%         ddEE(nrG, jj(nrG), kat)= dE;
%         ddists_chebyE(nrG, jj(nrG), kat) = dCZ;
% end
% end

function [d] = dictionary2(lang)
    PL = 1;
    EN = 2;
    dict(PL).t1 = "Wykres zbiorczy, jakieDist = %d, %s";
    dict(PL).me = "M (Euklides)";
    dict(PL).ee = "E (Euklides)";
    dict(PL).sc = "SygKat = ";

    dict(EN).t1 = "";%Summary chart, whichDist = %d, %s";
    dict(EN).me = "";%M (Euclidean)";
    dict(EN).ee = "";%E (Euclidean)";
    dict(EN).sc = "";%SigCat = ";
    d = dict(lang);
end
