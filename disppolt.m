tic; n = 0; inx = 1:length(dEM(1,:));
% nag = ["między", "wew"]
% fprintf(1,'\n\teuc. Max\tCity\tCheby\tEnerg\n')
% global n;
offset = size(Syg,1); old = [];
% for(i = 1:offset)

lingua = dictionary2(lang);

% flagaMaxima = 0;
if (flagaMaxima)
%     nrF = 3001;
    distMCtmp = distMC; distMEtmp = distME; distM_chebytmp = distM_cheby;
    dCx = dCM; dEx = dEM;  dCzx = dists_chebyM;
else % function of distribution of freq in spectrum
%     nrF = 3002;
    distMCtmp = distEC; distMEtmp = distEE; distM_chebytmp = distE_cheby;
     dCx = dCE; dEx = dEE; dCzx = dists_chebyE;
end

if (bf==0) distMEtmp(:) = 0; distMCtmp(:) = 0; distM_chebytmp(:) = 0;
    % else distMEtmp = distME; distMCtmp = distMC; distM_chebytmp = distM_cheby;
end
figure(nrF)
% dists_chebyM = dists_cheby; Psyg; dEM = dE; dCM = dC;
for(j = 1:length(v))
    %     j = fileSegNr(i);
    nseg=find(fileSegNr==j);
    for (s = 1:length(nseg))
        zakres1 = [];
        zakres2 = [];
        for (p = 1:length(nseg)) % na nowo zakresy
            %         i =s;
            nrs = nseg(p);
            %         wyglWidma(j,k).maxAf = max(wyglWidma(j,k).Af);
            if (segMio(nseg(p)) == 1 ) kat = 1;
                zakres1 = [ zakres1 nrs];
            else kat = 2; zakres2 = [ zakres2 nrs];
            end
        end
        % ze względu na Psyg
        %     zakres1 = 1:v(j).segLen;%1:max([v(:).segLen]);%;;
        %     zakres2 = v(j).segLen+1:2*v(j).segLen;%zakres1+length(zakres1);%;
        %     n = i;
        %     SygKat(i) = k;
        %     SygKat(i) % BR
        %         if(fileSegMio(nrs)==txBR)
        n = n+1;
        pom=int32(s/2); %if old == pom continue; end;
        old = pom;
        zakres1 = (pom);
        zakres2 = (pom);
%         if(4<jakieDist)
%             dEx = ddEM(1,1,1);
%             dCx = ddCM(1,1,1);
%             dCzx = ddists_chebyM(1,1,1);
%         end
        % n = i;
        i=nseg(s);
        switch(SygKat(i))
            %         distEE
            %     distEC
            %     distE_cheby
            case 1
                k = SygKat(i); %nrPoz(k) =
                plotDistance(nrF, bf, n+offset*(k-1), j, zakres1, k, Psyg, dCx, dEx, dCzx);
%                 plotDistance(nrF, bf, n+offset*(k-1), j, zakres2, k, Psyg, dEx, dCx, dCzx);
            case 2
                k = SygKat(i);
%                 plotDistance(nrF, bf, n+offset*(k-1), j, zakres1, k, Psyg, ...
%                     dEx+distMEtmp(1,1), dCx+distMCtmp(1,1), dCzx+distM_chebytmp(1,1));
                plotDistance(nrF, bf, n+offset*(k-1), j, zakres2, k, Psyg, ...
                    dCx+distMCtmp(1,1), dEx+distMEtmp(1,1), dCzx+distM_chebytmp(1,1));
            case 3 %+2,1
                k = SygKat(i);
                plotDistance(nrF, bf, n+offset*(k-1), j, zakres1, k, Psyg, ...
                    dCx+distMCtmp(2,1), dEx+distMEtmp(2,1), dCzx+distM_chebytmp(2,1));
%                 plotDistance(nrF, bf, n+offset*(k-1), j, zakres2, k, Psyg, ...
%                     dEx+distMEtmp(2,1), dCx+distMCtmp(2,1), dCzx+distM_chebytmp(2,1));
            case 4% 2,2+1,2
                k = SygKat(i);
%                 plotDistance(nrF, bf, n+offset*(k-1), j, zakres1, k, Psyg, ...
%                     dEx+distMEtmp(2,2)+distMEtmp(1,2), dCx+distMCtmp(2,2)+distMCtmp(1,2), dCzx+distM_chebytmp(2,2)+distM_chebytmp(1,2));
                 plotDistance(nrF, bf, n+offset*(k-1), j, zakres2, k, Psyg, ...
                    dCx+distMCtmp(2,2)+distMCtmp(1,2), dEx+distMEtmp(2,2)+distMEtmp(1,2), dCzx+distM_chebytmp(2,2)+distM_chebytmp(1,2));
        end
        % [ i n j k pom s]
    end
    %         continue;
end


txCalcus(1) = lingua.t1;
txCalcus(2) = lingua.t2;
txCalcus(3) = lingua.t3;
txCalcus(4) = lingua.t4;
figure(nrF); sgtitle(lingua.sgt); xlabels = lingua.xlabels;
if (flagaMaxima) txDist = lingua.m; else txDist = lingua.e; end
Nbf = bf;
for i = bf+1:Nbf+4
    subplot(2,4,i); axis('tight');
    if(i==2) title(lingua.tdw);  subtitle(txDist); end
    if(i==2+4) title(lingua.tdcc); end;
    if i > 4 i = i-4; end;
    xlabel(xlabels(mod(i,5)));
    hold off; 
    if (i==4) subtitle( txCalcus(jakieDist) ); end
end
% figure(nrF+1); hold off; title("mięśnie / treningi"); xlabel("Odległość City")
% figure(nrF+2); hold off; title("mięśnie / treningi"); xlabel("Odległość Eukidesa")
% figure(nrF+3); legend(["BR/pośredni", "BB/podchwyt"]); hold off; title("mięśnie / treningi"); xlabel("Odległość Chebysheva")
% figure(nrF+1); hold off;

% figPSW

% dEsyg/liczba danych=moc

toc;

function [d] = dictionary2(lang)
    PL = 1;
    EN = 2;

    dict(PL).xlabels = ["Moc"; "Odległość Manhattan"; "Odległość Euklidesa"; "Odległość Czebyszewa";]; 
    dict(PL).m = "Maxima";
    dict(PL).e = "Energia";    
    dict(PL).t1 = "CentrWidm-wyglWidma/maxAf";
    dict(PL).t2 = "CC-wyglWidma/maxAf";
    dict(PL).t3 = "CentrWidm-wyglWidma/Ps";
    dict(PL).t4 = "CCE-wyglWidma";
    dict(PL).sgt = "mięśnie / treningi (PS/BR-k, PS/BB-r PC/BR-b PC/BB-g)";
    dict(PL).tdw = "Odległości wewnątrzgrupowe";
    dict(PL).tdcc = "Odległości od centroidu centroidów";

    dict(EN).xlabels = ["Power"; "Manhattan Distance"; "Euclidean Distance"; "Chebyshev Distance";];    
    dict(EN).m = "Maxima";
    dict(EN).e = "Energy"; 
    dict(EN).t1 = "CentrSpectra-SmoothSpectra/maxAf";
    dict(EN).t2 = "CC-SmoothSpectra/maxAf";
    dict(EN).t3 = "CentrSpectra-SmoothSpectra/Ps";
    dict(EN).t4 = "CCE-SmoothSpectra";
    dict(EN).sgt = "muscles / trainings (IM/BR-k, IM/BB-r SP/BR-b SP/BB-g)";
    dict(EN).tdw = "Intra-group distances";
    dict(EN).tdcc = "Distances from the centroid of the centroids";

    d = dict(lang);
end
