% 0 od. centoidów J i Max,
% 1 od centroidów J i sum(Af)
% 2 od CC J  i Max
% 3 od CC i cenroidów J
nrs = 0; nf=2; %nrF = nrF+1;
% global Psyg, dEM, dCM, dists_chebyM; %ISTOTNE
tic;
ddEM = []; ddCM = []; ddists_chebyM = [];
if(1)
    jj= zeros(1,2); 
   
    for(j = 1:length(v)) % grupa training
        nseg=find(fileSegNr==j);
        for (i = 1:length(nseg))
            k = i;
            nrs = nseg(i);
            if(fileSegMio(nrs)==txBR)
                if v(j).infoTraining == 1
                    c=1; nrG=1;
                    %                 d=CC(c,:)-wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
                end
                if v(j).infoTraining == 2
                    c=2+1; nrG=2;
                    %                d=CC(c,:)-wyglWidma(j,i).Af/wyglWidma(j,i).maxAf;
                end
            end
            if(fileSegMio(nrs)==txBB)
                if v(j).infoTraining == 1
                    c=2; nrG=1;
                end
                if v(j).infoTraining == 2
                    c=2+2; nrG = 2;
                end
            end
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
                    Afw = CCE(c,:);
                    Af = CentrWidm(j,kat).AfE';
                case 5
                    Afw = CC(c,:);
                    Af = CentrWidm(j,kat).AfM';
                case 6
                    Afw = CCE(c,:);
                    Af = CentrWidm(j,kat).AfE';
            end
            d=Afw-Af; % wzorcowe

            dE=sqrt(sum(d.^2));
            dC=sum(abs(d));
            dCZ = max(abs(d));

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

% save centroids.mat CentrWidm dEM dCM dists_chebyM dEsyg Psyg dEM dCM dists_chebyM CC CCE Psr dC dE dists_cheby dEE dCE dists_chebyE dE2E dC2E dists_cheby2E dEsyg mx
toc;