tic;
% if (flagaMaxima) czy robić dodawanie?
%     %     nrF = 3001;
%     distMCtmp = distMC; distMEtmp = distME; distM_chebytmp = distM_cheby;
%     dCx = dCM; dEx = dEM;  dCzx = dists_chebyM;
% else % function of distribution of freq in spectrum
%     %     nrF = 3002;
%     distMCtmp = distEC; distMEtmp = distEE; distM_chebytmp = distE_cheby;
%     dCx = dCE; dEx = dEE; dCzx = dists_chebyE;
% end
z=[];
jj= zeros(1,2);

for(j = 1:length(v)) % grupa training
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg))
%         k = i;
        nrs = nseg(i);
        c = SygKat(nrs);  % index
        nrG = v(j).infoTraining; % Gestu
%         if(fileSegMio(nrs)==txBR)
%             if v(j).infoTraining == 1 c=1; nrG=1; end
%             if v(j).infoTraining == 2 c=2+1; nrG=2; end
%         end
%         if(fileSegMio(nrs)==txBB)
%             if v(j).infoTraining == 1 c=2; nrG=1; end
%             if v(j).infoTraining == 2 c=2+2; nrG = 2; end
%         end
        kat = segment(nrs).miesien;
        jj(nrG) = jj(nrG)+1;
        switch (jakieDist)            
            case 5
                dCx = ddCM(nrG, jj(nrG), kat);
                dEx = ddEM(nrG, jj(nrG), kat);
                dCzx= ddists_chebyM(nrG, jj(nrG), kat);
            case 6
                dCx = ddCE(nrG, jj(nrG), kat);
                dEx = ddEE(nrG, jj(nrG), kat);
                dCzx= ddists_chebyE(nrG, jj(nrG), kat);
        end
%         z=[z;nrG, jj(nrG), kat];
        zakres1 = 1;
        k = SygKat(i); 
%         [j,k, jj]
        switch(k)
            case 1, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx);
            case 2, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx);
            case 3, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx);
            case 4, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx);
             
            otherwise
                SygKat(i)
        end
    end
end

txCalcus(5) = "CC-CentrWidm";
txCalcus(6) = "CCE-CentrWidm";
Nbf = bf; if(flagaMaxima) tx = "Maxsima"; else tx = "Energia"; end
for i = bf+1:Nbf+4
    if( ~length(find(i==[1,5])) ) % pomijaj rysowanie
        subplot(2,4,i); axis('tight'); hold off; 
    else
        continue;
    end
    if(i==2) title("Odległości wewnątrzgrupowe"); end
    if(i==2+4) title("Odległości od centroidu centroidów"); end; %subtitle( tx ); end;
%     if(i==1+4) subtitle( "Maxsima" ); end;
%     if(i==1+4) subtitle( "Energia" ); end;
    if i > 4 i = i-4; end; %podpisy na górze i dole
    xlabel(xlabels(mod(i,5)));    
    if (i==4) subtitle( txCalcus(jakieDist) ); end
end
toc;