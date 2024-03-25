tic;

lingua = dictionary2(lang);

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

% offset = 0;
for(j = 1:length(v)) % grupa training
    nseg=find(fileSegNr==j);
    
    for (i = 1:length(nseg))
%         k = i;
        nrs = nseg(i);
        c = SygKat(nrs);  % index
        nrG = v(j).infoTraining; % Gestu
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
        m = [];
        if(kat-1) m = "^"; else m = "*"; end;% mięśeń marker
        if(kat-1) m = "."; else m = "."; end;% mięśeń marker
        if( find(wybrJ==j) )
            m = '*';
        end
        zakres1 = 1;
        k = SygKat(i); 
%         [j,k, jj]
        switch(k)
            case 1, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx, m);
            case 2, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx, m);
            case 3, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx, m);
            case 4, plotDistance(nrF, bf, nrs+offset*(k-1), 1, zakres1, k, [], dCx, dEx, dCzx, m);
             
            otherwise
                SygKat(i)
        end
    end
end
%% 

txCalcus(5) = lingua.t5;
txCalcus(6) = lingua.t6;
Nbf = bf; if(flagaMaxima) tx = lingua.m; else tx = lingua.e; end
for i = bf+1:Nbf+4
    if( ~length(find(i==[1,5])) ) % pomijaj rysowanie
        subplot(2,4,i); axis('tight'); hold off; 
    else
        continue;
    end
    if(i==3) title(lingua.tdw);  subtitle(" "); end
    if(i==3+4) title(lingua.tdcc); subtitle(" "); end; %subtitle( tx ); end;
%     if(i==1+4) subtitle( "Maxsima" ); end;
%     if(i==1+4) subtitle( "Energia" ); end;
  
    if i > 4 i = i-4; end; %podpisy na górze i dole
    xlabel(xlabels(mod(i,5)));    
    if (i==4) title( txCalcus(jakieDist) ); end; %if(jakieDist==6) subtitle(" ");end; end % for non overlay text plot scale
    if( i==2) subtitle(tx); end
end
toc;

function [d] = dictionary2(lang)
    PL = 1;
    EN = 2;

    dict(PL).t5 = "CC-CentrWidm";
    dict(PL).t6 = "CCE-CW";
    dict(PL).tdw = "Odległości wewnątrzgrupowe";
    dict(PL).tdcc = "Odległości od centroidu centroidów";
    dict(PL).m = "Maxima";
    dict(PL).e = "Energia";
    
    dict(EN).t5 = "CC-CentrSpectra";
    dict(EN).t6 = "CCE-CS";
    dict(EN).tdw = "Intra-group distances";
    dict(EN).tdcc = "Distances from the centroid of the centroids";
    dict(EN).m = "Maxima";
    dict(EN).e = "Energy";

    d = dict(lang);
end