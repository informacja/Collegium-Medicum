
clear  dEM dCM dists_chebyM %Psyg
nrs = 0; nf=2; 
% global Psyg, dEM, dCM, dists_chebyM;
for(j = 1:length(v)) % grupa training
%     nk = (find(fileSegNr==j))
%     for(g = 1:2)
%         for (i=1:v(j).segLen)%length(nk))
    nseg=find(fileSegNr==j);
    for (i = 1:length(nseg)) 
        nrs = nseg(i);
        k=i;
        wyglWidma(j,k).maxAf = max(wyglWidma(j,k).Af);        
        
        if (segMio(nrs) == 1 ) kat = 1; else kat = 2; end
            if( kat == 1 )
                if v(j).infoTraining == 1
                    c=1;
                    d=CC(c,:)-wyglWidma(j,k).Af/wyglWidma(j,k).maxAf; 
                end
                if v(j).infoTraining == 2
                   c=2+1;
                   d=CC(c,:)-wyglWidma(j,k).Af/wyglWidma(j,k).maxAf; 
                end
            end
            if( fileSegMio(nrs) == txBB )
                if v(j).infoTraining == 1
                    c=2;
                    d=CC(c,:)-wyglWidma(j,k).Af/wyglWidma(j,k).maxAf; 
                end
                if v(j).infoTraining == 2
                   c=2+2;
                   d=CC(c,:)-wyglWidma(j,k).Af/wyglWidma(j,k).maxAf; 
                end
            end
            dEM(j,k)=sqrt(sum(d.^2)); 
            dCM(j,k)=sum(abs(d))/100; 
            dists_chebyM(j,k) = max(abs(d));  %uwaga przesunięcie przecinka
            Psyg(j,k) = Esyg(j,k)/SygRawLen(nrs); % unormowany
            d=CentrWidm(j, v(j).infoTraining).AfE-wyglWidma(j,k).Af'/length(wyglWidma(j,k).Af); 
            dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d))/100; dists_cheby(j,k) = max(abs(d));  %uwaga przesunięcie przecink
            d2=CentrWidm(j, v(j).infoTraining).Af2M-wyglWidma(j,k).Af2'; 
            dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2))/100; dists_cheby2(j,k) = max(abs(d2));  %2-mocy
            dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
        end % odległość w grupie
%         baza = baza+v(j).segLen;
%     end % mio
%     figure(nrF+5), subplot(111),hold on; plot(abs(d)); title("abs(d)"); 
  
end