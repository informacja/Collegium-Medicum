
clear Psyg dEM dCM dists_chebyM
nrs = 0; nf=2; 
% global Psyg, dEM, dCM, dists_chebyM;
for(j = 1:length(v)) % grupa training
%     nk = (find(fileSegNr==j))
%     for(g = 1:2)
%         for (i=1:v(j).segLen)%length(nk))
    for (k = 1:length(find(fileSegNr==j))) 
            nrs = nrs + 1;
            d = [];
            if( fileSegMio(nrs) == txBR )
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
            d=CentrWidm(j).AfE-wyglWidma(j,k).Af'/Psyg(j,k); 
            dE(j,k)=sqrt(sum(d.^2)); dC(j,k)=sum(abs(d))/100; dists_cheby(j,k) = max(abs(d),[],1);  %uwaga przesunięcie przecink
            d2=CentrWidm(j).Af2M-wyglWidma(j,k).Af2'; 
            dE2(j,k)=sqrt(sum(d2.^2)); dC2(j,k)=sum(abs(d2))/100; dists_cheby2(j,k) = max(abs(d2),[],1);  %2-mocy
            dEsyg(j,k)=abs(Psyg(j,k)-Psr(j));
        end % odległość w grupie
%         baza = baza+v(j).segLen;
%     end % mio
%     figure(nrF+5), subplot(111),hold on; plot(abs(d)); title("abs(d)"); 
  
end