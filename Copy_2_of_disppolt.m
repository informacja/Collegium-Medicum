 tic; n = 0; inx = 1:length(dEM(1,:));
% nag = ["między", "wew"] 
% fprintf(1,'\n\teuc. Max\tCity\tCheby\tEnerg\n')
% global n;
offset = size(Syg,1); old = [];
% for(i = 1:offset)
figure(nrF+2)
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
 
% n = i;
i=nseg(s);
    switch(SygKat(i))
        case 1         
        k = SygKat(i); plotDistance(nrF, bf, n, j, zakres1, k, Psyg, dEM, dCM, dists_chebyM)
        plotDistance(nrF, bf, n+offset, j, zakres2, k, Psyg, dEM, dCM, dists_chebyM);
        case 2
        k = SygKat(i);             
            plotDistance(nrF, bf, n, j, zakres1, k, Psyg, dEM, dCM, dists_chebyM);   
            plotDistance(nrF, bf, n+offset, j, zakres2, k, Psyg, dEM, dCM, dists_chebyM);  
         case 3
        k = SygKat(i);                
            plotDistance(nrF, bf, n, j, zakres1, k, Psyg, dEM, dCM, dists_chebyM);  
            plotDistance(nrF, bf, n+offset, j, zakres2, k, Psyg, dEM, dCM, dists_chebyM);
         case 4
            k = SygKat(i);                
            plotDistance(nrF, bf, n, j, zakres1, k, Psyg, dEM, dCM, dists_chebyM);
            plotDistance(nrF, bf, n+offset, j, zakres2, k, Psyg, dEM, dCM, dists_chebyM);   
    end

% [ i n j k pom s]
    end
%         continue;
end

figure(nrF); sgtitle("mięśnie / treningi (BR-k, BB-r Posredni-b Podchwyt-g)"); xlabels = ["Energia"; "Odległość Manhatan"; "Odległość Eukidesa"; "Odległość Chebysheva";];
for i = 1:8 subplot(2,4,i); axis('tight'); 
    if(i==2) title("Odległości wewnątrzgrupowe"); end
    if(i==2+4) title("Odległości od centroidu centroidów"); end;
    if i > 4 i = i-4; end; 
    xlabel(xlabels(mod(i,5))); 

    hold off; end
% figure(nrF+1); hold off; title("mięśnie / treningi"); xlabel("Odległość City")
% figure(nrF+2); hold off; title("mięśnie / treningi"); xlabel("Odległość Eukidesa")
% figure(nrF+3); legend(["BR/pośredni", "BB/podchwyt"]); hold off; title("mięśnie / treningi"); xlabel("Odległość Chebysheva")
% figure(nrF+1); hold off;

% figPSW

% dEsyg/liczba danych=moc

toc;