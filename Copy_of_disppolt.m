 tic; n = 0; inx = 1:length(dEM(1,:));
% nag = ["między", "wew"] 
% fprintf(1,'\n\teuc. Max\tCity\tCheby\tEnerg\n')
for(i = 1:4)
   nrs = 0; 
   %n=0;%separacja
%    if i == 1 || i == 4
%             n = 0;
%         else% if i == 3 || i == 2
%             n = length(v);
%         end
    for(j = 1:length(v)) % grupa

            n = n+1;
    %     fprintf(1,'\ngr.%-2d',j)
        %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
    %     for (k = 1:length(find(fileSegNr==j))) 
            nrs = nrs + 1;
    %         fprintf(1,';  %6.3f %.3f %.3f %.3f',dEM(j,k),dCM(j,k),dists_chebyM(j,k),dEsyg(j,k)/SygRawLen(nrs)); %] [dE(2,1);dC(2,1);dists_cheby(2,1);;dEsyg(2,1)
    %     end
    %     figure(nrF); plot(dEM(j,:),'k.'); hold on; plot(dCM(j,:),'r*'); hold off; nrF = nrF+1;
    %     figure(nrF); stem(dEM(j,:)); hold on; stem(dCM(j,:)); hold off; nrF = nrF+1;
    zakres1 = 1:v(j).segLen;%1:max([v(:).segLen]);%;;
    zakres2 = v(j).segLen+1:2*v(j).segLen;%zakres1+length(zakres1);%;

    if i == 1 % BR 
        if(fileSegMio(nrs)==txBR)
            k = 1;
            if v(j).infoTraining == 1
                plotDistance(nrF, bf, n, j, zakres1, k, Psyg, dEM, dCM, dists_chebyM)

%         figure(30), plot(j,dists_chebyM(j,nrs),'k.'); hold on;
            end
            if v(j).infoTraining == 2
                plotDistance(nrF, bf, n, j, zakres2, k, Psyg, dEM, dCM, dists_chebyM);
            end
        end
    end
    if i == 2 %BB
        if(fileSegMio(n)==txBB)
            k = 2;
            if v(j).infoTraining == 1                
                plotDistance(nrF, bf, n, j, zakres1, k, Psyg, dEM, dCM, dists_chebyM);   
            end
            if v(j).infoTraining == 2                
                plotDistance(nrF, bf, n, j, zakres2, k, Psyg, dEM, dCM, dists_chebyM);  
            end
        end
    end

    if i == 3 % podchwyt
        if v(j).infoTraining == 1
            k = 1;
            if(fileSegMio(nrs)==txBR)                
                plotDistance(nrF, bf, n, j, zakres1, i, Psyg, dEM, dCM, dists_chebyM);  
            end
            if(fileSegMio(nrs)==txBB)                
                plotDistance(nrF, bf, n, j, zakres2, i, Psyg, dEM, dCM, dists_chebyM);
            end
        end
    end

     if i == 4 % posredni
        if v(j).infoTraining == 2
            k = 2;
            if(fileSegMio(nrs)==txBR)                
                plotDistance(nrF, bf, n, j, zakres1, i, Psyg, dEM, dCM, dists_chebyM);
            end
            if(fileSegMio(nrs)==txBB)                 
                plotDistance(nrF, bf, n, j, zakres2, i, Psyg, dEM, dCM, dists_chebyM);   
            end
        end
      end

    if i == 1
%         figure(nrF+1), plot(inx, dEM(j,:)), hold on; inx = inx + length(dEM(j,:));

    end

    
        %     if i == 3 % branchioradialis
%         if v(j). == 2
%                 figure(nrF); hold on; plot(nrs, dEM(j,:),'k.');   title("Euclid Max T1")
%         end
%     end
%     if i == 4 % biceps
%         if v(j).infoTraining == 2
%                 figure(nrF); hold on; plot(nrs, dEM(j,:),'k.');   title("Euclid Max T1")
%         end
%     end
%     else %if v(j).infoTraining == 2
%             figure(nrF+1); hold on; plot(dEM(j,:),'k.');  title("Euclid Max T2")
%         end
    %     plot(dCM(j,:),'r*');  
    %     plot(sqrt(dists_chebyM(j,:)),'g*'); 
    %     plot(sqrt(dEM(j,:)),'c*'); hold off; 
    %     legend('dEM','dCM','dists chebyM', 'dEM');  nrF = nrF+1;
%         figure(nrF+3); plot(xcorr(dEM(j,:))); title("Cross corelation");% nrF = nrF+1; %xlabel()
    
    % figure(nrF); hold on; plot(xcorr(v(j).dataB,v(j).dataR)); title(v(j).infoRecord); xlabel('xcorr(v(j).dataB,v(j).dataR)'); %nrF = nrF+1; 
    % figure(nrF); plot(xcorr(v(j).dataB(1:1000),v(j).dataR(1:1000))); title("Cross corelation"); nrF = nrF+1; 
    
    %     dEM(j,:),dataB
    % 
    %     figure(nrF); plot(dEsyg(j,:),'g^');  hold off; nrF = nrF+1;
    end
end

figure(nrF); sgtitle("mięśnie / treningi (BR-k, BB-r Posredni-b Podchwyt-g)"); xlabels = ["Energia"; "Odległość Manhatan"; "Odległość Eukidesa"; "Odległość Chebysheva";];
for i = 1:8 subplot(2,4,i); axis('tight'); 
    if(i==2) title("Odległości wewnątrzgrupowe"); end
    if(i==2+4) title("Odległości od centroidów"); end;
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