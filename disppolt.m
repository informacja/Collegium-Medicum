nrs = 0;
% nag = ["między", "wew"] 
fprintf(1,'\n\teuc. Max\tCity\tCheby\tEnerg\n')
for(j = 1:length(v)) % grupa
%     fprintf(1,'\ngr.%-2d',j)
    %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
%     for (k = 1:length(find(fileSegNr==j))) 
%         nrs = nrs + 1;
%         fprintf(1,';  %6.3f %.3f %.3f %.3f',dEM(j,k),dCM(j,k),dists_chebyM(j,k),dEsyg(j,k)/SygRawLen(nrs)); %] [dE(2,1);dC(2,1);dists_cheby(2,1);;dEsyg(2,1)
%     end
%     figure(nrF); plot(dEM(j,:),'k.'); hold on; plot(dCM(j,:),'r*'); hold off; nrF = nrF+1;
%     figure(nrF); stem(dEM(j,:)); hold on; stem(dCM(j,:)); hold off; nrF = nrF+1;
% %     figure(nrF); plot(dEM(j,:),'k*'); hold on; 
% %     plot(dCM(j,:),'r*');  
% %     plot(sqrt(dists_chebyM(j,:)),'g*'); 
% %     plot(sqrt(dEM(j,:)),'c*'); hold off; 
% %     legend('dEM','dCM','dists chebyM', 'dEM');  nrF = nrF+1;
% %     figure(nrF); plot(xcorr(dEM(j,:))); title("Cross corelation"); nrF = nrF+1; %xlabel()

figure(nrF); hold on; plot(xcorr(v(j).dataB,v(j).dataR)); title(v(j).infoRecord); xlabel('xcorr(v(j).dataB,v(j).dataR)'); %nrF = nrF+1; 
% figure(nrF); plot(xcorr(v(j).dataB(1:1000),v(j).dataR(1:1000))); title("Cross corelation"); nrF = nrF+1; 

%     dEM(j,:),dataB
% 
%     figure(nrF); plot(dEsyg(j,:),'g^');  hold off; nrF = nrF+1;
end


% dEsyg/liczba danych=moc