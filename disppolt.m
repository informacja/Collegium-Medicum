nrs = 0;
% nag = ["między", "wew"] 
fprintf(1,'\n\teuc. Max\tCity\tCheby\tEnerg')
for(j = 1:length(v)) % grupa
%     fprintf(1,'\ngr.%-2d',j)
    %dE(j)=0; dC(j)=0; % norma Euklidesowa, City (Manhatan)
%     for (k = 1:length(find(fileSegNr==j))) 
%         nrs = nrs + 1;
%         fprintf(1,';  %6.3f %.3f %.3f %.3f',dEM(j,k),dCM(j,k),dists_chebyM(j,k),dEsyg(j,k)/SygRawLen(nrs)); %] [dE(2,1);dC(2,1);dists_cheby(2,1);;dEsyg(2,1)
%     end
    figure(500); plot(dEM(j,:),'k.'); hold on; plot(dCM(j,:),'r*'); hold off;
    figure(501); stem(dEM(j,:)); hold on; stem(dCM(j,:)); hold off;
    figure(502); plot(dEM(j,:),'k*'); hold on; 
    plot(dCM(j,:),'r*');  
    plot(sqrt(dists_chebyM(j,:)),'g*'); 
    plot(sqrt(dEsyg(j,:)),'c*'); hold off; legend;  
    figure(503); plot(dEsyg(j,:),'g^');  hold off; 

    figure(502); plot(dEM(j,:),'k*'); hold on; 
    plot(dCM(j,:),'r*');  
    plot(sqrt(dists_chebyM(j,:)),'g*'); 
    plot(sqrt(dEsyg(j,:)),'c*'); hold off; legend;
end

xcorr(dEM(j,:))

% dEsyg/liczba danych=moc