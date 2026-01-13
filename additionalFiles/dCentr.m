% p = 8;
% a = v(p).dataB;
% b = v(p).dataR; korelacja nie jest informatywna
% a = segment(1).data;
% b = segment(2).data;
% N = length(a);
%     n = 0:N-1;
% 
%     [r,lags] = xcorr(a',b', 'normalized');
%     n0 = find(lags==0);
%     figure(1), plot(lags(N:2*N-1),r(N:2*N-1)); axis('tight');%stem(lags,r)
% figure(2), plot(a); hold on; plot(b); hold off;
%     return
tic; 
txTraining = [txPr txPc];
nrs = 0; dMM=[]; d=[]; dists_cheby=[]; k = 0;
i1=[1 3 1 2]; i2=[2 4 3 4]; nk=1;
for i = 1:4 %2:3  % 12,34, 13,24
    if(i==3) k=0; end,
    k = k +1; 
    if(i>2) nk=2; end 
    for(nband=1:twoOrOne)
        if(nband-1)
            d=CC(i1(i),nflim+1:end)-CC(i2(i),nflim+1:end);
        else
            d=CC(i1(i),1:nflim)-CC(i2(i),1:nflim);
        end
        % d=CC(i,:)-CC(i+1,:); % różnica posr(R-B), roz. podchwyt(R-B)
        Sb(nband).distME(nk,k)=sqrt(sum(d.^2));
        Sb(nband).distMC(nk,k)=sum(abs(d));
        Sb(nband).distM_cheby(nk,k) = max(abs(d));
        if(nband-1)
            d=CCE(i1(i),nflim+1:end)-CCE(i2(i),nflim+1:end);
        else
            d=CCE(i1(i),1:nflim)-CCE(i2(i),1:nflim);
        end
        % d=CCE(i,:)-CCE(i+1,:); % różnica posr(R-B), roz. podchwyt(R-B)
        Sb(nband).distEE(nk,k)=sqrt(sum(d.^2));
        Sb(nband).distEC(nk,k)=sum(abs(d));
        Sb(nband).distE_cheby(nk,k) = max(abs(d));
        if(i<3)
            fprintf(1, "%s (BR-BB) distM(1,%d)\tM=%f, E=%f, C=%f\n", txTraining(k),k, ...
                Sb(nband).distME(1,k), Sb(nband).distMC(1,k), Sb(nband).distM_cheby(1,k))
            fprintf(1, "%s (BR-BB) distE(1,%d)\tM=%f, E=%f, C=%f (Energy)\n", txTraining(k),k, ...
                Sb(nband).distEE(1,k), Sb(nband).distEC(1,k), Sb(nband).distE_cheby(1,k))
        end
    end
% end
% k=0;
% txTraining = [txBR txBB];
% for i = 1:2 % 13,24
%     k = k +1;
%     if(nband-1)
%         d=CC(i,nflim+1:end)-CC(i+2,nflim+1:end);
%     else
%         d=CC(i,1:nflim)-CC(i+2,1:nflim);
%     end
%     Sb(nband).distME(2,k)=sqrt(sum(d.^2)); 
%     Sb(nband).distMC(2,k)=sum(abs(d)); 
%     Sb(nband).distM_cheby(2,k) = max(abs(d)); 
%     if(nband-1)
%         d=CCE(i,nflim+1:end)-CCE(i+1,nflim+1:end);
%     else
%         d=CCE(i,1:nflim)-CCE(i+1,1:nflim);
%     end    
%     Sb(nband).distEE(2,k)=sqrt(sum(d.^2)); 
%     Sb(nband).distEC(2,k)=sum(abs(d)); 
%     Sb(nband).distE_cheby(2,k) = max(abs(d)); 
    if(i>2)
        fprintf(1, "%s (NT-SP) distM(2,%d)\tM=%f, E=%f, C=%f\n", txTraining(k),k, ...
            Sb(nband).distME(2,k), Sb(nband).distMC(2,k), Sb(nband).distM_cheby(2,k))
        fprintf(1, "%s (NT-SP) distE(2,%d)\tM=%f, E=%f, C=%f (Energy)\n", txTraining(k),k, ...
            Sb(nband).distEE(2,k), Sb(nband).distEC(2,k), Sb(nband).distE_cheby(2,k))
    end
end
 fprintf(1, "STD %s \tM=%f, E=%f, C=%f\n", txTraining(k), std(Sb(nband).distME(:)), std(Sb(nband).distMC(:)), std(Sb(nband).distM_cheby(:)))
 fprintf(1, "STD %s \tM=%f, E=%f, C=%f Energy\n", txTraining(k), std(Sb(nband).distEE(:,k)), std(Sb(nband).distEC(:)), std(Sb(nband).distE_cheby(:)))


% distM/E... (1,1:2) = odległość R-B dla S i C
% distM/E... (2,1:2) = odległość S-C dla R i B
% nrs = 0; dMM=[]; d=[]; dists_cheby=[];
% for i = 1:2
%     nrs = nrs +1;
%     d=CC(i,:)-CC(i+2,:); % różnica posr(R-B), roz. podchwyt(R-B)
%     dE(nrs)=sqrt(sum(d.^2)); 
%     dC(nrs)=sum(abs(d))/100; 
%     dists_cheby(nrs) = max(abs(d));  
%     dMM(nrs,:) = dE(nrs);
%     fprintf(1, "%s (BR-BB) dist. M=%f, E=%f, C=%f\n", txTraining(nrs),dE(nrs), dC(nrs), dists_cheby(nrs))
% end
%----------------------d----------------------------------------------------          
toc;
