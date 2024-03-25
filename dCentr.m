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
for i = 1:2:3
    k = k +1;
    d=CC(i,:)-CC(i+1,:); % różnica posr(R-B), roz. podchwyt(R-B)
    distME(1,k)=sqrt(sum(d.^2)); 
    distMC(1,k)=sum(abs(d)); 
    distM_cheby(1,k) = max(abs(d)); 
    d=CCE(i,:)-CCE(i+1,:); % różnica posr(R-B), roz. podchwyt(R-B)
    distEE(1,k)=sqrt(sum(d.^2)); 
    distEC(1,k)=sum(abs(d)); 
    distE_cheby(1,k) = max(abs(d)); 
    fprintf(1, "%s (BR-BB) distM(1,%d) M=%f, E=%f, C=%f\n", txTraining(k),k, ...
        distME(1,k), distMC(1,k), distM_cheby(1,k))
    fprintf(1, "%s (BR-BB) distE(1,%d) M=%f, E=%f, C=%f\n", txTraining(k),k,distEE(1,k), distEC(1,k), distE_cheby(1,k))
end
k=0;
for i = 1:2
    k = k +1;
    d=CC(i,:)-CC(i+2,:); % różnica posr(R-B), roz. podchwyt(R-B)
    distME(2,k)=sqrt(sum(d.^2)); 
    distMC(2,k)=sum(abs(d)); 
    distM_cheby(2,k) = max(abs(d)); 
    d=CCE(i,:)-CCE(i+2,:); % różnica posr(R-B), roz. podchwyt(R-B)
    distEE(2,k)=sqrt(sum(d.^2)); 
    distEC(2,k)=sum(abs(d)); 
    distE_cheby(2,k) = max(abs(d)); 
    fprintf(1, "%s (PS-PC) distM(2,%d) M=%f, E=%f, C=%f\n", txTraining(k),k, ...
        distME(2,k), distMC(2,k), distM_cheby(2,k))
    fprintf(1, "%s (PS-PC) distE(2,%d) M=%f, E=%f, C=%f\n", txTraining(k),k, distEE(2,k), distEC(2,k), distE_cheby(2,k))
end
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
