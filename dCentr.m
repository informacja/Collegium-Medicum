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
txTraining = ["Pośredni" "Podchwyt"];
nrs = 0; dMM=[]; d=[]; dists_cheby=[];
for i = 1:2:3
    nrs = nrs +1;
    d=CC(i,:)-CC(i+1,:); % posr, podchyt
    dE(nrs)=sqrt(sum(d.^2)); 
    dC(nrs)=sum(abs(d))/100; 
    dists_cheby(nrs) = max(abs(d));  
    dMM(nrs,:) = dE(nrs);
    fprintf(1, "%s (BR-BB) dist. M=%f, E=%f, C=%f\n", txTraining(nrs),dE(nrs), dC(nrs), dists_cheby(nrs))
end
%----------------------d----------------------------------------------------          
toc;
