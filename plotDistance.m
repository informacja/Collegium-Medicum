function plotDistance(nrF, bf, n, j, zakres, k, Psyg, dEM, dCM, dists_chebyM,m)
    if(nargin < 11) m = []; end 
    if(m=='*') jM=2; else jM=1; end
            if isempty(m) kol = [ "k."; "r."; "b."; "g."];
            else kol = [ strcat('k',m); strcat('r',m); strcat('b',m); strcat('g',m); ];
            end
        %     figure(nrF); 
      for(jj=1:jM)

        if( ~isempty(Psyg) )
                subplot(2,4,1+bf); hold on; plot(n, Psyg(j, zakres), kol(k,:));
            end
             subplot(2,4,2+bf); hold on; plot(n, dEM(j, zakres), kol(k,:));
            subplot(2,4,3+bf); hold on; plot(n, dCM(j, zakres), kol(k,:));
            subplot(2,4,4+bf); hold on; plot(n, dists_chebyM(j, zakres), kol(k,:));
            kol(k,:)="ko"; % dla kolejnego obiegu pętli
     end
        
%     dists_chebyM(j, zakres)
%     subplot(2,4,4+bf); hold on; plot(n, dists_chebyM(j, zakres), kol(k)); 
end

