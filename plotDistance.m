function plotDistance(nrF, bf, n, j, zakres, k, Psyg, dEM, dCM, dists_chebyM, m)
    if(nargin < 11) m = []; end 
    global colList;
    if(m=='*') jM=2; else jM=1; end
            if isempty(m) kol = [ "k."; "r."; "b."; "g."];
            else kol = [ strcat('k',m); strcat('r',m); strcat('b',m); strcat('g',m); ];
            end
        %     figure(nrF);
      global cntCol;
      global lastCirclePos;
      for(jj=1:jM)

        if( ~isempty(Psyg) )
            subplot(2,4,1+bf); hold on; plot(n, Psyg(j, zakres), kol(k,:));
        end
            subplot(2,4,2+bf); hold on; plot(n, dEM(j, zakres), kol(k,:));
            subplot(2,4,3+bf); hold on; plot(n, dCM(j, zakres), kol(k,:));
            subplot(2,4,4+bf); hold on; h = plot(n, dists_chebyM(j, zakres), kol(k,:));
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % if(
            kol(k,:)=strcat(colList(cntCol),'o');%"ko"; % dla kolejnego obiegu pętli, rysuj kółko
      end
      if(bf == 4 && k == 4)
          if(jM == 2)
              % kkk = char(kol(k,:));
              % if(length(kkk) > 1)
                  % if( kkk(2) == 'o')
                    lastCirclePos = [ lastCirclePos; [n, dists_chebyM(j, zakres)]];
                    zakres
                    % [nrF, bf, n, j, zakres, k, Psyg, dEM, dCM, dists_chebyM, m]
                  % end
              % end
          end
      end
        
%     dists_chebyM(j, zakres)
%     subplot(2,4,4+bf); hold on; plot(n, dists_chebyM(j, zakres), kol(k)); 
end

