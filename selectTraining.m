wybrJ = []; global fromSzU;

if(isempty(fromSzU))
    selected = "A.M.21.02.2025M58";
    fromSzU = 1; % for one cycle
else 
    selected = fromSzU;
end

for(i = 1:numel(fromSzU))
    for (j=1:length(v)) % wykluczenie z liczenia centroidu kategorii
        if(Qualisys)
            % if (v(j).infoRecord == 'L_Vastus Lateralis') wybrJ = [wybrJ j]; end
            % if (v(j).infoRecord == """22 podchwyt NORAXON ELEKTRODY""") wybrJ = [wybrJ j]; end
            % wybrJ = [wybrJ 1];
        else
            if (v(j).infoRecord == """22 pośredni NORAXON ELEKTRODY """) wybrJ = [wybrJ j]; end
            if (v(j).infoRecord == """22 podchwyt NORAXON ELEKTRODY""") wybrJ = [wybrJ j]; end
            % if (string(v(j).infoRecord) == fromSzU(i).i) wybrJ = [wybrJ j]; end
            % if (v(j).infoRecord == 'A.M.21.02.2025M58') wybrJ = [wybrJ j]; end
    
            
        end
        % if (v(j).infoRecord == """11 pośredni""") wybrJ = [wybrJ j]; end
        % if (v(j).infoRecord == """11 podchwyt""") wybrJ = [wybrJ j]; end 
    %     if (v(j).infoRecord == """35 pośredni """) wybrJ = [wybrJ j]; end
    %     if (v(j).infoRecord == """35 podchwyt""") wybrJ = [wybrJ j]; end 
    end
end
for( j = 1:length(wybrJ) )
    t = [1:length(v(wybrJ(j)).dataB)]/fpom;
%     if(exist(nrF+j,'figure') 
% close(nrF+j);
% end
    figure(nrF+j), if(~isempty(v(wybrJ(j)).dataR)) plot(t, v(wybrJ(j)).dataR); end
    hold on; plot(t, v(wybrJ(j)).dataB); subtitle(v(wybrJ(j)).infoRecord); 
    xlabel("Czas [s]"); %legend( v(wybrJ(j)).infoRDisp, v(wybrJ(j)).infoBDisp );
    ylabel(sprintf("Amplituda [%s]",string(Yunits))); axis('tight');
%     figPW("nomargin")
%     figPW(1)
end