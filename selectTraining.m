wybrJ = [];
for (j=1:length(v)) % wykluczenie z liczenia centroidu kategorii
    if (v(j).infoRecord == """22 pośredni NORAXON ELEKTRODY """) wybrJ = [wybrJ j]; end
    if (v(j).infoRecord == """22 podchwyt NORAXON ELEKTRODY""") wybrJ = [wybrJ j]; end
    % if (v(j).infoRecord == """11 pośredni""") wybrJ = [wybrJ j]; end
    % if (v(j).infoRecord == """11 podchwyt""") wybrJ = [wybrJ j]; end 
%     if (v(j).infoRecord == """35 pośredni """) wybrJ = [wybrJ j]; end
%     if (v(j).infoRecord == """35 podchwyt""") wybrJ = [wybrJ j]; end 
end

for( j = 1:length(wybrJ) )
    t = [1:length(v(wybrJ(j)).dataR)]/fpom;
%     if(exist(nrF+j,'figure') 
% close(nrF+j);
% end
    figure(nrF+j), plot(t, v(wybrJ(j)).dataR); hold on; plot(t, v(wybrJ(j)).dataB); title(v(wybrJ(j)).infoRecord); 
    xlabel("Czas [s]"); legend( v(wybrJ(j)).infoRDisp, v(wybrJ(j)).infoBDisp );
    ylabel(sprintf("Amplituda [%s]",string(Yunits))); axis('tight');
%     figPW("nomargin")
%     figPW(1)
end