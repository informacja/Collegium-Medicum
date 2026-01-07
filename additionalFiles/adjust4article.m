function adjust4article(s)
    if(isfield(s,"figNrGridSize")) figNrGridSize = s.figNrGridSize; end
    
    if(~isempty(ismember(findall(0,'type','figure'),groot))) % jeśli są otwarte figury, to je eksportuj     
        F = findobj('Type', 'figure'); postfix = char('a'-1);
        for(i=1:numel(F))
            if(exist("figNrList", "var") && isempty(find(figNrList==F(i)))) continue; end
            figure(F(i))   

            set(gcf,'WindowStyle','normal') % undock
            width = 13; % inches
            hight =  6.5000;
            
            if(gcf().Number==25) 
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width/2 hight*2]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==41) 
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==42) 
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==4003) 
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end

            return
            if(exist("figNrGridSize","var"))
                rowNr = find(gcf().Number==cell2mat(figNrGridSize(:,1)));
                if(~isempty(rowNr))
                    dim = figNrGridSize(rowNr,2:3);
                    set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                    set(gcf,'Position', [0 0 width*2*dim{1} hight*2*dim{2}]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
                end
            end
        end
    end
end