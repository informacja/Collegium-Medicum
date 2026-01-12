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
            
            fontSize = 10; % for Legend

            if(gcf().Number==25) 
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width/2 hight*2]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==41) 
                set(figure(41), 'CurrentAxes', gcf().Children(1)); axis auto
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize);
                set(figure(41), 'CurrentAxes', gcf().Children(3)); axis auto
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize);
                set(figure(41), 'CurrentAxes', gcf().Children(5)); axis auto; hold on; chLen = numel(gcf().CurrentAxes.Children);
                x=gcf().CurrentAxes.Children(1).XData; y=gcf().CurrentAxes.Children(1).YData; plot(x,y,'k.');
                pos = 2; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'r.');
                pos = 3; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'b.');
                x=gcf().CurrentAxes.Children(4).XData; y=gcf().CurrentAxes.Children(4).YData; plot(x,y,'g.');
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize);
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width*1.4 hight*1.4]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==42)
                set(figure(42), 'CurrentAxes', gcf().Children(1)); axis auto
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize);
                set(figure(42), 'CurrentAxes', gcf().Children(3)); axis auto
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize);
                set(figure(42), 'CurrentAxes', gcf().Children(5)); axis auto; hold on; chLen = numel(gcf().CurrentAxes.Children);
                x=gcf().CurrentAxes.Children(1).XData; y=gcf().CurrentAxes.Children(1).YData; plot(x,y,'k.');
                pos = 2; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'r.');
                pos = 3; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'b.');
                x=gcf().CurrentAxes.Children(4).XData; y=gcf().CurrentAxes.Children(4).YData; plot(x,y,'g.');
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize);
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width*1.4 hight*1.4]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
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