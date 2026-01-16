function adjust4article(s)
    if(isfield(s,"figNrGridSize")) figNrGridSize = s.figNrGridSize; end
    
    if(~isempty(ismember(findall(0,'type','figure'),groot))) % jeśli są otwarte figury, to je eksportuj     
        F = findobj('Type', 'figure'); postfix = char('a'-1);
        for(i=1:numel(F))
            if(exist("figNrList", "var") && isempty(find(figNrList==F(i)))) continue; end
            figure(F(i))   
            set(0, 'currentfigure', F(i));

            set(gcf,'WindowStyle','normal') % undock
            width = 13; % inches
            hight =  6.5000;
            
            fontSize = 10; % for Legend

            if(gcf().Number==25) 
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width/2 hight*2]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==41) m=1.4; delete(findobj(gcf, 'type', 'Legend'));
                set(figure(41), 'CurrentAxes', gcf().Children(1)); axis tight
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(figure(41), 'CurrentAxes', gcf().Children(3)); axis tight
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(figure(41), 'CurrentAxes', gcf().Children(5)); axis tight; hold on; chLen = numel(gcf().CurrentAxes.Children);
                x=gcf().CurrentAxes.Children(1).XData; y=gcf().CurrentAxes.Children(1).YData; plot(x,y,'k.');
                pos = 2; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'r.');
                pos = 3; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'b.');
                x=gcf().CurrentAxes.Children(4).XData; y=gcf().CurrentAxes.Children(4).YData; plot(x,y,'g.');
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width*1.4 hight*1.4]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==42)  m = 1.4; delete(findobj(gcf, 'type', 'Legend'));
                set(figure(42), 'CurrentAxes', gcf().Children(1)); axis tight
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(figure(42), 'CurrentAxes', gcf().Children(3)); axis tight
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(figure(42), 'CurrentAxes', gcf().Children(5)); axis tight; hold on; chLen = numel(gcf().CurrentAxes.Children);
                x=gcf().CurrentAxes.Children(1).XData; y=gcf().CurrentAxes.Children(1).YData; plot(x,y,'k.');
                pos = 2; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'r.');
                pos = 3; x=gcf().CurrentAxes.Children(pos).XData; y=gcf().CurrentAxes.Children(pos).YData; plot(x,y,'b.');
                x=gcf().CurrentAxes.Children(4).XData; y=gcf().CurrentAxes.Children(4).YData; plot(x,y,'g.');
                legend(["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northwest','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(gcf,'Units','inches');                        % jednostka wyiarowania okna
                set(gcf,'Position', [0 0 width*1.4 hight*1.4]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end
            if(gcf().Number==509) 
                delete(findobj(gcf, 'type', 'Legend'));
                set(figure(509), 'CurrentAxes', gcf().Children(2)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(figure(509), 'CurrentAxes', gcf().Children(4)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(figure(509), 'CurrentAxes', gcf().Children(6)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(figure(509), 'CurrentAxes', gcf().Children(8)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
           
            end
            if(gcf().Number==512) delete(findobj(gcf, 'type', 'Legend'));
                set(figure(512), 'CurrentAxes', gcf().Children(2)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(figure(512), 'CurrentAxes', gcf().Children(4)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(figure(512), 'CurrentAxes', gcf().Children(6)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(figure(512), 'CurrentAxes', gcf().Children(8)); legend(gcf().CurrentAxes.Children( [end-1 end-2 end-3] ), ["Gaussian" "Doubly exponential" "Maxwell–Boltzmann"], 'Location','northeast','FontSize',fontSize);
                set(gcf,'Units','inches');                        % jednostka wymiarowania okna
                set(gcf,'Position', [0 0 width hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
           
            end
            if(gcf().Number==4003) m = 1; delete(findobj(gcf, 'type', 'Legend'));
                set(figure(4003), 'CurrentAxes', gcf().Children(6)); legend(gcf().CurrentAxes.Children( legendHelper(gcf().CurrentAxes) ), ["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','bestoutside','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                % set(figure(4003), 'CurrentAxes', gcf().Children(3)); legend(gcf().CurrentAxes.Children( legendHelper(gcf().CurrentAxes) ), ["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northeast','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                % set(figure(4003), 'CurrentAxes', gcf().Children(5)); legend(gcf().CurrentAxes.Children( legendHelper(gcf().CurrentAxes) ), ["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northeast','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                % set(figure(4003), 'CurrentAxes', gcf().Children(7)); legend(gcf().CurrentAxes.Children( legendHelper(gcf().CurrentAxes) ), ["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northeast','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                % set(figure(4003), 'CurrentAxes', gcf().Children(9)); legend(gcf().CurrentAxes.Children( legendHelper(gcf().CurrentAxes) ), ["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northeast','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                % set(figure(4003), 'CurrentAxes', gcf().Children(11));legend(gcf().CurrentAxes.Children( legendHelper(gcf().CurrentAxes) ), ["NT/BR";"NT/BB";"SP/BR";"SP/BB"],'Location','northeast','FontSize',fontSize); ylim(gca, [0 gca().YLim(2)*m]);
                set(gcf,'Units','inches');                        % jednostka wymiarowania oknanortheast
                set(gcf,'Position', [0 0 width hight]); % wymiary okna: (x,y,dx,dy), (x,y)- lewy dolny
            end

            continue;
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

function [nr]=legendHelper(axis)
blackClr = [0 0 0];
redClr  = [1 0 0];
blueClr = [0 0 1];
greenClr = [0 1 0];
k = []; r = []; b = []; g = [];

i=1;
    while(isempty(k)||isempty(r)||isempty(b)||isempty(g))
        if(length(find(axis.Children(i).Color==blackClr))==3) k = i; end
        if(length(find(axis.Children(i).Color==redClr))==3) r = i; end
        if(length(find(axis.Children(i).Color==blueClr))==3) b = i; end
        if(length(find(axis.Children(i).Color==greenClr))==3) g = i; end
        i=i+1;
        if(i > length(axis.Children)) break; end
    end
    nr=[k r b g];
% delete(findobj(gcf, 'type', 'Legend'));
end