% https://stackoverflow.com/questions/26084069/how-to-produce-unicode-characters-with-matlab-latex-interpreter
% rmdir("figBase",'s')
parametry = -1:1;
for(ppp = parametry)   
    Parseval = ppp;
    main
    if(ppp~=end(parametry))
        rmdir(sourcePath,'s')
        clear all;
        close all;    
    end
end

return

fromPath = "figBase/";
% załaduj z pliku zamiast liczyć
figFilenames = dir(strcat(fromPath,"*.fig"));
saveAllOpenedFigures = 1;
if(~isempty(figFilenames)) % jeśli w katalogu są pliki *.fig
    close all force
    for (i = 1:length(figFilenames))
        [filepath,name,ext] = fileparts(figFilenames(i).name);
        pos = find(name == '_', 1, 'last');
        % nr = str2num(name(pos+1:end));
        % hDest = figure(nr);
        a = openfig(strcat(figFilenames(i).folder,"/",figFilenames(i).name));
        % gcf.Number = nr;
        % copyobj(allchild(a),hDest); % occcures errors
        % close(a)
    end
    backup = 0;
    saveAllOpenedFigures = 1;
    exportPath = "fig4ArticleStyleLudwin/";
else
    main % policz    
    backup = 1;
    saveAllOpenedFigures = 1;
    exportPath = "figBase/"; % UWAGA: domyślny katalog zapisu *.fig
end

% find figure 1000 (last nr)
figHandles = findobj('Type', 'figure');
num = [];
for(i = 1:length(figHandles))
    num = [num; figHandles(i).Number];
end
num = sort(num);

fNrScale = [503 527];

if(saveAllOpenedFigures)
    toSave = num;
    cetroidsFigNr = 2; % for scaling
else
    toSave = [];    
    index = [];
    index = find(num==501);
    toSave = [toSave; num(index)]; 
    index = find(num==502);
    toSave = [toSave; num(index)]; 
    
    if(length(toSave) == 2) % make one figure
        
        figure(toSave(1))
        ax1=gca;
        figure(toSave(2))
        ax2=gca;
        
        nrFigure = 555;
        figure(nrFigure);
        tcl=tiledlayout("flow");
        ax1.Parent=tcl;
        ax1.Layout.Tile=1;
        ax2.Parent=tcl;
        ax2.Layout.Tile=2;
        toSave(nrFigure);
    end
    % 1000 ------Spectrum--------
    index1 = [];
    for(i = 1:length(num))
        if(floor(num(i)/1000) == floor(1))
            index1 = [index1; i];
        end
    end
    cetroidsFigNr=0;
    if(length(index1)>1)
        cetroidsFigNr = num(index1(end)); % it is sorted, so we get the higest numbers
        toSave = [toSave; cetroidsFigNr];
    end
    
    % 2000 ------Centroidy-------
    index2 = [];
    for(i = 1:length(num))
        if(floor(num(i)/1000) == floor(2))
            index2 = [index2; i];
        end
    end
    
    if(length(index2) < 3)
        toThere = length(index2);
    else    
        toThere = 2;
    end
    for i = 2:toThere
     centroFigNr = num(index2(i)); % it is sorted, so we get the higest numbers
     toSave = [toSave; centroFigNr];
    end
    
    if(length(index2) > 5)
     cetroidsFigNr = num(index2(length(index2)-2)); 
     toSave = [toSave; cetroidsFigNr];
    end
    
    % 3000 ------
    
    % 4000 ------
    for(i = 4001:4003)
        i = find(num==i);
        toSave = [toSave; num(i)]; 
    end
end

for (i = toSave')
    % if(i<5) continue;  end % for crashing matlab problem   
    figure(i);
    % if eg 1014 
    
    mnoznik = 1;
    if( (i == cetroidsFigNr | find(i==fNrScale)) ) mnoznik = 1.6; end
    if(i==527) mnoznik = 2.5; end;
    figPW("path", exportPath, "exportPDF", 1,"openFolder",1,"saveCopyFig", ...
        backup,"skipSaveAs",1, "scale", mnoznik,"styleLudwin",1);
end