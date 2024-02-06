exportPath = "figBaseNEW/"; 

figNrList = [
    251
    253
    % 501
    % 502 % done by shringking
    503
    509
    1013
    1023
    2005
    2069
    2078
    4001
    4002
    4003
];

figHandles = findobj('Type', 'figure');
num = [];
for(i = 1:length(figHandles))
    num = [num; figHandles(i).Number];
end
num = sort(num); % figNr

toSave = [];
for(i = figNrList')
    i = find(num==i);
    toSave = [toSave; num(i)]; 
end

for (i = toSave')   
    figure(i);   
    figPW("path", exportPath, "exportPDF", 0,"openFolder",1,"saveCopyFig", 1,"skipSaveAs", 1, "TNR", 0);
end