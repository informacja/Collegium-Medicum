function [v,k,fromSzU] = LoadNeurosoftData(dirname,k)
files = dir(fullfile(dirname,'**','*.wav'));
datafiles = fullfile({files.folder},{files.name});
cat = 0;
   % fromSzU(1).i = "";
for(d = datafiles)
    fullfilename = char(string(d));
    lastslash_pos = find(fullfilename == '/', 1, 'last');
    filename = fullfilename(lastslash_pos+1:end);

    fullpath =  fullfilename(1:lastslash_pos-1); 
    lastslash_pos = find(fullpath == '/', 1, 'last');
    patientInfo = fullpath(lastslash_pos+1:end); % .

    superfullpath =  fullfilename(1:lastslash_pos-1);
    lastslash_pos = find(superfullpath == '/', 1, 'last');
    categoryInfo = superfullpath(lastslash_pos+1:end); % ..
    if(filename(3)~='1') continue; end % check "1 '1' sEMG.wav"
    
    k = k + 1;
    v(k).infoRDisp = filename; %folder nadrzędny; 
    % v(k).data_nEMG
    [tmpData, freq] = audioread(fullfilename);
    fpom = 2000;
    tmpData = downsample(tmpData,int32(freq/fpom));

    if(categoryInfo(1)=="m")
        cat = cat+1;
        v(k).dataB = tmpData; % Biceps
        v(k).infoTraining = 3;
        v(k).infoBName = "Biceps brachii";
        v(k).segLen = length(tmpData);
    elseif(categoryInfo(1)=="n")
        cat = cat+1;
        v(k).dataB = tmpData; % tibialis anterior
        v(k).infoTraining = 2;
        v(k).infoBName = "Tibialis anterior";
        v(k).segLen = length(tmpData);
    else
        cat = cat+1;
        v(k).data_sEMG = tmpData;
    end
    v(k).infoRecord = patientInfo;
    % % v(k).t = [];

end
for(i=1:k) % todo Unique
    % fromSzU(i).i = v(i).infoRecord;
    tab(i) = string(v(i).infoRecord);
end
aa=  unique(tab);
for(i = 1:numel(aa))
fromSzU(i).i = aa(i);
end
end