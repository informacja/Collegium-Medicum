function list = dir2(varargin)
    p = inputParser;
    p.FunctionName = mfilename('name');
    addOptional(p, 'name', '.', @ischar);
    addParameter(p, 'OnlyDirectories', false, @islogical);
    if numel(varargin) == 2 
        varargin = [{'.'}, varargin]; 
    end
    parse(p, varargin{:});
    list = dir(p.Results.name);
    if p.Results.OnlyDirectories
        dirFlags = [list.isdir];
        list = list(dirFlags); 
    end
    list = list(arrayfun(@(x) ~strcmp(x.name(1),'.'), list)); 
end