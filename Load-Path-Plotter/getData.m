function [loadedVar] = getData(general, dataClass, varargin)
    if length(varargin) > 0
        varToLoad = varargin{1,1};
        if nargin > 3
            body = string(varargin{1,2});
        end
    end
    switch dataClass
        case "e"
            curDir = general.dirs.workingDir + general.dirs.prepPathE...
                     + body + general.files.matExt;
        case "n"
            curDir = general.dirs.workingDir + general.path.nodes;
        case "s"
            curDir = general.dirs.workingDir + general.path.stress;
        case "g"
            curDir = general.dirs.workingDir + general.path.general;
        case "u"
            curDir = general.dirs.workingDir + general.path.utilities;
    end
    if nargin < 3
        loadedVar = load(curDir);
    else
        loadedVar = load(curDir, varToLoad);
        loadedVar = loadedVar.(varToLoad);
    end
end