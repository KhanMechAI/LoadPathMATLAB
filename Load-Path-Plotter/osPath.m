function pathSeparator = osPath()
    %osPath - Checks the operating system to then ensuren the correct path separators are used for paths.
    %
    % Syntax: pathSeparator = osPath()
    pathSeparator = "/";
    if ispc
        pathSeparator = "\";
        system(strjoin(['taskkill /fi "WINDOWTITLE eq ', modelName,'.pdf"'],''));
    end
end