function [constStruct] = loadConstants(loadVar, pathSeparator)
    %loadConstants - Loads the package specific constants to fascilitate importing
    %
    % Syntax: [constStruct] = loadConstants(package)
        CONSTANTS_FILE = pwd + pathSeparator + "constants.mat";
        constStruct = load(CONSTANTS_FILE, loadVar);
        constStruct = constStruct.(loadVar);
end