A = cell(0,1);
B = cell(0,4);
C = cell(0,19);
opt = {'MultipleDelimsAsOne',true};
F_NAME = '/MATLAB-Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/Simulation Files/ds.dat'
NODE_FORMAT = '%d32%f%f%f'
ELEMENT_FORMAT = repmat('%d32',[1,19])
REGEXP_NODES_ELEMS = '/com,\*+\s(?<dataType>Nodes|Elements)'
MODEL_DATA = '\<(?!the|for\>)(?<modelData>[\w-\d])+'
STR_NODES = '%[/com,*********** Nodes]'
% START_NODES = 
fid = fopen(F_NAME,'rt');
str = fgetl(fid);
while ischar(str)
    [match, nonMatch] = regexp(str, REGEXP_NODES_ELEMS, 'names', 'split')
    if ~isempty(match)
        [submatch, nonMatch] = regexp(nonMatch{2}, MODEL_DATA, 'names')
        

        fgetl(fid)
        dataArray = textscan(fid,outputFormat,opt{:});
    end
    str = fgetl(fid);
end
fclose(fid);
function [arrayLength] = getArrayLength(fileId, fullPath)
    %getArrayLength - This function returns the a scalar used to preallocate memory for arrays used when preprocessing data
    %
    % Syntax: [dataArray] = myFun(input)
    %
    % Long description

    fgetl(fid)
    caseType = match.dataType
    switch caseType
        case 'Nodes'
            %leaving as a case switch for future
        case 'Elements'
            fgetl(fid)
        otherwise
            return
    end

    str = fgetl(fid)
    str = split(str,',')
    arrayLength = str2num(str(end))-1
    outputPath = [fullPath strjoin({submatch.modelData},'_')]
    switch caseType
        case 'Nodes'
            coords = zeros(arrayLength,count(NODE_FORMAT,'%'))
            save(outputPath, 'coords', '-v7.3')
        case 'Elements'
            connectivity = zeros(arrayLength,count(ELEMENT_FORMAT,'%'))
            save(outputPath, 'connectivity', '-v7.3')
        otherwise
            return
    end
end



