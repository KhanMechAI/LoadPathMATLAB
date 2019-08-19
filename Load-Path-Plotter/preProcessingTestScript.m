function [] = preProcessingTestScript()   
    %TODO: mkdir a preprocessing folder for saving data
    A = cell(0,1);
    B = cell(0,4);
    C = cell(0,19);
    opt = {'MultipleDelimsAsOne',true};
    F_NAME = '/MATLAB Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/Simulation Files/ds.dat'
    OUT_PATH = '/MATLAB Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/_output_data/'
    NODE_FORMAT = '%d32%f%f%f'
    ELEMENT_FORMAT = repmat('%d32',[1,19])
    REGEXP_NODES_ELEMS = '/com,\*+\s(?<dataType>Nodes|Elements)'
    MODEL_DATA = '\<(?!the|for\>)(?<modelData>[\w-\d])+'
    STR_NODES = '%[/com,*********** Nodes]'
    % START_NODES = 
    fileId = fopen(F_NAME,'rt')
    str = fgetl(fileId);
    while ischar(str)
        [match, nonMatch] = regexp(str, REGEXP_NODES_ELEMS, 'names', 'split');
        if ~isempty(match)
            [submatch, nonMatch] = regexp(nonMatch{2}, MODEL_DATA, 'names')
            [matObj, variable] = initialiseMatFile(fileId, OUT_PATH);
        end
        str = fgetl(fileId);
    end
    fclose(fileId);
        function [matFileObject, storedVariable] = initialiseMatFile(fileId, fullPath)
            %getArrayLength - This function returns the a scalar used to preallocate memory for arrays used when preprocessing data
            %
            % Syntax: [dataArray] = myFun(input)
            %
            % Long description
            %TODO: ADD IDENTIFIER TO OUTPUT FILE NAMES TO DISTINGUISH BETWEEN NODES AND ELEMENTS
            caseType = match.dataType;
            
            switch caseType
                case 'Nodes'
                    outputFormat = NODE_FORMAT;
                case 'Elements'
                    outputFormat = ELEMENT_FORMAT;
                    fgetl(fileId)
                otherwise
                    return
            end
        
            str = fgetl(fileId);
            str = split(str,',');
            arrayLength = str2num(str{end})-1;
            outputPath = [fullPath strjoin({submatch.modelData},'_') '.mat'];
            storedVariable = 'coords';
            switch caseType
                case 'Nodes'
                    coords = zeros(arrayLength,count(outputFormat,'%'));
                case 'Elements'
                    storedVariable = 'connectivity';
                    connectivity = zeros(arrayLength,count(outputFormat,'%'));
                otherwise
                    return
            end
            save(outputPath, storedVariable, '-v7.3');
            matFileObject = matfile(outputPath, 'Writable', true);
            fgetl(fileId)
            matFileObject.(storedVariable) = textscan(fileId,outputFormat,opt{:});
        end
    end


