function [] = preProcessingTestScript()   
    %TODO: mkdir a preprocessing folder for saving data
    %TODO: Need to implement a memory size check for large simulations. Then a differnet approach with perhaps tall arrays can be used to manage big data simulations
    A = cell(0,1);
    B = cell(0,4);
    C = cell(0,19);
    opt = {'MultipleDelimsAsOne',true};
    F_NAME = '/MATLAB Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/Simulation Files/ds.dat';
    OUT_PATH = '/MATLAB Drive/Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/_output_data/';
    NODE_FORMAT = '%d32%f32%f32%f32';
    ELEMENT_FORMAT = repmat('%d32',[1,19]);
    REGEXP_NODES_ELEMS = '/com,\*+\s(?<dataType>Nodes|Elements)';
    MODEL_DATA = '\<(?!the|for\>)(?<modelData>[\w-\d])+';
    STR_NODES = '%[/com,*********** Nodes]';
    % START_NODES = 
    fileId = fopen(F_NAME,'rt');
    str = fgetl(fileId);
    while ischar(str)
        [match, nonMatch] = regexp(str, REGEXP_NODES_ELEMS, 'names', 'split');
        if ~isempty(match)
            [submatch, nonMatch] = regexp(nonMatch{2}, MODEL_DATA, 'names');
            [matObj] = initialiseMatFile(fileId, OUT_PATH);
        end
        str = fgetl(fileId);
    end
    fclose(fileId);
        function [matFileObject] = initialiseMatFile(fileId, fullPath)
            %getArrayLength - This function returns the a scalar used to preallocate memory for arrays used when preprocessing data
            %
            % Syntax: [dataArray] = myFun(input)
            %
            % Long description
            %TODO: ADD IDENTIFIER TO OUTPUT FILE NAMES TO DISTINGUISH BETWEEN NODES AND ELEMENTS
            outputPath = [fullPath strjoin({submatch.modelData},'_') '.mat'];
            
            caseType = match.dataType;
            caseType = nodeOrElement(caseType);

            if caseType
                fgetl(fileId);
                outputFormat = ELEMENT_FORMAT;
            elseif caseType > -1
                outputFormat = NODE_FORMAT;
            else
                return
            end
            
            [arrayLength] = getArrayLength(fileId);
            fgetl(fileId);
                
            if caseType
                [dataLabels, inputData] = readElement();
            else
                [dataLabels, inputData] = readNode();
            end

            matFileObject = readToMat(outputPath, dataLabels, inputData);

            function [dataLabels, inputData] = readElement()
                %readElement - Reads element connectivity data and transforms the data ready for saving
                %
                % Syntax: [dataLabels, inputData] = readElement()
                connectivity = zeros(arrayLength,count(outputFormat,'%'), 'int32');
            
                tmp = textscan(fileId,outputFormat,opt{:});
                endIdx = length(tmp{1}(1:end-1));
                connectivity(1:endIdx, :) = tmp{1}(1:end-1);
            
                dataLabels = ["connectivity"];
                inputData.connectivity = connectivity;
            end
            
            function [dataLabels, inputData] = readNode()
            %readNode - Reads node coordinate data and transforms the data ready for saving
            %
            % Syntax: [dataLabels, inputData] = readNode()
                coords = zeros(arrayLength, 3, 'single');
                nodeIdx = zeros(arrayLength, 1, 'int32');
                
                tmp = textscan(fileId,outputFormat,opt{:});
                endIdx = length(tmp{1}(1:end-1));

                nodeIdx(1:endIdx) = tmp{1}(1:end-1);
                coords(1:endIdx+1,:) = cell2mat({tmp{2:end}});
                coords = coords(1:end-1,:);

                dataLabels = ["nodeIdx", "coords"];
                inputData.nodeIdx = nodeIdx;
                inputData.coords = coords;
            end
        end
    end

function [retVal] = nodeOrElement(caseType)
    switch caseType
        case 'Nodes'
            retVal = 0;
        case 'Elements'
            retVal = 1;
        otherwise
            retVal = -1;
    end
    return
end


function [arrayLength] = getArrayLength(fileId)
%getArrayLength - Returns the maximum length of the array to be read
%
% Syntax: [arrayLength] = getArrayLength(fileId)
    str = fgetl(fileId); str = split(str,',');
    arrayLength = str2num(str{end});
end

function [matFileObject] = readToMat(outputPath, dataLabels, inputData)
%readToMat - Description
%
% Syntax: [outputFmt] = readToMat(input)
%
% Long description
%TODO: save as structure with variable names - lookup save function for help
    save(outputPath, '-v7.3', '-struct', 'inputData')
    matFileObject = matfile(outputPath, 'Writable', true);
    for k = 1:length(inputData)
        matFileObject.(dataLabels(k,1)) = inputData.(dataLabels(k,1));
    end
end

function [nameArray] = nameGen(inputData)
    %nameGen - This function was to concatenate a string to an array for
    %saving. Unused currently.
    nameArray(:,1) = "inputData." + inputData(:,:);
end


%PENDING DELETE
function [] = getNLines(fileId, n)
%getNLines - Iterates over 'n' lines of a file
%
% Syntax: [] = getNLines(input)
%
% Long description
    for k=1:n
        fgetl(fileId)
    end
end



