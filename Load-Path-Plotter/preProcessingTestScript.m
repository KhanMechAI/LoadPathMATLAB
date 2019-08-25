function [] = preProcessingTestScript()   
    %TODO: mkdir a preprocessing folder for saving data
    %TODO: Need to implement a memory size check for large simulations. Then a differnet approach with perhaps tall arrays can be used to manage big data simulations
    % START_NODES = 
    ANSYS = 1;
    STRAND7 = 2;
    softwareImport = 1;
    switch softwareImport
        case ANSYS
            readAnsys()
        case STRAND7
            readStrand7()
        otherwise
            return 
    end
end

function generateDirectories(working, directories)
    fn = fieldnames(directories);
    for k = 1:numel(fn)
        [a, b] = mkdir(working + directories.(fn{k}));
    end
end

function [output] = readAnsys(fileId)
    %readAnsys - ANSYS specific import
    %
    % Syntax: [output] = readAnsys(fileId)
    %
    % Long description
    %TODO: Write regex to match the end of element creation and terminate while loop
    atHome = 0;
    MATLAB_ONLINE = "/MATLAB Drive/";
    MATLAB_HOME = "/Users/khan/MATLAB-Drive/";
    if atHome
        filePath = MATLAB_HOME;
    else    
        filePath = MATLAB_ONLINE;
    end
    pathSep = "/";
    ANSYS = "ANSYS";
    [const] = loadConstants(ANSYS);
    TMP1 = filePath + "Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/Simulation Files/";
    TMP2 = filePath + "Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/";

    const = const.(ANSYS);
    fileName = TMP1 + const.files.ds;

    dirs = const.dirs;

    prepPath = TMP2 + const.dirs.prepPath;
    outPath = TMP2 + const.dirs.outPath;

    generateDirectories(TMP2, [dirs]);

    fmt = const.format;
    regex = const.regex;
   
    fileId = fopen(fileName,'rt');
    str = fgetl(fileId);
    prevStr = "";
    nodeBlockCount = 1;
    elementBlockCount = 1;
    while ischar(str)
        [match, nonMatch] = regexp(str, regex.block, 'names', 'split');
        if ~isempty(match)
            [blockData, nonMatch] = regexp(nonMatch{2}, regex.fields, 'names', 'split');
            switch match.fieldType
                case "n"
                    nodePrepPath = TMP2 + dirs.prepPathN + num2str(nodeBlockCount) + '.mat'
                    readNode(fileId, blockData, fmt.nodes, nodePrepPath);
                    nodeBlockCount = nodeBlockCount + 1;
                case "e"
                    elementPrepPath = TMP2 + dirs.prepPathE + num2str(elementBlockCount) + '.mat'
                    readElement(fileId, elementPrepPath, blockData, prevStr);
                    elementBlockCount = elementBlockCount + 1;
                case "end"
                    continue
            end
        end
        prevStr = str;
        str = fgetl(fileId);
    end
    fclose(fileId);
end

function readElement(fileId, prepPath, blockData, elementData)
    %readElement - Reads element connectivity data and transforms the data ready for saving
    %
    % Syntax: [dataLabels, inputData] = readElement()
    %TODO: Write function to handle elements with more than 8 nodes.

    [nLines] = getNLines(fileId, 2);

    elementFields = uint32(str2double(split(strtrim(nLines(2,:)))));

    nElements = uint32(str2double(blockData(end).fields));

    if ~isempty([blockData.solid])
        nodesPerElement = elementFields(9);
    else
        nodesPerElement = uint32(str2double(blockData(1).fields));
    end

    outputFormat = repmat('%u32',[1,19]);
    connectivity = zeros([nElements,nodesPerElement], 'uint32');
    elementIdx = zeros([nElements,1],'uint32');
    
    tmp = textscan(fileId, outputFormat, 'MultipleDelimsAsOne',true, 'CollectOutput', true);%Read in all data
    % endIdx = length(tmp{1}(1:end));%Get end index

    connectivity(:, :) = [elementFields(12:12 + nodesPerElement - 1)'; tmp{:}(1:end-1, 12:12 + nodesPerElement - 1)];%Read to preallocated array
    elementIdx(:, :) = [elementFields(11); tmp{:}(1:end-1, 11)];

    dataLabels = ["elementIdx", "connectivity"];
    inputData.connectivity = connectivity;
    inputData.elementIdx = elementIdx;

    matFileObject = readToMat(prepPath, dataLabels, inputData);
end

function readNode(fileId, blockData, outputFormat, prepPath)
%readNode - Reads node coordinate data and transforms the data ready for saving
%
% Syntax: [dataLabels, inputData] = readNode()
    maxNodes = uint32(str2double(blockData(end).fields));
    
    fgetl(fileId); %increment the cursor to node block
    
    coords = zeros([maxNodes, 3], 'single');
    nodeIdx = zeros([maxNodes, 1], 'uint32');

    tmp = textscan(fileId,outputFormat,'MultipleDelimsAsOne', true, 'CollectOutput', true);
    endIdx = length(tmp{1}(1, 1:end-1));

    nodeIdx(1:endIdx) = tmp{1}(1:endIdx);
    coords(1:endIdx, :) = tmp{2:end}(1:endIdx);

    dataLabels = ["nodeIdx", "coords"];
    inputData.nodeIdx = nodeIdx;
    inputData.coords = coords;

    matFileObject = readToMat(prepPath, dataLabels, inputData);
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
    arrayLength = str2double(str{end});
end

function [matFileObject] = readToMat(prepPath, dataLabels, inputData)
%readToMat - Description
%
% Syntax: [outputFmt] = readToMat(input)
%
% Long description
%TODO: save as structure with variable names - lookup save function for help
    save(prepPath, '-v7.3', '-struct', 'inputData')
    matFileObject = matfile(prepPath, 'Writable', true);
    for k = 1:length(inputData)
        matFileObject.(char(dataLabels(k,1))) = inputData.(char(dataLabels(k,1)));
    end
end

function [nameArray] = nameGen(inputData)
    %nameGen - This function was to concatenate a string to an array for
    %saving. Unused currently.
    nameArray(:,1) = "inputData." + inputData(:,:);
end

function [varargout] = getElementType(elementTypeString)
%getElementType - Parses the element type and returns the connectivity
%
% Syntax: [connectivity] = getElementType(elementTypeString)
%
% Long description
    elementType = split(elementTypeString,',')
    iBody = uint32(str2double(elementType(2)))
    element = uint32(str2double(elementType(3)))

    temp = {element, iBody}
    for k = 1:nargout
        varargout{k} = temp{k}
    end
end

function [varargout] = parseEblock(eblockString)
%parseEblock - This function processes the element type string and returns the format of the element definition in the ds.dat file and the element connectivity
%
% Syntax: [varargout] = parseEblock(eblockString)
%
% Long description
    eblock = split(eblockString,',')
    nElements = uint32(str2double(eblock(end)))
    isSolid = false
    if any(strcmp(eblock, 'solid'))
        isSolid = true
    end
    temp = {nElements, isSolid}
    for k = 1:nargout
        varargout{k} = temp{k}
    end
end

function [output] = getReadFormat(elementType, isSolid, nNodes, varargin)
%getReadFormat - Description
%
% Syntax: [output] = getReadFormat(elementType, )
%
% Long description
%TODO: Need to setup the read formats and also the other 1-12 fields and separate from the nodal connectivity. Need to store each field in its own variable in the mat file. thsi function should output the formats to then write into the mat file.
    if isSolid
        nNodesPerElement = varargin(9)
        nodeReadFormat = repmat('%u32',[1,nNodesPerElement]);
        solidCell
    end
    
end

function [nLines] = getNLines(fileId, n)
    %getNLines - Iterates over 'n' lines of a file returning each in a cell of a cell array
    %
    % Syntax: [] = getNLines(input)
    %
    % Long description
    nLines = strings(n,1);
    for k=1:n
        nLines(k,:) = fgetl(fileId);
    end
end

function [connectivity] = getConnectivity(elementType)
%getConnectivity - Returns the connectivity of a specific element type.
%
% Syntax: [connectivity] = getConnectivity(elementType)
%
% Long description
    
end

function [constStruct] = loadConstants(softwareImport)
    %loadConstants - Loads the package specific constants to fascilitate importing
    %
    % Syntax: [constStruct] = loadConstants(softwareImport)
        CONSTANTS_FILE = pwd + "/" + "constants.mat";
        constStruct = load(CONSTANTS_FILE, softwareImport);
    end
    
        
function [matFileObject] = initialiseMatFile(fileId, fullPath, blockType, outputFormat)
    %getArrayLength - This function returns the a scalar used to preallocate memory for arrays used when preprocessing data
    %
    % Syntax: [dataArray] = myFun(input)
    %
    % Long description
    %TODO: ADD IDENTIFIER TO OUTPUT FILE NAMES TO DISTINGUISH BETWEEN NODES AND ELEMENTS
    prepPath = [fullPath strjoin({submatch.modelData},'_') '.mat'];

    caseType = match.dataType;
    caseType = nodeOrElement(caseType);

    if caseType
        fgetl(fileId);
        outputFormat = formatElements;
    elseif caseType > -1
        outputFormat = formatNodes;
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

    matFileObject = readToMat(prepPath, dataLabels, inputData);

end