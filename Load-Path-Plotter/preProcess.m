function [] = preProcess(package, simDir, workingDir)   
    %TODO: Need to implement a memory size check for large simulations. Then a differnet approach with perhaps tall arrays can be used to manage big data simulations
    %TODO: Build out the preprocessing functions around connectivity and
    %element faces etc.
    %TODO: read in the seed file?
    [genConst] = loadConstants("GENERAL");
    ANSYS = genConst.package.ANSYS;
    STRAND7 = genConst.package.STRAND7;
    package = 1;

    atHome = 0;
    MATLAB_ONLINE = "/MATLAB Drive/";
    MATLAB_HOME = "/Users/khan/MATLAB-Drive/";
    if atHome
        filePath = MATLAB_HOME;
    else    
        filePath = MATLAB_ONLINE;
    end
    pathSep = "/";
    
    simDir = filePath + "Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Fine/Simulation Files/";
    workingDir = filePath + "Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Fine/";

    switch package
        case ANSYS
            readAnsys(simDir, workingDir)
        case STRAND7
            readStrand7()
        otherwise
            return 
    end
end

function generateDirectories(working, directories)
    %TODO:  Future imporvement to clear the directories prior to creation.
    fn = fieldnames(directories);
    for k = 1:numel(fn)
        [a, b] = mkdir(working + directories.(fn{k}));
    end
end

function [] = readAnsys(simDir, workingDir)
    %readAnsys - ANSYS specific import
    %
    % Syntax: [output] = readAnsys(fileId)
    ANSYS = "ANSYS";
            
    [const] = loadConstants(ANSYS);
    [genConst] = loadConstants("GENERAL");

    generateDirectories(workingDir, [genConst.dirs]);

    genConst.dirs.workingDir = workingDir;

    readDsDat()
    readNodalSol()
    generateData(genConst)
 

    function varargout = readDsDat()
        %readDsDat - Imports data from ds.dat file.
        %
        % Syntax: varargout = readDsDat()
        fileName = simDir + const.files.ds;

        fmt = const.format;
        regex = const.regex;
    
        fileId = fopen(fileName,'rt');
        str = fgetl(fileId);
        prevStr = "";
        elementBlockCount = 1;
        blockElementTypeMapping = {};
    
        while ~atEnd(str, regex.elemEnd)
            [match, nonMatch] = regexp(str, regex.block, 'names', 'split');
            if ~isempty(match)
                [blockData, nonMatch] = regexp(nonMatch{2}, regex.fields, 'names', 'split');
                switch match.fieldType
                    case "n"
                        readNode(fileId, blockData, const, genConst);
                    case "e"
                        [elementType, nonMatch] = regexp(prevStr, regex.elementType, 'names', 'split');
                        elementType = uint32(str2double(elementType.elementType));
                        readElement(fileId, const, genConst, blockData, elementType, elementBlockCount);
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
    
    function varargout = readNodalSol()
        %readNodalSol - reads in the stress data from nodalSolution.txt
        %
        % Syntax: varargout = readNodalSol()
        fileName = simDir + const.files.nodalSol;
        
        dirs = genConst.dirs;
        fmt = const.format;
        regex = const.regex;

        stressPrepPath = dirs.workingDir + dirs.prepPathS + genConst.files.stress;
    
        fileId = fopen(fileName,'rt');
        str = convertCharsToStrings(fgetl(fileId));
        lineTest = true;

        maxNodes = getData(genConst, "g", genConst.varNames.maxNodes);

        stress = zeros([maxNodes, 6], 'single');
        nodalStressIdx = zeros([maxNodes, 1], 'uint32');
        
        while lineTest
            [match, nonMatch] = regexp(str, regex.nodeSolBlock, 'names', 'noemptymatch');
            if ~isempty(match)
                if length(match)>1 && ~isempty(match{1,2}(:))
                    tmp = textscan(fileId, fmt.nodalStress,'MultipleDelimsAsOne', true, 'CollectOutput', true);
                    lineTest = false;
                end
            end
            str = convertCharsToStrings(fgetl(fileId));
        end

        endIdx = length(tmp{1}(1:end, 1));

        nodalStressIdx(1:endIdx) = tmp{1}(1:endIdx);
        stress(1:endIdx, :) = tmp{2:end}(1:endIdx, :);

        dataLabels = [genConst.varNames.nodalStressIdx, genConst.varNames.stress];
        inputData.stress = stress;
        inputData.nodalStressIdx = nodalStressIdx;

        readToMat(stressPrepPath, dataLabels, inputData, false);
    end

end

function [loadedVar] = getData(genConst, dataClass, varargin)
    varToLoad = varargin{1,1};
    if nargin > 3
        body = string(varargin{1,2});
    end
    switch dataClass
        case "e"
            curDir = genConst.dirs.workingDir + genConst.dirs.prepPathE...
                     + body + genConst.files.matExt;
        case "n"
            curDir = genConst.dirs.workingDir + genConst.path.nodes;
        case "s"
            curDir = genConst.dirs.workingDir + genConst.path.stress;
        case "g"
            curDir = genConst.dirs.workingDir + genConst.path.globalData;
    end
    loadedVar = load(curDir, varToLoad);
    loadedVar = loadedVar.(varToLoad);
end

function endElements = atEnd(str, pat)
    %atEnd - Returns true if at the end of the elements. False otherwise.
    %
    % Syntax: endElements = atEnd(str)
    endElements = false;

    [match, nonMatch] = regexp(str, pat, 'names', 'split');

    if ~isempty(match)
        if strcmp(match.sectionKey, "elem")
            if strcmp(match.context, "end")
                endElements = true;
            end
        end
    end
end

function readElement(fileId, const, genConst, blockData, elementType, elementBlockCount)
    %readElement - Reads element connectivity data and transforms the data ready for saving
    %
    % Syntax: [dataLabels, inputData] = readElement()
    %TODO: Write function to handle elements with more than 8 nodes.
    dirs = genConst.dirs;
    fmt = const.format;

    elementPrepPath = dirs.workingDir + dirs.prepPathE + num2str(elementBlockCount) + genConst.files.matExt;
    generalConstPath = dirs.workingDir + genConst.path.globalData;
    [nLines] = getNLines(fileId, 2);

    elementFields = uint32(str2double(split(strtrim(nLines(2,:)))));

    nElements = uint32(str2double(blockData(end).fields));

    if ~isempty([blockData.solid])
        nodesPerElement = elementFields(9);
    else
        nodesPerElement = uint32(str2double(blockData(1).fields));
    end

    connectivity = zeros([nElements,nodesPerElement], 'uint32');
    elementIdx = zeros([nElements,1],'uint32');
    
    tmp = textscan(fileId, fmt.elements, 'MultipleDelimsAsOne',true, 'CollectOutput', true);%Read in all data

    connectivity(:, :) = [elementFields(12:12 + nodesPerElement - 1)'; tmp{:}(1:end-1, 12:12 + nodesPerElement - 1)];%Read to preallocated array
    elementIdx(:, :) = [elementFields(11); tmp{:}(1:end-1, 11)];

    dataLabels = [genConst.varNames.elementIdx, genConst.varNames.connectivity];
    inputData.connectivity = connectivity;
    inputData.elementIdx = elementIdx;

    readToMat(elementPrepPath, dataLabels, inputData, false);
    
    dataLabels = [genConst.varNames.elementData];

    globalData.elementData = [nElements, nodesPerElement, elementType];
   
    readToMat(generalConstPath, dataLabels, globalData, true);
end

function readNode(fileId, blockData, const, genConst)
    %readNode - Reads node coordinate data and transforms the data ready for saving
    %
    % Syntax: [dataLabels, inputData] = readNode()
    dirs = genConst.dirs;
    fmt = const.format;
    nodePrepPath = dirs.workingDir + dirs.prepPathN + genConst.files.nodes;
    generalConstPath = dirs.workingDir + genConst.path.globalData;
    maxNodes = uint32(str2double(blockData(end).fields));
    
    fgetl(fileId); %increment the cursor to node block
    
    genConst.local.maxNodes = maxNodes;

    coords = zeros([maxNodes, 3], 'single');
    nodeIdx = zeros([maxNodes, 1], 'uint32');

    tmp = textscan(fileId, fmt.nodes,'MultipleDelimsAsOne', true, 'CollectOutput', true);
    endIdx = length(tmp{1, 1}(:));

    nodeIdx(1:endIdx) = tmp{1}(1:endIdx);
    coords(1:endIdx, :) = tmp{2:end}(1:endIdx, :);

    dataLabels = [genConst.varNames.nodeIdx, genConst.varNames.coords];
    inputData.nodeIdx = nodeIdx;
    inputData.coords = coords;

    readToMat(nodePrepPath, dataLabels, inputData, false);

    dataLabels = [genConst.varNames.maxNodes];
    globalData.maxNodes = maxNodes;
    
    readToMat(generalConstPath, dataLabels, globalData, false);

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
end


function [arrayLength] = getArrayLength(fileId)
%getArrayLength - Returns the maximum length of the array to be read
%
% Syntax: [arrayLength] = getArrayLength(fileId)
    str = fgetl(fileId); str = split(str,',');
    arrayLength = str2double(str{end});
end

function varargout = readToMat(prepPath, dataLabels, inputData, appendData)
%readToMat - Reads data to matfile
%
% Syntax: [matFileObject] = readToMat(prepPath, dataLabels, inputData)
%TODO: Need to rework to be able to handle any input variable or struct. Currently only a var named 'inputData' will work.
    if appendData
        save(prepPath, '-v7.3', '-struct', 'inputData', '-append')
    else
        save(prepPath, '-v7.3', '-struct', 'inputData')
    end
    matFileObject = matfile(prepPath, 'Writable', true);
    for k = 1:length(inputData)
        matFileObject.(char(dataLabels(k,1))) = inputData.(char(dataLabels(k,1)));
    end

    nout = max(nargout, 1) - 1;

    varargout{1} = [];
    for k = 1:nout
        varargout{k} = matFileObject;
    end
end

function [varargout] = getElementType(elementTypeString)
    %getElementType - Parses the element type and returns the connectivity
    %
    % Syntax: [connectivity] = getElementType(elementTypeString)

    elementType = split(elementTypeString,',')
    iBody = uint32(str2double(elementType(2)))
    element = uint32(str2double(elementType(3)))

    temp = {element, iBody}
    for k = 1:nargout
        varargout{k} = temp{k}
    end
end

function [nLines] = getNLines(fileId, n)
    %getNLines - Iterates over 'n' lines of a file returning each in a cell of a cell array
    %
    % Syntax: [nLines] = getNLines(input)
    nLines = strings(n,1);
    for k=1:n
        nLines(k,:) = fgetl(fileId);
    end
end

function [connectivity] = getConnectivity(elementType)
    %getConnectivity - Returns the connectivity of a specific element type.
    %
    % Syntax: [connectivity] = getConnectivity(elementType)

    
end

function [constStruct] = loadConstants(loadVar)
    %loadConstants - Loads the package specific constants to fascilitate importing
    %
    % Syntax: [constStruct] = loadConstants(package)
        CONSTANTS_FILE = pwd + "/" + "constants.mat";
        constStruct = load(CONSTANTS_FILE, loadVar);
        constStruct = constStruct.(loadVar);
end

function [faces] = faceDef(elementType)
    I = 1; J = 2; K = 3; L = 4; M = 5; N = 6; O = 7; P = 8;
    Q = 9; R = 10; S = 11; T = 12; U = 13; V = 14; W = 15; X = 16;
    Y = 17; Z = 18; A = 19; B = 20;
    switch elementType
        
        case 185
            faces = [...
                [J,I,L,K];...
                [I,J,N,M];...
                [J,K,O,N];...
                [K,L,P,O];...
                [L,I,M,P];...
                [M,N,O,P];...
                ];

        case 186
            
            faces = [...
                [J,Q,I,T,L,S,K,R,];...
                [I,Q,J,Z,N,U,M,Y,];...
                [J,R,K,A,O,V,N,Z,];...
                [K,S,L,B,P,W,O,A,];...
                [L,T,I,Y,M,X,P,B,];...
                [M,U,N,V,O,W,P,X,];...
                ];

        case 187
            faces = [...
                [J;I;K;],...
                [I;J;L;],...
                [J;K;L;],...
                [K;I;L;],...
                ];
        % case 181
        %     % TODO: Need a routine to construct the shell elements with some thickness.
        %     faces = [...
        %         [I;J;K;L;],...
        %         [I;J;K;L;],...
        %         [J;I;],...
        %         [K;J;],...
        %         [L;K;],...
        %         [I;L;],...
        %         ];
    end
end

function generateData(genConst)

    elementData = getData(genConst, "g", genConst.varNames.elementData);
    % nodeIdx = getData(genConst, "n", genConst.varNames.nodeIdx);
    maxNodes = getData(genConst, "g", genConst.varNames.maxNodes);
    nBodies = size(elementData, 1);
    maxElements = 8;
    nodalConnectivity = zeros(maxNodes, maxElements);

    for k = 1:nBodies
        faceArray = faceDef(elementData(k,3));
        [nFaces, nNodesPerFace] = size(faceArray);
        bodyFaceArray = zeros([elementData(k,1)*nFaces, nNodesPerFace], 'uint32');
        connectivity = zeros([elementData(k,1), elementData(k,2)], 'uint32');
        connectivity(:)  = getData(genConst,"e", genConst.varNames.connectivity, k);
        bodyFaceArray(:) = reshape(connectivity(:,faceArray),elementData(k,1)*nFaces, nNodesPerFace);

        
        [~, ind, idx] = unique(sort(bodyFaceArray, 2), 'rows');
        a = histcounts(idx, length(ind)) < 2;
        tmp = bodyFaceArray(ind,:);
        b = tmp(a,:);
        plotFaces(genConst, b)
        % duplicate indices
        duplicate_ind = setdiff(1:size(bodyFaceArray, 1), ind);
        % duplicate values
        duplicate_value = bodyFaceArray(duplicate_ind, 3);
        
        [~, uidx, b] = unique(sort(bodyFaceArray, 2), 'rows');
        bodySurfaceFaces = bodyFaceArray(uidx,:);
        for n = 1:maxNodes
            indices = find(any(connectivity==n,2));
            nodalConnectivity(n,1:length(indices)) = indices;
        end
        nodalConnectivity( ~any(nodalConnectivity,2), : ) = [];  %rows
        nodalConnectivity( :, ~any(nodalConnectivity,1) ) = [];  %columns
    end
end

function [coords] = getCoords(genConst, nodeNumber)
    %getCoords - Description
    %
    % Syntax: [coords] = getCoords(nodeNumber)
    %
    % Long description
    nodes = matfile(genConst.dirs.workingDir + genConst.path.nodes);
    idx = nodes.nodeIdx(nodeNumber, :);
    coords = nodes.coords(idx,:);
end

function plotFaces(genConst, faceArray)
    %plotFaces - Description
    %
    % Syntax: plotFaces(faceArray)
    %
    % Long description
    nodes = reshape(faceArray,[],1);
    [uNodes, ~, ic] = unique(nodes);
    coords = getCoords(genConst, uNodes);
    faceCoords = coords(ic, :);
    [nRows, nCols] = size(faceArray);
    faceCoords = reshape(faceCoords, nRows, nCols,[]);
    tmp = permute(faceCoords, [2 1 3]);

    triangulatedFaceCoords = [tmp([1:3], :,:), tmp([1,3,4], :,:)];
    
    triIntersect([0, 0, 0], triangulatedFaceCoords)
    
    patch('XData',X,'YData',Y,'ZData',Z, 'EdgeColor','black','FaceColor','none', 'EdgeAlpha', 0.3)
    
end

function [in] = triIntersect(point, faceVerticies)
    %triIntersect - Description
    %
    % Syntax: [in] = triIntersect(point, faceVerticies)
    %
    % Long description
    tol = 0.000001;
    distantPoint = [1e8 1e8 1e8];

    V0 = squeeze(faceVerticies(1, :, :));
    V1 = squeeze(faceVerticies(2, :, :));
    V2 = squeeze(faceVerticies(3, :, :));
    u = V1 - V0;
    v = V2 - V0;
    n = cross(u,v);

    if all(~n)
        in = false;
        return
    end
    ray = distantPoint - point;
    w0 = point - V0;
    a = -dot(n, w0);
    b = dot(n, ray);
    
    if abs(b) < tol
        if a == 0
            in = true;
            return 
        else
            in = false;
            return
        end
    end

    r = a/b;
    if r < 0
        in = false;
        return
    end

    I = point + r * ray;

    uu = dot(u,u);
    uv = dot(u,v);
    vv = dot(v,v);
    w = I - V0;
    wu = dot(w,u);
    wv = dot(w,v);
    D = uv * uv - uu * vv;

    s = (uv * wv - vv * wu) / D;

    if (s < 0.0) || (s > 1.0)
        in = false;
        return
    end
    
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)
        in = false;
        return 
    end
    in = true;
end