function [] = preProcess(package, simulationDir, workingDir)   
    %TODO: Need to implement a memory size check for large simulations. Then a differnet approach with perhaps tall arrays can be used to manage big data simulations
    %TODO: Build out the preprocessing functions around connectivity and
    %element faces etc.
    %TODO: read in the seed file?
    atHome = 1;
    MATLAB_ONLINE = "/MATLAB Drive/";
    MATLAB_HOME = "/Users/khan/MATLAB-Drive/";
    if atHome
        filePath = MATLAB_HOME;
    else    
        filePath = MATLAB_ONLINE;
    end
    pathSep = "/";

    [general] = loadConstants("general", pathSep);
    ANSYS = general.package.ANSYS;
    STRAND7 = general.package.STRAND7;
    package = 1;

    general.dirs.simulationDir = filePath + "Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/Simulation Files/";
    general.dirs.workingDir = filePath + "Load-Path-Plotter/LoadPathMATLAB/Load-Path-Plotter/Examples/Example10 - Notched Plate Coarse/";

    general.pathSep = pathSep;


    switch package
        case ANSYS
            readAnsys(general)
        case STRAND7
            readStrand7()
        otherwise
            return 
    end


end

function generateDirectories(general)
    %TODO:  Future imporvement to clear the directories prior to creation.
    fn = fieldnames(general.dirs(:));
    for k = 1:numel(fn)
        [a, b] = mkdir(general.dirs.workingDir + general.dirs.(fn{k}));
    end
end

function [] = readAnsys(general)
    %readAnsys - ANSYS specific import
    %
    % Syntax: [output] = readAnsys(fileId)
    ANSYS = "ANSYS";
    pathSep = general.pathSep;
    [const] = loadConstants(ANSYS, pathSep);

    generateDirectories(general);

    readDsDat()
    readNodalSol()
    generateData(general)
 

    function varargout = readDsDat()
        %readDsDat - Imports data from ds.dat file.
        %
        % Syntax: varargout = readDsDat()
        fileName = general.dirs.simulationDir + const.files.ds;

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
                        readNode(fileId, blockData, const, general);
                    case "e"
                        [elementType, nonMatch] = regexp(prevStr, regex.elementType, 'names', 'split');
                        elementType = uint32(str2double(elementType.elementType));
                        readElement(fileId, const, general, blockData, elementType, elementBlockCount);
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
        fileName = general.dirs.simulationDir + const.files.nodalSol;
        
        dirs = general.dirs;
        fmt = const.format;
        regex = const.regex;

        stressPrepPath = dirs.workingDir + dirs.prepPathS + general.files.stress;
    
        fileId = fopen(fileName,'rt');
        str = convertCharsToStrings(fgetl(fileId));
        lineTest = true;

        maxNodes = getData(general, "g", general.varNames.maxNodes);

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

        dataLabels = [general.varNames.nodalStressIdx, general.varNames.stress];
        inputData.stress = stress;
        inputData.nodalStressIdx = nodalStressIdx;

        saveToMat(stressPrepPath, dataLabels, inputData, false);
    end

    dataLabels = [general.varNames.general];
    inputData.general = general;
    generalOutPath = general.dirs.workingDir + general.path.general;

    saveToMat(generalOutPath, dataLabels, inputData, false);
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

function readElement(fileId, const, general, blockData, elementType, elementBlockCount)
    %readElement - Reads element connectivity data and transforms the data ready for saving
    %
    % Syntax: [dataLabels, inputData] = readElement()
    %TODO: Write function to handle elements with more than 8 nodes.
    dirs = general.dirs;
    fmt = const.format;

    elementPrepPath = dirs.workingDir + dirs.prepPathE + num2str(elementBlockCount) + general.files.matExt;
    generalConstPath = dirs.workingDir + general.path.general;
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

    dataLabels = [general.varNames.elementIdx, general.varNames.connectivity];
    inputData.connectivity = connectivity;
    inputData.elementIdx = elementIdx;

    saveToMat(elementPrepPath, dataLabels, inputData, false);
    
    dataLabels = [general.varNames.elementData];

    general.elementData = [nElements, nodesPerElement, elementType];
   
    saveToMat(generalConstPath, dataLabels, general, true);
end

function readNode(fileId, blockData, const, general)
    %readNode - Reads node coordinate data and transforms the data ready for saving
    %
    % Syntax: [dataLabels, inputData] = readNode()
    dirs = general.dirs;
    fmt = const.format;
    nodePrepPath = dirs.workingDir + dirs.prepPathN + general.files.nodes;
    generalConstPath = dirs.workingDir + general.path.general;
    maxNodes = uint32(str2double(blockData(end).fields));
    
    fgetl(fileId); %increment the cursor to node block
    
    general.local.maxNodes = maxNodes;

    coords = zeros([maxNodes, 3], 'single');
    nodeIdx = zeros([maxNodes, 1], 'uint32');

    tmp = textscan(fileId, fmt.nodes,'MultipleDelimsAsOne', true, 'CollectOutput', true);
    endIdx = length(tmp{1, 1}(:));

    nodeIdx(1:endIdx) = tmp{1}(1:endIdx);
    coords(1:endIdx, :) = tmp{2:end}(1:endIdx, :);

    dataLabels = [general.varNames.nodeIdx, general.varNames.coords];
    inputData.nodeIdx = nodeIdx;
    inputData.coords = coords;

    saveToMat(nodePrepPath, dataLabels, inputData, false);

    dataLabels = [general.varNames.maxNodes];
    general.maxNodes = maxNodes;
    
    saveToMat(generalConstPath, dataLabels, general, false);

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

function varargout = saveToMat(prepPath, dataLabels, inputData, appendData)
    %saveToMat - Reads data to matfile
    %
    % Syntax: [matFileObject] = saveToMat(prepPath, dataLabels, inputData)

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

function generateData(general)

    elementData = getData(general, "g", general.varNames.elementData);
    % nodeIdx = getData(general, "n", general.varNames.nodeIdx);
    maxNodes = getData(general, "g", general.varNames.maxNodes);
    nBodies = size(elementData, 1);
    maxElements = 8;
    nodalConnectivity = zeros(maxNodes, maxElements);

    for k = 1:nBodies
        faceArray = faceDef(elementData(k,3));
        [nFaces, nNodesPerFace] = size(faceArray);
        bodyFaceArray = zeros([elementData(k,1)*nFaces, nNodesPerFace], 'uint32');
        connectivity = zeros([elementData(k,1), elementData(k,2)], 'uint32');
        connectivity(:)  = getData(general,"e", general.varNames.connectivity, k);
        bodyFaceArray(:) = reshape(connectivity(:,faceArray),elementData(k,1)*nFaces, nNodesPerFace);

        
        [~, ind, idx] = unique(sort(bodyFaceArray, 2), 'rows');
        a = histcounts(idx, length(ind)) < 2;
        tmp = bodyFaceArray(ind,:);
        b = tmp(a,:);
        plotFaces(general, b)
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

function [coords] = getCoords(general, nodeNumber)
    %getCoords - Description
    %
    % Syntax: [coords] = getCoords(nodeNumber)
    %
    % Long description
    nodes = matfile(general.dirs.workingDir + general.path.nodes);
    idx = nodes.nodeIdx(nodeNumber, :);
    coords = nodes.coords(idx,:);
end

function [faceCoords] =  getOuterFaces(general, faceArray)
    % TODO - Save face coords to matfile
    uNodes = reshape(faceArray,[],1);
    [uNodes, ~, ic] = unique(uNodes);
    coords = getCoords(general, uNodes);
    faceCoords = coords(ic, :);
    [nRows, nCols] = size(faceArray);
    faceCoords = reshape(faceCoords, nRows, nCols,[]);
    faceCoords = permute(faceCoords, [2 1 3]);
    % TODO - It is currently element 185 specific - make general. i.e. I'm making triangles from the 4 nodes that define the face. Needs to be adapted for elements with more or less nodes.
    generateLocatingVariables(general, [faceCoords([1:3], :,:), faceCoords([1,3,4], :,:)]);
end


function [varargout] = generateLocatingVariables(general, faceVerticies)
    %generateLocatingVariables - generates and stores the preprocessed data for locating the point in mesh.
    % varargout: true if error in data, false otherwise
    %
    % Syntax: [varargout] = generateLocatingVariable(faceVerticies)
    V0 = squeeze(faceVerticies(1, :, :));
    V1 = squeeze(faceVerticies(2, :, :));
    V2 = squeeze(faceVerticies(3, :, :));
    u = V1 - V0;
    v = V2 - V0;
    n = cross(u, v, 2);
    uu = dot(u,u, 2);
    uv = dot(u,v, 2);
    vv = dot(v,v, 2);
    D = uv .* uv - uu .* vv;

    errorInData = true;
    if all(~n, 'all')
        errorInData  =false;
    end

    dataLabels = [general.varNames.V0, general.varNames.u, general.varNames.v, general.varNames.n, general.varNames.uu, general.varNames.uv, general.varNames.vv, general.varNames.D, general.varNames.infPoint];

    inputData.V0 = V0;
    inputData.u = u;
    inputData.v = v;
    inputData.n = n;
    inputData.uu = uu;
    inputData.uv = uv;
    inputData.vv = vv;
    inputData.D = D;
    inputData.infPoint = general.constants.infPoint;

    outPath = general.dirs.workingDir + general.path.utilities;

    saveToMat(outPath, dataLabels, inputData, false);

    varargout{1} = [];
    for k = 1:nargout
        varargout{k} = errorInData;
    end
end

function plotFaces(general, faceArray)
    %plotFaces - Description
    %
    % Syntax: plotFaces(faceArray)

    faceCoords = getOuterFaces(general, faceArray);

    triangulatedFaceCoords = [faceCoords([1:3], :,:), faceCoords([1,3,4], :,:)];
    hold on

    patch('XData',faceCoords(:, :, 1),'YData',faceCoords(:, :, 2),'ZData',faceCoords(:, :, 3), 'EdgeColor','black','FaceColor','none', 'EdgeAlpha', 0.3)

    testPoint = [1 1 1];

    plot3(testPoint(1), testPoint(2), testPoint(3), '*', 'MarkerFaceColor', 'red')
    hold off

    utilities = getData(general, "u");
    triIntersect(testPoint, utilities)

end

function [in] = triIntersect(point, utilities)
    %triIntersect - Description
    %
    % Syntax: [in] = triIntersect(point, faceVerticies)
    %
    % Adapted from: http://geomalgorithms.com/a06-_intersect-2.html
    in = false;

    ray = utilities.infPoint - point;

    w0 = point - utilities.V0;

    a = -dot(utilities.n, w0, 2);
    b = utilities.n * ray.';

    r = a ./ b;

    I = point + (r * ray);

    w = I - utilities.V0;
    wu = dot(w,utilities.u, 2);
    wv = dot(w,utilities.v, 2);

    s = (utilities.uv .* wv - utilities.vv .* wu) ./ utilities.D ;

    t = (utilities.uv .* wu - utilities.uu .* wv) ./ utilities.D ;

    if sum((t > 0.0 & (s + t) > 1.0) & (r > 0) & (s > 0.0) & (s < 1.0)) > 0
        in = true;
    end
end

