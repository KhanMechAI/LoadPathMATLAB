function [general] = preProcess(general)   
    %TODO: Need to implement a memory size check for large simulations. Then a differnet approach with perhaps tall arrays can be used to manage big data simulations
    %TODO: Build out the preprocessing functions around connectivity and
    %element faces etc.
    %TODO: read in the seed file?

    pathSep = general.local.pathSep;

    ANSYS = general.package.ANSYS;
    STRAND7 = general.package.STRAND7;
    package = 1;


    switch package
        case ANSYS
            [bodyData, nodalData, stressData] = readAnsys();
        case STRAND7
            readStrand7()
        otherwise
            return
    end

    generateData(bodyData, nodalData, stressData, general)

    function [bodyData, nodalData, stressData] = readAnsys()
        %readAnsys - ANSYS specific import
        %
        % Syntax: [output] = readAnsys(fileId)
        ANSYS = "ANSYS";
        pathSep = general.local.pathSep;
        [const] = loadConstants(ANSYS, pathSep);

        generateDirectories(general);

        [bodyData, nodalData] = readDsDat();
        [stressData] = readNodalSol();

        function [bodyData, nodalData] = readDsDat()
            %readDsDat - Imports data from ds.dat file.
            %
            % Syntax: varargout = readDsDat()
            fileName = general.local.dirs.simulationDir + const.files.ds;

            regex = const.regex;

            fileId = fopen(fileName,'rt');
            str = fgetl(fileId);
            prevStr = "";
            elementBlockCount = 1;
            bodyData{1}.elementBlockData = [];
            bodyData{1}.connectivity = [];
            bodyData{1}.elementIdx = [];
            nodalData.coordinates = [];
            nodalData.nodeIdx = [];

            while ~atEnd(str, regex.elemEnd)
                [match, nonMatch] = regexp(str, regex.block, 'names', 'split');
                if ~isempty(match)
                    [blockData, ~] = regexp(nonMatch{2}, regex.fields, 'names', 'split');
                    switch match.fieldType
                        case "n"
                            [nodalData.coordinates, nodalData.nodeIdx] = readNode();
                        case "e"
                            [elementBlockData, ~] = regexp(prevStr, regex.elementBlockData, 'names', 'split');
                            elementType = uint32(str2double(elementBlockData.elementType));
                            iBody = uint32(str2double(elementBlockData.iBody));
                            [bodyData{iBody}.connectivity, bodyData{iBody}.elementIdx, bodyData{iBody}.elementBlockData] = readElement();
                            elementBlockCount = elementBlockCount + 1;
                        case "end"
                            continue
                    end
                end
                prevStr = str;
                str = fgetl(fileId);
            end
            fclose(fileId);

            function [connectivity, elementIdx, elementBlockData] = readElement()
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

                general.local.body{iBody}.elementData(elementBlockCount,:) = [nElements, nodesPerElement, elementType];

                connectivity = zeros([nElements,nodesPerElement], 'uint32');
                elementIdx = zeros([nElements,1],'uint32');
                
                tmp = textscan(fileId, const.format.elements, 'MultipleDelimsAsOne',true, 'CollectOutput', true);%Read in all data
            
                connectivity(:, :) = [elementFields(12:12 + nodesPerElement - 1)'; tmp{:}(1:end-1, 12:12 + nodesPerElement - 1)];%Read to preallocated array
                elementIdx(:, :) = [elementFields(11); tmp{:}(1:end-1, 11)];

                iBodyiBlockNodes = unique(reshape(connectivity,1,[]));
                general.local.body{iBody}.elementData(elementBlockCount,:) = [nElements, nodesPerElement, elementType];
                elementBlockData = [nElements, nodesPerElement, elementType];


                elementPrepPath = general.local.dirs.workingDir + general.dirs.prepPathE + num2str(elementBlockCount) + general.files.matExt;

                dataLabels = ["elementIdx", "connectivity"];
                inputData.connectivity = connectivity;
                inputData.elementIdx = elementIdx;

                saveToMat(elementPrepPath, dataLabels, inputData, false);
            end

            function [coords, nodeIdx] = readNode()
                %readNode - Reads node coordinate data and transforms the data ready for saving
                %
                % Syntax: [dataLabels, inputData] = readNode()

                maxNodes = uint32(str2double(blockData(end).fields));

                fgetl(fileId); %increment the cursor to node block

                general.local.maxNodes = maxNodes;

                coords = NaN([maxNodes, 3], 'double');
                nodeIdx = zeros([maxNodes, 1], 'uint32');

                tmp = textscan(fileId, const.format.nodes,'MultipleDelimsAsOne', true, 'CollectOutput', true);
                endIdx = length(tmp{1, 1}(:));

                nodeIdx(1:endIdx) = tmp{1}(1:endIdx);
                coords(1:endIdx, :) = tmp{2:end}(1:endIdx, :);

                nodePrepPath = general.local.dirs.workingDir + general.dirs.prepPathN + general.files.nodes;
                dataLabels = ["nodeIdx", "coords"];
                inputData.nodeIdx = nodeIdx;
                inputData.coords = coords;

                saveToMat(nodePrepPath, dataLabels, inputData, false);
            end

        end

        function [stressData] = readNodalSol()
            %readNodalSol - reads in the stress data from nodalSolution.txt
            %
            % Syntax: varargout = readNodalSol()
            fileName = general.local.dirs.simulationDir + const.files.nodalSol;

            stressPrepPath = general.local.dirs.workingDir + general.dirs.prepPathS + general.files.stress;

            fileId = fopen(fileName,'rt');
            str = convertCharsToStrings(fgetl(fileId));
            lineTest = true;

            stress = zeros([general.local.maxNodes, 6], 'double');
            nodalStressIdx = zeros([general.local.maxNodes, 1], 'uint32');

            while lineTest
                [match, nonMatch] = regexp(str, const.regex.nodeSolBlock, 'names', 'noemptymatch');
                if ~isempty(match)
                    if length(match)>1 && ~isempty(match{1,2}(:))
                        tmp = textscan(fileId, const.format.nodalStress,'MultipleDelimsAsOne', true, 'CollectOutput', true);
                        lineTest = false;
                    end
                end
                str = convertCharsToStrings(fgetl(fileId));
            end

            endIdx = length(tmp{1}(1:end, 1));

            nodalStressIdx(1:endIdx) = tmp{1}(1:endIdx);
            stress(1:endIdx, :) = tmp{2:end}(1:endIdx, :);

            dataLabels = ["nodalStressIdx", "stress"];
            inputData.stress = stress;
            inputData.nodalStressIdx = nodalStressIdx;

            saveToMat(stressPrepPath, dataLabels, inputData, false);

            stressData.stress = stress;
            stressData.nodalStressIdx = nodalStressIdx;
        end

        dataLabels = ["general"];
        inputData.general = general;
        generalOutPath = general.local.dirs.workingDir + general.path.general;

        saveToMat(generalOutPath, dataLabels, inputData, false);
    end
end

function generateData(bodyData, nodalData, stressData, general)
    nBodies = size(bodyData{:}, 1);
    locatingData{nBodies}.elementFaceIndexingArray = [];
    locatingData{nBodies}.faceMappingVector = [];
    locatingData{nBodies}.faceCoords = [];
    locatingData{nBodies}.bodyHullIdx = [];
    locatingData{nBodies}.kdTreeElements = [];
    locatingData{nBodies}.rayCastVars = [];
    for k = 1:nBodies
        bodyData{k}.elementInterpFunc = cell(bodyData{k}.elementBlockData(1),1);

        % DISABLED DURING DEVELOPEMENT
        
        [locatingData{k}.elementFaceIndexingArray, ...
        locatingData{k}.faceMappingVector, ...
        locatingData{k}.faceCoords, ...
        locatingData{k}.bodyHullIdx] = generateFaces(bodyData{k}.elementBlockData, bodyData{k}.connectivity);
        
        [bodyData{k}.elementInterpFunc] = generateInterpFunction();

        locatingData{k}.kdTreeElements = getElementCentroids();
        % TODO - It is currently element 185 specific - make general. i.e.
        % I'm making triangles from the 4 nodes that define the face. 
        % Needs to be adapted for elements with more or less nodes.
        locatingData{k}.rayCastVars = generateLocatingVariables(general, [locatingData{k}.faceCoords([1:3], :,:), locatingData{k}.faceCoords([1,3,4], :,:)]);
    end
    inputData.bodyData = bodyData;
    inputData.locatingData = locatingData;
    dataLabels = ["bodyData","locatingData"];
    outPath = general.local.dirs.workingDir + general.path.utilities;
    saveToMat(outPath, dataLabels, inputData, false);
    function [elementFaceIndexingArray, faceMappingVector, faceCoords, bodyHullIdx] = generateFaces(elementData, connectivity)
            faceDefinition = faceDef(elementData(3));

            [nFaces, nNodesPerFace] = size(faceDefinition);
            faceArray =  zeros([elementData(1)*nFaces, nNodesPerFace], 'uint32');

            faceArray(:,:) = reshape(connectivity(:,faceDefinition),elementData(1)*nFaces, nNodesPerFace);

            [sortedFaceArray, ~] = sort(faceArray, 2);

            [~, ind, faceMappingVector] = unique(sortedFaceArray, 'rows');

            elementFaceIndexingArray = reshape(1:elementData(1)*nFaces, elementData(1), nFaces);

            bodyHullIdx = histcounts(faceMappingVector, length(ind)) < 2;

            filteredFaces = faceArray(ind,:);

            faceCoords = getFaceCoords(filteredFaces);
        end
    function [faceCoords] = getFaceCoords(faceArray)
        %getFaceCoords - Returns the coordinates of faces defined by node indicies. 
        %                Faces to be defined as either clockwise or counter clockwise.
        % Syntax: [faceCoords] = getFaceCoords(general, faceArray)
        [nRows, nCols] = size(faceArray);
        faceOrderVector = zeros(nRows*nCols, 1);
        faceOrderVector(:) = reshape(faceArray,1,nRows*nCols); % get all node from all faces in vector
        [uNodes, faceOrder, faceIdx] = unique(faceOrderVector); % remove duplicates for efficiency
        faceCoords = zeros(size(uNodes,1),3);
        faceCoords(:,:) = getCoords(nodalData, uNodes); %get the coords of the nodes
        faceCoords = faceCoords(faceIdx, :);
        faceCoords = reshape(faceCoords, nRows, nCols,3);
        faceCoords = permute(faceCoords, [2 1 3]);
    end
    function [interpolationFunctions] = generateInterpFunction()
        %generateInterpFunction - Description
        %
        % Syntax: [interpolationFunctions] = generateInterpFunction(input)
        %
        % Long description
        myCluster = parcluster('local');
        nWorkers = myCluster.NumWorkers;
        interpolationFunctions = cell(bodyData{k}.elementBlockData(1), 1);
        parfor (j = 1:bodyData{k}.elementBlockData(1), nWorkers)
            localNodes = bodyData{k}.connectivity(j,:);
            localNodeCoords = getCoords(nodalData, localNodes);
            localStress = stressData.stress(localNodes,:);
            interpolationFunctions{j} = interpFunction(localStress, localNodeCoords);
        end
    end
    function kdTreeElements = getElementCentroids()
        elementCentroids = zeros([bodyData{k}.elementBlockData(1), 3]);
        myCluster = parcluster('local');
        nWorkers = myCluster.NumWorkers;
        parfor (j = 1:bodyData{k}.elementBlockData(1), nWorkers)
            localNodes = bodyData{k}.connectivity(j,:);
            localNodeCoords = getCoords(nodalData, localNodes);
            elementCentroids(j,:) = mean(localNodeCoords,1);
        end
        kdTreeElements = createns(elementCentroids,'NSMethod','kdtree');
    end
end

function generateDirectories(general)
    %TODO:  Future imporvement to clear the directories prior to creation.
    fn = fieldnames(general.dirs(:));
    for k = 1:numel(fn)
        [a, b] = mkdir(general.local.dirs.workingDir + general.dirs.(fn{k}));
    end
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

function [retVal] = interpFunction(stress, nodalCoordinates)
    %Natural interpolation method is used to form a stress function to then
    %compute the paths.
    Fxx = scatteredInterpolant(nodalCoordinates(:,1), nodalCoordinates(:,2), nodalCoordinates(:,3), stress(:,1), 'natural');
    Fyy = scatteredInterpolant(nodalCoordinates(:,1), nodalCoordinates(:,2), nodalCoordinates(:,3), stress(:,2), 'natural');
    Fzz = scatteredInterpolant(nodalCoordinates(:,1), nodalCoordinates(:,2), nodalCoordinates(:,3), stress(:,3), 'natural');
    Fxy = scatteredInterpolant(nodalCoordinates(:,1), nodalCoordinates(:,2), nodalCoordinates(:,3), stress(:,4), 'natural');
    Fyz = scatteredInterpolant(nodalCoordinates(:,1), nodalCoordinates(:,2), nodalCoordinates(:,3), stress(:,5), 'natural');
    Fxz = scatteredInterpolant(nodalCoordinates(:,1), nodalCoordinates(:,2), nodalCoordinates(:,3), stress(:,6), 'natural');
    retVal = {Fxx, Fyy, Fzz, Fxy, Fyz, Fxz};
end
function [coords] = getCoords(nodalData, nodeNumber)
    %getCoords - returns coordinates of node number(s).
    %
    % Syntax: [coords] = getCoords(nodeNumber)
    [~, idx, ~] = intersect(nodalData.nodeIdx, nodeNumber);
    coords = nodalData.coordinates(idx, :);
end

function [rayCastVars] = generateLocatingVariables(general, faceVerticies)
    % generateLocatingVariables - generates and stores the preprocessed data for locating the point in mesh.
    % varargout: true if error in data, false otherwise
    %
    % Syntax: [varargout] = generateLocatingVariable(faceVerticies)
    V0 = squeeze(faceVerticies(1, :, :));
    V1 = squeeze(faceVerticies(2, :, :));
    V2 = squeeze(faceVerticies(3, :, :));
    u = V1 - V0;
    v = V2 - V0;
    n = cross(u, v, 2);

    if all(~n, 'all')
        rayCastVars  = false;
        return
    end

    uu = dot(u,u, 2);
    uv = dot(u,v, 2);
    vv = dot(v,v, 2);
    D = (uv).^2 - (uu .* vv);

    rayCastVars.V0 = V0;
    rayCastVars.u = u;
    rayCastVars.v = v;
    rayCastVars.n = n;
    rayCastVars.uu = uu;
    rayCastVars.uv = uv;
    rayCastVars.vv = vv;
    rayCastVars.D = D;
    rayCastVars.infPoint = general.constants.infPoint;
    rayCastVars.tol = general.constants.tol;
end