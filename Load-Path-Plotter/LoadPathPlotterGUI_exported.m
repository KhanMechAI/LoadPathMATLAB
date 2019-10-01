classdef LoadPathPlotterGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        LoadPathGeneratorUIFigure  matlab.ui.Figure
        GridLayout                 matlab.ui.container.GridLayout
        LeftPanel                  matlab.ui.container.Panel
        ModelDataPanel             matlab.ui.container.Panel
        ModelFileEditFieldLabel    matlab.ui.control.Label
        ModelFileEditField         matlab.ui.control.EditField
        modelDataButton            matlab.ui.control.Button
        stressDataButton           matlab.ui.control.Button
        StressFileEditFieldLabel   matlab.ui.control.Label
        StressFileEditField        matlab.ui.control.EditField
        PreProcessButton           matlab.ui.control.Button
        BuildModelButton           matlab.ui.control.Button
        SeedsPanel                 matlab.ui.container.Panel
        SeedFileEditFieldLabel     matlab.ui.control.Label
        SeedFileEditField          matlab.ui.control.EditField
        seedFileButton             matlab.ui.control.Button
        SeedTable                  matlab.ui.control.Table
        Button                     matlab.ui.control.Button
        Button_2                   matlab.ui.control.Button
        PathdataPanel              matlab.ui.container.Panel
        PathDataTable              matlab.ui.control.Table
        CalculateButton            matlab.ui.control.Button
        PlotButton                 matlab.ui.control.Button
        CenterPanel                matlab.ui.container.Panel
        UIAxes                     matlab.ui.control.UIAxes
        RightPanel                 matlab.ui.container.Panel
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end


    properties (Access = private)
        modelBuiltFlag = false;
        stressImportedFlag = false;
        calculatedFlag = false;
        preProcessedFlag = false; % Description
        modelDataImportedFlag = false; % Description
    end
    %GUI vars
    properties (Access = private)
        selectedRows = [];
        tmpSeedTable
    end
    
    %Regex
    properties (Access = private)
        regexElementBlockData = "(?:et,)(?<iBody>\d+),(?<elementType>\d+)";
        regexBlock = "(?<fieldType>n|e)block";
        regexFields = "(?<fields>\d+)|(?<solid>solid)";
        regexElemEnd = "/wb,(?<sectionKey>\w+),(?<context>\w+)";
        regexNodeSolBlock = ["\s+(?<sectionKey>NODE)\s+","(?<stressFields>\<S\w{1,2})\s+"];
    end
    
    %Formats
    properties (Access = private)
        formatNodes = '%u32%f32%f32%f32';
        formatElements = repmat('%u32',[1,19]);
        formatNodalStress = '%u32%f32%f32%f32%f32%f32%f32';
    end
    
    %Directories
    properties
        outPath = "_output_data/"; % Path for outputs
        prepPath = "_prep_data/"; % Path for all preprocessed data 
        db = "db.mat";
    end
    
    % Constants
    properties (Access = public)
        infPoint;
        tol = 0.0001;
        stressTensor = [1 4 6;...
                        4 2 5;...
                        6 5 3];
    end
    
    % LPP Vars
    properties (Access = public)
        modelData
        locatingData
        nodalData
        stressData
        kdTreeCentroids
        interpolationFunctions
        rayCast
        maxNodes = 0;
        nBodies = 0;
        paths
    end
    
    properties
        pathData
        pathHandles
        currentPath = 1;
        defaultPathData = [1000 0.5]
    end
    
    methods (Access = private)
        
        function getFile(app, caller, pathOutputHandle)
            getFileSettings = getFileHelper(app, caller);
            [fileName,path,~] = uigetfile(getFileSettings{1}, 'Select a File');
            if ~path
                pathOutputHandle.Value = "File path...";
            else
                pathOutputHandle.Value = string(path) + string(fileName);
            end
        end
        
        function results = getFileHelper(app, caller)
            switch caller
                case 1
                    results = {'*.dat','ANSYS Files (*.dat)'};
                case 2
                    results = {'*.txt','Stress File (*.txt)'};
                case 3
                    results = {'*.txt','Seed File (*.txt)'};
            end
        end
        
        function success = readDsDat(app)
            %readDsDat - Imports data from ds.dat file.
            %
            % Syntax: varargout = readDsDat()
            success = 1;
            fileName = string(app.ModelFileEditField.Value);

            fileId = fopen(fileName,'rt');
            str = fgetl(fileId);
            prevStr = "";
            elementBlockCount = 1;

            while ~atEnd(app, str, app.regexElemEnd)
                [match, nonMatch] = regexp(str, app.regexBlock, 'names', 'split');
                if ~isempty(match)
                    [blockData, ~] = regexp(nonMatch{2}, app.regexFields, 'names', 'split');
                    switch match.fieldType
                        case "n"
                            readNode();
                        case "e"
                            [elementBlockData, ~] = regexp(prevStr, app.regexElementBlockData, 'names', 'split');
                            elementType = uint32(str2double(elementBlockData.elementType));
                            iBody = uint32(str2double(elementBlockData.iBody));
                            readElement(iBody);
                            elementBlockCount = elementBlockCount + 1;
                        case "end"
                            continue
                    end
                end
                prevStr = str;
                str = fgetl(fileId);
            end
            fclose(fileId);
            app.nBodies = elementBlockCount - 1;
            if app.nBodies < 1
                success = 0;
            end

            function readElement(iBody)
                %readElement - Reads element connectivity data and transforms the data ready for saving
                %
                % Syntax: [dataLabels, inputData] = readElement()
                %TODO: Write function to handle elements with more than 8 nodes.


                [nLines] = getNLines(app, fileId, 2);

                elementFields = uint32(str2double(split(strtrim(nLines(2,:)))));

                nElements = uint32(str2double(blockData(end).fields));

                if ~isempty([blockData.solid])
                    nodesPerElement = elementFields(9);
                else
                    nodesPerElement = uint32(str2double(blockData(1).fields));
                end

                connectivity = zeros([nElements,nodesPerElement], 'uint32');
                elementIdx = zeros([nElements,1],'uint32');
                
                tmp = textscan(fileId, app.formatElements, 'MultipleDelimsAsOne',true, 'CollectOutput', true);%Read in all data
                connectivity(:, :) = [elementFields(12:12 + nodesPerElement - 1)'; tmp{:}(1:end-1, 12:12 + nodesPerElement - 1)];%Read to preallocated array
                elementIdx(:, :) = [elementFields(11); tmp{:}(1:end-1, 11)];

                app.modelData{iBody}.elementData = [nElements, nodesPerElement, elementType];
                app.modelData{iBody}.connectivity = connectivity;
                app.modelData{iBody}.elementIdx = elementIdx;

            end

            function [coords, nodeIdx] = readNode()
                %readNode - Reads node coordinate data and transforms the data ready for saving
                %
                % Syntax: [dataLabels, inputData] = readNode()

                app.maxNodes = uint32(str2double(blockData(end).fields));

                fgetl(fileId); %increment the cursor to node block

                coords = NaN([app.maxNodes, 3], 'double');
                nodeIdx = zeros([app.maxNodes, 1], 'uint32');

                tmp = textscan(fileId, app.formatNodes,'MultipleDelimsAsOne', true, 'CollectOutput', true);
                endIdx = length(tmp{1, 1}(:));

                nodeIdx(1:endIdx) = tmp{1}(1:endIdx);
                coords(1:endIdx, :) = tmp{2:end}(1:endIdx, :);
                
                [~, farthestNode] = max(sum(coords,2));
                
                app.infPoint = coords(farthestNode, :) + 1;

                app.nodalData.nodeIdx = nodeIdx;
                app.nodalData.coordinates = coords;
            end
        end
                
        function success = readNodalSol(app)
            %readNodalSol - reads in the stress data from nodalSolution.txt
            %
            % Syntax: varargout = readNodalSol()
            success = 1;
            fileName = string(app.StressFileEditField.Value);
        
            fileId = fopen(fileName,'rt');
            str = convertCharsToStrings(fgetl(fileId));
            lineTest = true;
        
            stress = zeros([app.maxNodes, 6], 'double');
            nodalStressIdx = zeros([app.maxNodes, 1], 'uint32');
            tmp = {};
            while lineTest
                [match, ~] = regexp(str, app.regexNodeSolBlock, 'names', 'noemptymatch');
                if ~isempty(match)
                    if length(match)>1 && ~isempty(match{1,2}(:))
                        tmp = textscan(fileId, app.formatNodalStress,'MultipleDelimsAsOne', true, 'CollectOutput', true);
                        lineTest = false;
                    end
                end
                str = convertCharsToStrings(fgetl(fileId));
            end
            if isempty(tmp)
                success = 0;
                return
            end
            endIdx = length(tmp{1}(1:end, 1));
            nodalStressIdx(1:endIdx) = tmp{1}(1:endIdx);
            stress(1:endIdx, :) = tmp{2:end}(1:endIdx, :);
            
            app.stressData.stress = stress;
            app.stressData.nodalStressIdx = nodalStressIdx;
        end
        
        function endElements = atEnd(~, str, pat)
    
            %atEnd - Returns true if at the end of the elements. False otherwise.
            %
            % Syntax: endElements = atEnd(str)
            endElements = false;
        
            [match, ~] = regexp(str, pat, 'names', 'split');
        
            if ~isempty(match)
                if strcmp(match.sectionKey, "elem")
                    if strcmp(match.context, "end")
                        endElements = true;
                    end
                end
            end
        end
        
        function [nLines] = getNLines(~, fileId, n)
            %getNLines - Iterates over 'n' lines of a file returning each in a cell of a cell array
            %
            % Syntax: [nLines] = getNLines(input)
            nLines = strings(n,1);
            for k=1:n
                nLines(k,:) = fgetl(fileId);
            end            
        end
        
        function generateFaces(app)
            for k = 1:app.nBodies 
                elementData = app.modelData{k}.elementData;
                connectivity = app.modelData{k}.connectivity;
                
                faceDefinition = faceDef(app, app.modelData{k}.elementData(3));
    
                [nFaces, nNodesPerFace] = size(faceDefinition);
                faceArray =  zeros([elementData(1)*nFaces, nNodesPerFace], 'uint32');
    
                faceArray(:,:) = reshape(connectivity(:,faceDefinition),elementData(1)*nFaces, nNodesPerFace);
    
                [sortedFaceArray, ~] = sort(faceArray, 2);
    
                [~, ind, app.locatingData{k}.faceMappingVector] = unique(sortedFaceArray, 'rows');
    
                app.locatingData{k}.elementFaceIndexingArray = reshape(1:elementData(1)*nFaces, elementData(1), nFaces);
    
                app.locatingData{k}.hullIdx = histcounts(app.locatingData{k}.faceMappingVector, length(ind)) < 2;
    
                filteredFaces = faceArray(ind,:);
    
                app.locatingData{k}.faceCoords = getFaceCoords(app, filteredFaces);
            end
        end
        
        function [faceCoords] = getFaceCoords(app, faceArray)
            %getFaceCoords - Returns the coordinates of faces defined by node indicies. 
            %                Faces to be defined as either clockwise or counter clockwise.
            % Syntax: [faceCoords] = getFaceCoords(general, faceArray)
            [nRows, nCols] = size(faceArray);
            faceOrderVector = zeros(nRows*nCols, 1);
            faceOrderVector(:) = reshape(faceArray,1,nRows*nCols); % get all node from all faces in vector
            [uNodes, ~, faceIdx] = unique(faceOrderVector); % remove duplicates for efficiency
            faceCoords = zeros(size(uNodes,1),3);
            faceCoords(:,:) = getCoords(app, uNodes); %get the coords of the nodes
            faceCoords = faceCoords(faceIdx, :);
            faceCoords = reshape(faceCoords, nRows, nCols,3);
            faceCoords = permute(faceCoords, [2 1 3]);
        end
        function faceCoords = getBodyFaceCoords(app, body)
           faceCoords = app.locatingData{body}.faceCoords(:,app.locatingData{body}.hullIdx,:);
        end
        function [coords] = getCoords(app, nodeNumber)
            %getCoords - returns coordinates of node number(s).
            %
            % Syntax: [coords] = getCoords(nodeNumber)
            [~, idx, ~] = intersect(app.nodalData.nodeIdx, nodeNumber);
            coords = app.nodalData.coordinates(idx, :);
        end
        
        function [faces] = faceDef(~, elementType)
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
        
       
        function generateData(app)
            %generateInterpFunction - Description
            %
            % Syntax: [interpolationFunctions] = generateInterpFunction(input)

            for k = 1:app.nBodies
                interpolationFunctions = cell(app.modelData{k}.elementData(1), 1);
                elementCentroids = zeros([app.modelData{k}.elementData(1), 3]);
                connectivity = app.modelData{k}.connectivity;
                elementData = app.modelData{k}.elementData;

                for j = 1:app.modelData{k}.elementData(1)
                    localNodes = connectivity(j,:);
                    localNodeCoords = getCoords(app, localNodes);
                    localStress = app.stressData.stress(localNodes,:);
                    interpolationFunctions{j} = interpFunction(app, localStress, localNodeCoords);
                    elementCentroids(j,:) = mean(localNodeCoords,1);
                end
%                 parfor (j = 1:elementData(1), nWorkers)
%                     localNodes = connectivity(j,:);
%                     localNodeCoords = getCoords(app, localNodes);
%                 end
%                 parfor (j = 1:elementData(1), nWorkers)
%                     localNodes = app.modelData{k}.connectivity(j,:);
%                     localNodeCoords = getCoords(app, localNodes);
%                     localStress = app.stressData.stress(localNodes,:);
%                     interpolationFunctions{j} = interpFunction(app, localStress, localNodeCoords);
%                     elementCentroids(j,:) = mean(localNodeCoords,1);
%                 end
                %TODO - Make more general for other element types
                %TODO - put error catch logic for degenerate faces
                app.interpolationFunctions{k} = interpolationFunctions;
                app.rayCast{k} = generateLocatingVariables(app, [app.locatingData{k}.faceCoords(1:3, :,:), app.locatingData{k}.faceCoords([1,3,4], :,:)]);
                app.kdTreeCentroids{k} = createns(elementCentroids,'NSMethod','kdtree');
            end
        end
        
        function [retVal] = interpFunction(~, stress, nodalCoordinates)
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
        
        function [rayCast] = generateLocatingVariables(app, faceVerticies)
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
                rayCast = false;
                return
            end
        
            uu = dot(u,u, 2);
            uv = dot(u,v, 2);
            vv = dot(v,v, 2);
            D = (uv).^2 - (uu .* vv);
        
            rayCast.V0 = V0;
            rayCast.u = u;
            rayCast.v = v;
            rayCast.n = n;
            rayCast.uu = uu;
            rayCast.uv = uv;
            rayCast.vv = vv;
            rayCast.D = D;
            rayCast.infPoint = app.infPoint;
            rayCast.tol = app.tol;
        end
        
        function loadSeedsToTable(app)
            fileName = string(app.SeedFileEditField.Value);
            tempTable = readtable(fileName);
            tempTable.Properties.VariableNames = {'X' 'Y' 'Z' };
            Active = logical(ones(size(tempTable,1), 1));
            if ~isempty(app.SeedTable.Data.Seed)
                lastSeed = app.SeedTable.Data.Seed(end)+1;
            else
               lastSeed = 1; 
            end
            
            Seed = int32([lastSeed:lastSeed+size(tempTable(:,1),1)-1])';
            if isempty(Seed)
                return
            end
            tempTable.Active = Active;
            tempTable.Seed = Seed;
            app.SeedTable.Data = [app.SeedTable.Data; tempTable];
        end

        function activatePreProcess(app)
            if app.modelBuiltFlag && app.stressImportedFlag
                app.PreProcessButton.Enable = true;
            end
        end
        function updatePathData(app)
            app.pathData(app.currentPath,:) = [app.PathLengthEditField.Value, app.StepSizeEditField.Value];     
            app.pathData
        end
        
        function togglePath(app, pathDir)
            switch pathDir
                case 1
                    app.XButton.Enable = app.XToggle.Value;
                    value = app.XToggle.Value;
                    if value && app.currentPath == 1
                        changePathDataDisplay(app, true);
                    elseif ~value && app.currentPath == 1
                        changePathDataDisplay(app, false);
                    end
                case 2
                    app.YButton.Enable = app.YToggle.Value;
                    value = app.YToggle.Value;
                    if value && app.currentPath == 2
                        changePathDataDisplay(app, true);
                    elseif ~value && app.currentPath == 2
                        changePathDataDisplay(app, false);
                    end
                case 3
                    app.ZButton.Enable = app.ZToggle.Value;
                    value = app.ZToggle.Value;
                    if value && app.currentPath == 3
                        changePathDataDisplay(app, true);
                    elseif ~value && app.currentPath == 3
                        changePathDataDisplay(app, false);
                    end
            end
        end
        
        function changePathDataInput(app, pathToChange)
            switch pathToChange
                case 1
                    app.pathDataLabel.Text = 'X path data';
                case 2
                    app.pathDataLabel.Text = 'Y path data';
                case 3
                    app.pathDataLabel.Text = 'Z path data';
            end
        end
        function changePathDataDisplay(app, toggle)
            app.PathLengthEditField.Enable = toggle;
            app.StepSizeEditField.Enable = toggle;
        end
    end
    
    methods (Access = public)
        function plotPoints(app, point)
            hold(app.UIAxes,"on")
            plot3(app.UIAxes, point.X, point.Y, point.Z, '+', 'MarkerFaceColor', 'red')
            hold(app.UIAxes,"off")
        end
        
        function [minIntensity, maxIntensity] = calculateMinMaxIntensity(app)
            minIntensity = inf;
            maxIntensity = -inf;
            for d = [1 2 3]
                if ~app.PathDataTable.Data.Active(d)
                    continue
                end
                
                tmpMax = max(max(app.paths{d}(:,4,:)));
                tmpMin = min(min(app.paths{d}(:,4,:)));
                
                if tmpMin < minIntensity
                    minIntensity = tmpMin;
                end
                
                if tmpMax > maxIntensity
                    maxIntensity = tmpMax;
                end
            end
        end
        
        function togglePaths(app)
            for d = [1 2 3]
                if ~app.PathDataTable.Data.Active(d)
                    continue
                end
                for k = 1:app.SeedTable.Data.Seed(end)
                    line = app.pathHandles{d}{k};
                    if app.SeedTable.Data.Active(k)
                        line.Visible = 'on';
                    else
                        line.Visible = 'off';
                    end
                end
            end
        end
        
        function plotPaths(app)
            pathHandles = cell([1, 3]);

            for d = [1 2 3]
                if ~app.PathDataTable.Data.Active(d)
                    continue
                end
                nPaths = size(app.paths{d},3);
                paths = reshape(permute(app.paths{d},[1 3 2]), size(app.paths{d},1)*nPaths,[]);
                pathLength = app.PathDataTable.Data.PathLength(d);

                pathHandles{d} = cell([nPaths,1]);

                for k=1:nPaths
                    startFirst = ((k-1)*pathLength) + 1;
                    finishFirst =  ((k-1)*pathLength)+pathLength/2;
                    finishSecond = k*pathLength;
                    firstNan = find(~isnan(paths(startFirst:finishFirst,1)),1);
                    lastNan = find(isnan(paths(finishFirst+1:finishSecond,1)),1);
                    
                    pathHandles{d}{k} = patch(app.UIAxes,...
                                            paths(startFirst+firstNan+1:finishFirst+lastNan,1),...
                                            paths(startFirst+firstNan+1:finishFirst+lastNan,2),...
                                            paths(startFirst+firstNan+1:finishFirst+lastNan,3),...
                                            paths(startFirst+firstNan+1:finishFirst+lastNan,4),...
                                            'EdgeColor',...
                                            'interp'...
                                            );
                end
            end
            cb = colorbar(app.UIAxes);
            cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), 10);
            app.pathHandles = pathHandles;
            togglePaths(app)
        end
        
        function plotBody(app, varargin)
            if nargin == 1 
                for k = 1:app.nBodies
                    bodyFaceCoords = getBodyFaceCoords(app, k);
                    plotFaces(app, bodyFaceCoords, 'black', 0.2);
                end
            else
                bodyFaceCoords = getBodyFaceCoords(app, varargin{1});
                plotFaces(app, bodyFaceCoords, 'black', 0.2);
            end
        end
        
        function plotFaces(app, faceArray, colour, alp)
            %plotFaces - Description
            %
            % Syntax: plotFaces(faceArray)
            hold(app.UIAxes,"on")
            patch(app.UIAxes,'XData',faceArray(:, :, 1),'YData',faceArray(:, :, 2),'ZData',faceArray(:, :, 3), 'EdgeColor',colour,'FaceColor','none', 'EdgeAlpha', alp)
            hold(app.UIAxes,"off")
        end
        
        function plotSeeds(app)
            seeds = app.SeedTable.Data(app.SeedTable.Data.Active,{'X','Y', 'Z'});
            plotPoints(app, seeds)
        end
        
        function plotRay(app, varargin)
            point = zeros([2,3]);
            point(1,:) = [1e2 1e2 1e2];
            if nargin > 1
                point(2,:) = varargin{1};
            else
                point(2,:) = p0;
            end
            hold(app.UIAxes,"on")
            plot3(app.UIAxes, point(:,1), point(:,2), point(:,3))
            hold(app.UIAxes,"off")
        end
        function buttonInteractivityHelper(app)
            if ~app.modelDataImportedFlag
                buttonEnabledArray = [false, false, false, false];
            elseif ~app.modelBuiltFlag
                buttonEnabledArray = [true, false, false, false];
            elseif ~app.stressImportedFlag
                buttonEnabledArray = [true, false, false, false];
            elseif ~app.preProcessedFlag
                buttonEnabledArray = [true, true, false, false];
            elseif size(app.SeedTable.Data(:,1), 1) < 1
                buttonEnabledArray = [true, true, false, false];
            elseif ~app.calculatedFlag
                buttonEnabledArray = [true, true, true, false];
            else
                buttonEnabledArray = [true, true, true, true];
            end
            
            buttonArray = {app.BuildModelButton,...
                           app.PreProcessButton,...
                           app.CalculateButton,...
                           app.PlotButton,...
            };
            for k = 1:numel(buttonArray)
                button = buttonArray{k};
                button.Enable = buttonEnabledArray(k);
            end
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            Seed = 1;
            X = 0;
            Y = 0;
            Z = 0;
            Active = true;
            app.tmpSeedTable = table(Seed, X, Y, Z, Active);
            app.SeedTable.Data = app.tmpSeedTable;
            app.UIAxes.DataAspectRatio = [1 1 1];
                       
            Path = {'X', 'Y', 'Z'}';
            PathLength = [1000 1000 1000]';
            StepSize = [0.5 0.5 0.5]';
            Active = [true false false]';
            Calculated = [false false false]';
            app.PathDataTable.Data = table(Path, PathLength, StepSize, Active, Calculated);
        end

        % Button pushed function: stressDataButton
        function stressDataButtonPushed(app, event)
            getFile(app, 2, app.StressFileEditField);
            app.stressImportedFlag = true;
            app.preProcessedFlag = false;
            app.calculatedFlag = false;
            buttonInteractivityHelper(app)
        end

        % Button pushed function: seedFileButton
        function seedFileButtonPushed(app, event)

            getFile(app, 3, app.SeedFileEditField);
            loadSeedsToTable(app)
            plotSeeds(app)
            app.calculatedFlag = false;
            buttonInteractivityHelper(app)

        end

        % Button pushed function: modelDataButton
        function modelDataButtonPushed(app, event)
            getFile(app, 1, app.ModelFileEditField);
            app.modelDataImportedFlag = true;
            app.modelBuiltFlag = false;
            app.stressImportedFlag = false;
            app.preProcessedFlag = false;
            app.calculatedFlag = false;
            buttonInteractivityHelper(app)
            
        end

        % Button pushed function: Button
        function addRowsToTable(app, event)
            oldData = app.SeedTable.Data;
%             app.seedCounter = app.seedCounter + 1;
            if ~isempty(oldData.Seed)
                newData = {oldData.Seed(end)+1, 0, 0, 0, true};
            else
                newData = {1, 0, 0, 0, true};
            end
            app.SeedTable.Data = [oldData; newData];
            app.calculatedFlag = false;
            buttonInteractivityHelper(app)
        end

        % Button pushed function: Button_2
        function removeRowsFromTable(app, event)
            rowsToDelete = app.selectedRows;
            if ~isempty(rowsToDelete)
                app.SeedTable.Data(rowsToDelete,:) = [];
                nSeedsToAdd = size(app.SeedTable.Data(:,1),1);
                app.SeedTable.Data(:,1) = table(int32([1:nSeedsToAdd]'));
                app.selectedRows = [];
            end
            buttonInteractivityHelper(app)
        end

        % Cell selection callback: SeedTable
        function SeedTableCellSelection(app, event)
            if ~strcmp(event.EventName, 'ButtonPushed')
                app.selectedRows = unique(event.Indices(:,1));
            end
        end

        % Button pushed function: BuildModelButton
        function BuildModelButtonPushed(app, event)
            success = readDsDat(app);
            if success
                app.locatingData{app.nBodies}.elementFaceIndexingArray = [];
                app.locatingData{app.nBodies}.faceMappingVector = [];
                app.locatingData{app.nBodies}.faceCoords = [];
                app.locatingData{app.nBodies}.hullIdx = [];
                generateFaces(app);
                app.UIAxes.cla();
                plotBody(app)
                plotSeeds(app)
                app.modelBuiltFlag = true;
                activatePreProcess(app)
            else
                app.modelBuiltFlag = false;
                errordlg('File Error','Empty model, check your model data file carefully.');
            end
            buttonInteractivityHelper(app)
        end

        % Button pushed function: PreProcessButton
        function PreProcessButtonPushed(app, event)

            success = readNodalSol(app);
            if success
                app.kdTreeCentroids{app.nBodies} = [];
                app.interpolationFunctions{app.nBodies} = {};
                app.rayCast{app.nBodies} = {};
                generateData(app);
                app.preProcessedFlag = true;
            else
                app.preProcessedFlag = false;
                errordlg('File Error','Input data error, check your stress and model data carefully.');
            end
            buttonInteractivityHelper(app)
        end

        % Callback function
        function XToggleValueChanged(app, event)
            togglePath(app, 1);
        end

        % Callback function
        function YToggleValueChanged(app, event)
           togglePath(app, 2);
        end

        % Callback function
        function ZToggleValueChanged(app, event)
            togglePath(app, 3);
        end

        % Callback function
        function PathLengthEditFieldValueChanged(app, event)
            updatePathData(app)
        end

        % Callback function
        function StepSizeEditFieldValueChanged(app, event)
            updatePathData(app)
        end

        % Callback function
        function XButtonPushed(app, event)
            app.currentPath = 1;
            changePathDataDisplay(app, true)
            changePathDataInput(app, 1);
        end

        % Callback function
        function YButtonPushed(app, event)
            app.currentPath = 2;
            changePathDataDisplay(app, true)
            changePathDataInput(app, 2);
        end

        % Callback function
        function ZButtonPushed(app, event)
            app.currentPath = 3;
            changePathDataDisplay(app, true)
            changePathDataInput(app, 3);
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            % put some condition checking code in first
            app.paths = CalculateLoadPaths(app, app.SeedTable.Data);
            app.calculatedFlag = true;
            buttonInteractivityHelper(app)
            plotPaths(app);
        end

        % Display data changed function: SeedTable
        function SeedTableDisplayDataChanged(app, event)
            newDisplayData = app.SeedTable.DisplayData;
        end

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
            app.UIAxes.cla();
            plotBody(app)
            plotSeeds(app)
            plotPaths(app)
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.LoadPathGeneratorUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {776, 776, 776};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {776, 776};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {464, '1x', 12};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create LoadPathGeneratorUIFigure and hide until all components are created
            app.LoadPathGeneratorUIFigure = uifigure('Visible', 'off');
            app.LoadPathGeneratorUIFigure.AutoResizeChildren = 'off';
            app.LoadPathGeneratorUIFigure.Position = [100 100 1360 776];
            app.LoadPathGeneratorUIFigure.Name = 'Load Path Generator';
            app.LoadPathGeneratorUIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.LoadPathGeneratorUIFigure);
            app.GridLayout.ColumnWidth = {464, '1x', 12};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create ModelDataPanel
            app.ModelDataPanel = uipanel(app.LeftPanel);
            app.ModelDataPanel.Title = 'Model Data';
            app.ModelDataPanel.Position = [11 585 440 160];

            % Create ModelFileEditFieldLabel
            app.ModelFileEditFieldLabel = uilabel(app.ModelDataPanel);
            app.ModelFileEditFieldLabel.Position = [13 94 61 22];
            app.ModelFileEditFieldLabel.Text = 'Model File';

            % Create ModelFileEditField
            app.ModelFileEditField = uieditfield(app.ModelDataPanel, 'text');
            app.ModelFileEditField.FontColor = [0.8 0.8 0.8];
            app.ModelFileEditField.Position = [90 94 183 22];
            app.ModelFileEditField.Value = 'File path...';

            % Create modelDataButton
            app.modelDataButton = uibutton(app.ModelDataPanel, 'push');
            app.modelDataButton.ButtonPushedFcn = createCallbackFcn(app, @modelDataButtonPushed, true);
            app.modelDataButton.Position = [290 94 100 22];
            app.modelDataButton.Text = 'Choose File...';

            % Create stressDataButton
            app.stressDataButton = uibutton(app.ModelDataPanel, 'push');
            app.stressDataButton.ButtonPushedFcn = createCallbackFcn(app, @stressDataButtonPushed, true);
            app.stressDataButton.Position = [290 59 100 22];
            app.stressDataButton.Text = 'Choose File...';

            % Create StressFileEditFieldLabel
            app.StressFileEditFieldLabel = uilabel(app.ModelDataPanel);
            app.StressFileEditFieldLabel.Position = [13 59 62 22];
            app.StressFileEditFieldLabel.Text = 'Stress File';

            % Create StressFileEditField
            app.StressFileEditField = uieditfield(app.ModelDataPanel, 'text');
            app.StressFileEditField.FontColor = [0.8 0.8 0.8];
            app.StressFileEditField.Position = [90 59 183 22];
            app.StressFileEditField.Value = 'File path...';

            % Create PreProcessButton
            app.PreProcessButton = uibutton(app.ModelDataPanel, 'push');
            app.PreProcessButton.ButtonPushedFcn = createCallbackFcn(app, @PreProcessButtonPushed, true);
            app.PreProcessButton.Enable = 'off';
            app.PreProcessButton.Position = [231 8 100 22];
            app.PreProcessButton.Text = 'Pre Process';

            % Create BuildModelButton
            app.BuildModelButton = uibutton(app.ModelDataPanel, 'push');
            app.BuildModelButton.ButtonPushedFcn = createCallbackFcn(app, @BuildModelButtonPushed, true);
            app.BuildModelButton.Enable = 'off';
            app.BuildModelButton.Position = [89 9 100 22];
            app.BuildModelButton.Text = 'Build Model';

            % Create SeedsPanel
            app.SeedsPanel = uipanel(app.LeftPanel);
            app.SeedsPanel.Title = 'Seeds';
            app.SeedsPanel.Position = [11 227 441 343];

            % Create SeedFileEditFieldLabel
            app.SeedFileEditFieldLabel = uilabel(app.SeedsPanel);
            app.SeedFileEditFieldLabel.Position = [9 275 56 22];
            app.SeedFileEditFieldLabel.Text = 'Seed File';

            % Create SeedFileEditField
            app.SeedFileEditField = uieditfield(app.SeedsPanel, 'text');
            app.SeedFileEditField.FontColor = [0.8 0.8 0.8];
            app.SeedFileEditField.Position = [89 275 191 22];
            app.SeedFileEditField.Value = 'File path...';

            % Create seedFileButton
            app.seedFileButton = uibutton(app.SeedsPanel, 'push');
            app.seedFileButton.ButtonPushedFcn = createCallbackFcn(app, @seedFileButtonPushed, true);
            app.seedFileButton.Position = [299 275 100 22];
            app.seedFileButton.Text = 'Choose File...';

            % Create SeedTable
            app.SeedTable = uitable(app.SeedsPanel);
            app.SeedTable.ColumnName = {'Seeds'; 'X'; 'Y'; 'Z'; 'Active'};
            app.SeedTable.ColumnWidth = {60, 60, 60, 60};
            app.SeedTable.RowName = {''};
            app.SeedTable.ColumnSortable = true;
            app.SeedTable.ColumnEditable = true;
            app.SeedTable.CellSelectionCallback = createCallbackFcn(app, @SeedTableCellSelection, true);
            app.SeedTable.DisplayDataChangedFcn = createCallbackFcn(app, @SeedTableDisplayDataChanged, true);
            app.SeedTable.Position = [6 17 371 234];

            % Create Button
            app.Button = uibutton(app.SeedsPanel, 'push');
            app.Button.ButtonPushedFcn = createCallbackFcn(app, @addRowsToTable, true);
            app.Button.Position = [391 231 27 22];
            app.Button.Text = '+';

            % Create Button_2
            app.Button_2 = uibutton(app.SeedsPanel, 'push');
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app, @removeRowsFromTable, true);
            app.Button_2.Position = [391 201 27 22];
            app.Button_2.Text = '-';

            % Create PathdataPanel
            app.PathdataPanel = uipanel(app.LeftPanel);
            app.PathdataPanel.Title = 'Path data';
            app.PathdataPanel.Position = [11 50 440 166];

            % Create PathDataTable
            app.PathDataTable = uitable(app.PathdataPanel);
            app.PathDataTable.ColumnName = {'Path'; 'Path  length'; 'Step size'; 'Active'; 'Calculated'};
            app.PathDataTable.RowName = {};
            app.PathDataTable.ColumnEditable = [false true true true false];
            app.PathDataTable.Position = [9 12 420 121];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.LeftPanel, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Enable = 'off';
            app.CalculateButton.Position = [84 17 100 22];
            app.CalculateButton.Text = 'Calculate';

            % Create PlotButton
            app.PlotButton = uibutton(app.LeftPanel, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.Enable = 'off';
            app.PlotButton.Position = [261 17 99 23];
            app.PlotButton.Text = 'Plot';

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.CenterPanel);
            title(app.UIAxes, 'Review Model')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.DataAspectRatio = [1 1 1];
            app.UIAxes.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes.Position = [6 6 854 751];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Show the figure after all components are created
            app.LoadPathGeneratorUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = LoadPathPlotterGUI_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.LoadPathGeneratorUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.LoadPathGeneratorUIFigure)
        end
    end
end