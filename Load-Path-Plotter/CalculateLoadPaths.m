function allPaths = CalculateLoadPaths(app, SeedTable)
    %CalculateLoadPaths - Description
    %
    % Syntax: paths = CalculateLoadPaths(input)
    
    %p0 is initial seed point. Projection multiplier is used to 'jump' gaps
    %between parts in the model. It will project the path from onee part to
    %another. This procedure is a big source of time, some thought needs to
    %be given to how to optimise this routine.

    
    % [general] = getData(general, dataClass, varargin)
    activePaths = [1 2 3];
    allPaths = cell([3,1]);
    for pathDim = activePaths(app.PathDataTable.Data.Active)
        switch pathDim
            case 1
                V =@(stress, shearxy,shearxz) [stress, shearxy, shearxz];
            case 2
                V =@(stress, shearxy, shearyz) [shearxy, stress, shearyz];
            case 3
                V =@(stress, shearxz, shearyz) [shearxz, shearyz, stress];
        end
        pathLength = app.PathDataTable.Data.PathLength(pathDim);
        stepSize = app.PathDataTable.Data.StepSize(pathDim);
        paths = nan(pathLength,4,size(SeedTable,1));
        for seed=1:size(SeedTable,1)
            reversePath = false;
            p0 = [SeedTable.X(seed) SeedTable.Y(seed) SeedTable.Z(seed)];
            [startElement, startElementUtilities, startBody] = findCurrentElement(app, p0);
            if ~startElement
                fprintf('Seed %1$ not in model.', seed)
                continue
            end
            [F, Fs1, Fs2] = getInterpolationData(startElement, app.interpolationFunctions{startBody}, pathDim, reversePath);

            w = 1; %Start from center of array
            prevElement = startElement;
            currElementUtilities = startElementUtilities;
            currBody = startBody;
            while w < pathLength
                if w < pathLength/2
                    iterator = floor(pathLength/2) + w;
                else
                    iterator = pathLength - w;
                end
                try
                    paths(iterator, 1:3, seed) = p0;
                catch
                    aaa=1;
                end
                stress = F(p0(1), p0(2), p0(3));
                shear1 = Fs1(p0(1), p0(2), p0(3));
                shear2 = Fs2(p0(1), p0(2), p0(3));
                paths(iterator, 4, seed) = norm([stress shear1 shear2]);
                %Find dp1 at first interpolation.
                %Normalise the poining vector relative to the initial test point.
                %Calculate new points:

                %Runge-Kutta
                dp1 = V(stress, shear1, shear2)*stepSize/paths(iterator, 4, seed);
                p1 = p0 + dp1;

                dp2 = stressInterp(p1);
                p2 = p0 + 0.5*dp2;

                dp3 = stressInterp(p2);
                p3 = p0 + 0.5*dp3;

                dp4 = stressInterp(p3);

                p0 = p0 + 1/6 * (dp1 + 2*dp2 + 2*dp3 + dp4);

                %Locate element that the point is inside
                if ~triIntersect(p0, currElementUtilities)
                    [element, currElementUtilities, currBody] = findCurrentElement(app, p0);
                    if ~element
                        if reversePath
                            %Termninates loop
                            w = pathLength;
                        else
                            w = ceil(pathLength/2);
                            currBody = startBody;
                            currElementUtilities = startElementUtilities;
                            element = startElement;
                            reversePath = true;
                            p0 = [SeedTable.X(seed) SeedTable.Y(seed) SeedTable.Z(seed)];
                            [F, Fs1, Fs2] = getInterpolationData(element, app.interpolationFunctions{currBody}, pathDim, reversePath);
                        end
                    end
                end

                if element && prevElement ~= element
                    [F, Fs1, Fs2] = getInterpolationData(element, app.interpolationFunctions{currBody}, pathDim, reversePath);
                end
                w=w+1;
            end
        end
        allPaths{pathDim} = paths;
    end

    function [varargout] = findCurrentElement(app, p0, varargin)
        kNN = 8;
        body = 1;
        element = false;
        inElement = false;
        varargout = {element, [], body};
        if nargin >2
            element = varargin{1};
            body = varargin{2};
            elementUtilities = getUtilities(app.locatingData{body}, app.rayCast{body}, 'element', element);
            inElement = triIntersect(p0, elementUtilities);
            if inElement
                varargout = {element, elementUtilities, body};
            end
        end

        while body < app.nBodies + 1 && ~inElement
            testElements = knnsearch(app.kdTreeCentroids{body}, p0, 'K', kNN);
            k = 1;
            while k < kNN && ~inElement
                elementUtilities = getUtilities(app.locatingData{body}, app.rayCast{body}, 'element', testElements(k));
                inElement = triIntersect(p0, elementUtilities);
                if ~inElement
                    k = k +1;
                end
            end
            if ~inElement
                body = body + 1;
            else
                element = testElements(k);
                varargout = {element, elementUtilities, body};
            end
        end
    end
    function [dPoint] = stressInterp(p)
        if isnan(p)
           aaaa=1; 
        end
        stress = F(p(1), p(2), p(3));
        shear1 = Fs1(p(1), p(2), p(3));
        shear2 = Fs2(p(1), p(2), p(3));
        dPoint =  V(stress, shear1, shear2)*stepSize/paths(iterator, 4, seed);
        if isnan(dPoint)
           aaaa=1; 
        end
    end
    function plotElement(varargin)
        row=[];
        if nargin > 0
            element = varargin{:};
        end
        indexingArray = app.locatingData{currBody}.elementFaceIndexingArray;
        mappingVector = app.locatingData{currBody}.faceMappingVector;
        elidx = mappingVector(indexingArray(element,:));

        elementFaceCoords = app.locatingData{currBody}.faceCoords(:,elidx,:);
        app.plotFaces(elementFaceCoords, 'red', 1);
    end
    function plotRay(varargin)
        point = zeros([2,3]);
        point(1,:) = app.infPoint;
        if nargin > 0
            point(2,:) = varargin{1};
        else
            point(2,:) = p0;
        end
        gcf;
        hold(app.UIAxes,"on")
        plot3(app.UIAxes, point(:,1), point(:,2), point(:,3))
        hold(app.UIAxes,"off")
    end
end

function [F, Fs1, Fs2] = getInterpolationData(element, interpFunctions, pathDim, reversePath)
    %getInterpolationData - Description
    %
    % Syntax: [varagout] = getInterpolationData(point, general)
    switch pathDim
        case 1
            interpolantIdx = [1 4 6]; %Fxx Fxy Fxz
        case 2
            interpolantIdx = [2 4 5]; %Fyy Fxy Fyz
        case 3
            interpolantIdx = [3 6 5]; %Fzz Fyz Fxz
    end
    [F, Fs1, Fs2] = interpFunctions{element}{interpolantIdx};
    if reversePath
        F=@(x, y, z) -F(x, y, z);
        Fs1=@(x, y, z) -Fs1(x, y, z);
        Fs2=@(x, y, z) -Fs2(x, y, z);
    end
end


function [returnUtilities] = getUtilities(locatingData, rayCast, utilityType, varargin)
    %getUtilities - Description
    %
    % Syntax: [elementUtilities] = getUtilities(elementNo, utilities)
    mappingVector = locatingData.faceMappingVector;

    switch utilityType
        case 'body'
            idx = locatingData.hullIdx;
        case 'element'
            idx = locatingData.elementFaceIndexingArray(varargin{1},:);
    end

    idx = [mappingVector(idx)];
    idx = [idx; idx + size(rayCast.V0, 1)/2];

    returnUtilities.V0 = rayCast.V0(idx,:);
    returnUtilities.u = rayCast.u(idx,:);
    returnUtilities.v = rayCast.v(idx,:);
    returnUtilities.n = rayCast.n(idx,:);
    returnUtilities.uu = rayCast.uu(idx,:);
    returnUtilities.uv = rayCast.uv(idx,:);
    returnUtilities.vv = rayCast.vv(idx,:);
    returnUtilities.D = rayCast.D(idx,:);
    returnUtilities.infPoint = rayCast.infPoint;
    returnUtilities.tol = rayCast.tol;
end



function [in] = triIntersect(point, rayCastData)
    %triIntersect - Description
    %
    % Syntax: [in] = triIntersect(point, faceVerticies)
    %
    % Adapted from: http://geomalgorithms.com/a06-_intersect-2.html

    in = false;

    rayDir = rayCastData.infPoint - point;

    rayDir = repmat(rayDir, size(rayCastData.V0,1),1);

    w0 = point - rayCastData.V0;

    a = -dot(rayCastData.n, w0, 2);
    b = dot(rayCastData.n, rayDir, 2);

    r = a ./ b;

    test1Idx = r >= 0 & r <= 1;
    if ~any(test1Idx)
        return
    end

    I = point + (r(test1Idx, :) .* rayDir(test1Idx, :));

    w = I - rayCastData.V0(test1Idx, :);
    wu = dot(w, rayCastData.u(test1Idx, :), 2);
    wv = dot(w, rayCastData.v(test1Idx, :), 2);

    s = ((rayCastData.uv(test1Idx, :) .* wv) - (rayCastData.vv(test1Idx, :) .* wu))./ rayCastData.D(test1Idx, :) ;
    test2Idx = test1Idx(test1Idx) & (s >= 0 & s <= 1);
    test1Idx(test1Idx) = test2Idx;
    if ~any(test2Idx)
        return
    end

    t = ((rayCastData.uv(test1Idx, :) .* wu(test2Idx)) - (rayCastData.uu(test1Idx, :) .* wv(test2Idx))) ./ rayCastData.D(test1Idx, :) ;
    test3Idx = test2Idx(test2Idx) & (t >= 0 & (s(test2Idx) + t) <= 1);
    test2Idx(test2Idx) = test3Idx;
    test1Idx(test1Idx) = test3Idx;
    if ~any(test3Idx)
        return
    end

    intersections = sum(test3Idx);

    if rem(intersections, 2) ~= 0
        in = true;
    end
    if any(intersections, 'all') 
        suckEggsBoy = true;

    end
end




