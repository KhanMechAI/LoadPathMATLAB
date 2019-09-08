function [x_path, y_path,z_path, intensity] =  RunLibrary_rungekuttaNatInter3D(...
    general, xseed,yseed,zseed, reversePath, wb)

    %p0 is initial seed point. Projection multiplier is used to 'jump' gaps
    %between parts in the model. It will project the path from onee part to
    %another. This procedure is a big source of time, some thought needs to
    %be given to how to optimise this routine.

    utilities = getData(general, "u");
    % [general] = getData(general, dataClass, varargin)

    projectionMultiplier = 2;
    p0 = [xseed; yseed; zseed];

    %Anonymous functions to streamline the organisation of the path
    %direction.
    Vx =@(stress, shearxy,shearxz) [stress; shearxy; shearxz];
    Vy =@(stress, shearxy, shearyz) [shearxy; stress; shearyz];
    Vz =@(stress, shearxz, shearyz) [shearxz; shearyz; stress];

    switch lower(general.constants.pathDir)
        case 'x'
            V = Vx;
        case 'y'
            V = Vy;
        case 'z'
            V = Vz;
    end
    %Locate the seed point in the model globally.
    [currBody] = findCurrentBody(p0);
    if ~currBody
        fprintf('Seed Point (%f, %f, %f) not in solution domain\n', xseed, yseed, zseed);
        currBody=1;
        figure;
        plotBody
        plotPoint
        plotRay
        [currBody] = findCurrentBody(p0);
        x_path = [];
        y_path = [];
        z_path =[];
        intensity = [];
        return
    end

    [element] = findCurrentElement(p0, currBody);

    if element
        %Get and set stress function
        [F, Fs1, Fs2] = getInterpolationData(element, utilities.bodyData{currBody}.elementInterpFunc, general.constants.pathDir, reversePath);
    else
        fprintf('Seed Point (%f, %f, %f) not in solution domain\n', xseed, yseed, zseed);
        x_path = [];
        y_path = [];
        z_path =[];
        intensity = [];
        return
    end

    %Populating with NaN's prevents plotting errors later.
    p = NaN(3,general.constants.pathLength,'double');
    intensity = NaN(1,general.constants.pathLength,'double');
    w = 1;
    prevElement = element;
    currElementUtilities = getElementUtilities(element, utilities.locatingData{currBody});

    while element ~= false && w <= general.constants.pathLength
        %Terminate program if cancel button is pressed
        if getappdata(wb,'canceling')
            delete(wb)
            break
        end
        %Interpolate stress initially, and get the relative normalisation
        %value. If the element is unchanged, the same stress function can
        %be used.

        p(:,w) = p0;
        stress = F(p0(1), p0(2), p0(3));
        shear1 = Fs1(p0(1), p0(2), p0(3));
        shear2 = Fs2(p0(1), p0(2), p0(3));
        intensity(w) = norm([stress shear1 shear2]);
        %Find dp1 at first interpolation.
        %Normalise the poining vector relative to the initial test point.
        %Calculate new points:

        %Runge-Kutta
        dp1 = V(stress, shear1, shear2)*general.constants.stepSize/intensity(w);
        p1 = p0 + dp1;

        dp2 = stress_interp(p1);
        p2 = p0 + 0.5*dp2;

        dp3 = stress_interp(p2);
        p3 = p0 + 0.5*dp3;

        dp4 = stress_interp(p3);

        p0 = p0 + 1/6 * (dp1 + 2*dp2 + 2*dp3 + dp4);

        %Locate element that the point is inside
        if ~triIntersect(p0, currElementUtilities)
            [currBody] = findCurrentBody(p0, currBody);
            if ~currBody
                element = false;
                disp('not in body')
            else
                [element] = findCurrentElement(p0, currBody);
                if element
                    currElementUtilities = getElementUtilities(element, utilities.locatingData{currBody});
                else
                    disp('in body, but not in element')
                end
            end
        end

        if element && prevElement ~= element
            [F, Fs1, Fs2] = getInterpolationData(element, utilities.bodyData{currBody}.elementInterpFunc, general.constants.pathDir, reversePath);
        end
        w=w+1;
    end
    nancols = ~isnan(p(1,:));
    if nancols > 1
        %To keep plot inside domain
        nancols = nancols-1;
    end
    x_path = p(1,:);
    y_path = p(2,:);
    z_path = p(3,:);
    intensity = intensity(:)';
    
    if reversePath
        x_path = fliplr(x_path);
        y_path = fliplr(y_path);
        z_path = fliplr(z_path);
        intensity = fliplr(intensity);
    end

    function [d_point] = stress_interp(p)
        stress = F(p(1), p(2), p(3));
        shear1 = Fs1(p(1), p(2), p(3));
        shear2 = Fs2(p(1), p(2), p(3));
        d_point =  V(stress, shear1, shear2)*general.constants.stepSize/intensity(w);
    end

    function plotBody(varargin)
        if nargin == 0 
            for k = 1:length(utilities.locatingData)
                bodyHullIdx = utilities.locatingData{k}.bodyHullIdx;
                bodyFaceCoords = utilities.locatingData{k}.faceCoords(:,bodyHullIdx,:);
                plotFaces(bodyFaceCoords, 'black', 0.2);
            end
        else
            bodyHullIdx = utilities.locatingData{varargin}.bodyHullIdx;
            bodyFaceCoords = utilities.locatingData{varargin}.faceCoords(:,bodyHullIdx,:);
            plotFaces(bodyFaceCoords, 'black', 0.2);
        end
    end
    
    function plotElement()
    
        indexingArray = utilities.locatingData{currBody}.elementFaceIndexingArray;
        mappingVector = utilities.locatingData{currBody}.faceMappingVector;
        elidx = mappingVector(indexingArray(element,:));
    
        elementFaceCoords = utilities.locatingData{currBody}.faceCoords(:,elidx,:);
        plotFaces(elementFaceCoords, 'red', 1);
    end

    function plotPoint(varargin)
        if nargin > 0
            point = varargin{1};
        else
            point = p0;
        end
        gcf;
        hold on
        plot3(point(1), point(2), point(3), '*', 'MarkerFaceColor', 'red')
        hold off
    end

    function plotRay(varargin)
        point = zeros([2,3]);
        point(1,:) = [1e8 1e8 1e8]
        if nargin > 0
            point(2,:) = varargin{1};
        else
            point(2,:) = p0;
        end
        gcf;
        hold on
        plot3(point(:,1), point(:,2), point(:,3), '*', 'MarkerFaceColor', 'red')
        hold off
    end

    function [element] = findElementInBody(p0, currBody)
        %findElementInBody - Description
        %
        % Syntax: [element] = findElementInBody(p0, currBody)
        k = 1;
        inElement = false;
        while k < utilities.bodyData.elementBlockData(1) && ~inElement
            elementUtilities = getElementUtilities(testElements(k), utilities.locatingData{currBody});
            inElement = triIntersect(p0, elementUtilities);
            if ~inElement
                k = k +1;
            end
        end
        if inElement
            element = testElements(k);
        else
            element = false;
        end
    end

    function [element] = findCurrentElement(p0, currBody)
        try
            testElements = knnsearch(utilities.locatingData{currBody}.kdTreeElements,p0','K',8);
        catch
            aa=1;
        end
        k = 1;
        inElement = false;
        while k < size(testElements, 2) && ~inElement
            elementUtilities = getElementUtilities(testElements(k), utilities.locatingData{currBody});
            inElement = triIntersect(p0, elementUtilities);
            if ~inElement
                k = k +1;
            end
        end
        if inElement
            element = testElements(k);
        else
            element = false;
        end
    end

    function [body] = findCurrentBody(p0, varargin)
        if nargin > 1
            body = varargin{1};
        else
            body = 1;
        end
        [bodyUtilities] = getBodyHullUtilities(utilities.locatingData{body});
        inBody = triIntersect(p0, bodyUtilities);
        if inBody
            return
        else
            k = 1;
            while k < size(utilities.locatingData{:}, 2) && ~inBody
                if k ~= body
                    [bodyUtilities] = getBodyHullUtilities(elementNo, utilities.locatingData{k});
                    inBody = triIntersect(p0, bodyUtilities);
                end
                if ~inBody
                    k = k+1;
                end
            end
        end
        if inBody
            body = k;
        else
            body = false;
        end
    end
end

function [varargout] = projection(in, p0, Element)
    if ~in
        extension = 1;
        while ~in && extension < projectionMultiplier+1
            R = (p0 - p(:,w)) * extension * 2 + p0;
            [in, new_Element] = point_in_element(R, PartArr);
            extension = extension+1;
        end
        if in
            p0 = R;
            Element = new_Element;
        end
    end
    varargout = {in, p0, Element};
end
function [loadedVar] = getData(general, dataClass, varargin)
    if length(varargin) > 0
        varToLoad = varargin{1,1};
        if nargin > 3
            body = string(varargin{1,2});
        end
    end
    switch dataClass
        case "e"
            curDir = general.local.dirs.workingDir + general.dirs.prepPathE...
                     + body + general.files.matExt;
        case "n"
            curDir = general.local.dirs.workingDir + general.path.nodes;
        case "s"
            curDir = general.local.dirs.workingDir + general.path.stress;
        case "g"
            curDir = general.local.dirs.workingDir + general.path.general;
        case "u"
            curDir = general.local.dirs.workingDir + general.path.utilities;
    end
    if nargin < 3
        loadedVar = load(curDir);
    else
        loadedVar = load(curDir, varToLoad);
        loadedVar = loadedVar.(varToLoad);
    end
end
function [in] = triIntersect(point, rayCastData)
    %triIntersect - Description
    %
    % Syntax: [in] = triIntersect(point, faceVerticies)
    %
    % Adapted from: http://geomalgorithms.com/a06-_intersect-2.html

    point = point';

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

    % test1 = s >= 0 & s <= 1;
    % test2 = t >= 0 & (s + t) <= 1;
    % test3 = r >= 0;
    intersections = sum(test3Idx);

    if rem(intersections, 2) ~= 0
        in = true;
    end
    if any(intersections, 'all') 
        suckEggsBoy = true;

    end
end

function [F, Fs1, Fs2] = getInterpolationData(element, interpFunctions, pathDir, reversePath)
    %getInterpolationData - Description
    %
    % Syntax: [varagout] = getInterpolationData(point, general)
    switch pathDir
        case 'X'
            interpolantIdx = [1 4 6]; %Fxx Fxy Fxz
        case 'Y'
            interpolantIdx = [2 4 5]; %Fyy Fxy Fyz
        case 'Z'
            interpolantIdx = [3 6 5]; %Fzz Fyz Fxz
    end
    [F, Fs1, Fs2] = interpFunctions{element}{interpolantIdx};
    if reversePath
        F=@(x,y, z) -F(x,y, z);
        Fs1=@(x,y,z) -Fs1(x,y,z);
        Fs2=@(x,y,z) -Fs2(x,y,z);
    end
end

function [elementUtilities] = getElementUtilities(elementNo, utilities)
    %getElementUtilities - Description
    %
    % Syntax: [elementUtilities] = getElementUtilities(elementNo, utilities)
    indexingArray = utilities.elementFaceIndexingArray;
    mappingVector = utilities.faceMappingVector;
    elidx = indexingArray(elementNo,:);
    idx = [mappingVector(elidx)];
    idx = [idx; idx + size(utilities.rayCastVars.V0, 1)/2];
    elementUtilities.V0 = utilities.rayCastVars.V0(idx,:);
    elementUtilities.u = utilities.rayCastVars.u(idx,:);
    elementUtilities.v = utilities.rayCastVars.v(idx,:);
    elementUtilities.n = utilities.rayCastVars.n(idx,:);
    elementUtilities.uu = utilities.rayCastVars.uu(idx,:);
    elementUtilities.uv = utilities.rayCastVars.uv(idx,:);
    elementUtilities.vv = utilities.rayCastVars.vv(idx,:);
    elementUtilities.D = utilities.rayCastVars.D(idx,:);
    elementUtilities.infPoint = utilities.rayCastVars.infPoint;
    elementUtilities.tol = utilities.rayCastVars.tol;
end

function [bodyUtilities] = getBodyHullUtilities(utilities)
    %getElementUtilities - Description
    %
    % Syntax: [bodyUtilities] = getElementUtilities(elementNo, utilities)
    mappingVector = utilities.faceMappingVector;
    hullIdx = utilities.bodyHullIdx;
    idx = [mappingVector(hullIdx)];
    idx = [idx; idx + size(utilities.rayCastVars.V0, 1)/2];
    bodyUtilities.V0 = utilities.rayCastVars.V0(idx,:);
    bodyUtilities.u = utilities.rayCastVars.u(idx,:);
    bodyUtilities.v = utilities.rayCastVars.v(idx,:);
    bodyUtilities.n = utilities.rayCastVars.n(idx,:);
    bodyUtilities.uu = utilities.rayCastVars.uu(idx,:);
    bodyUtilities.uv = utilities.rayCastVars.uv(idx,:);
    bodyUtilities.vv = utilities.rayCastVars.vv(idx,:);
    bodyUtilities.D = utilities.rayCastVars.D(idx,:);
    bodyUtilities.infPoint = utilities.rayCastVars.infPoint;
    bodyUtilities.tol = utilities.rayCastVars.tol;
end



function plotFaces(faceArray, colour, alp)
    %plotFaces - Description
    %
    % Syntax: plotFaces(faceArray)

    % triangulatedFaceCoords = [faceArray([1:3], :,:), faceArray([1,3,4], :,:)];
    
    fig = gcf;
    hold on
    patch('XData',faceArray(:, :, 1),'YData',faceArray(:, :, 2),'ZData',faceArray(:, :, 3), 'EdgeColor',colour,'FaceColor','none', 'EdgeAlpha', alp)
    hold off
end