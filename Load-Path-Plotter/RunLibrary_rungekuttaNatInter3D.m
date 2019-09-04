function [x_path, y_path,z_path, intensity] =  RunLibrary_rungekuttaNatInter3D(...
    general, xseed,yseed,zseed, ReversePath, wb)

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

    [in, Element] = point_in_element(p0, utilities);

    if in
        %Get and set stress function
        [F, Fs1, Fs2] = setInterpFunc(Element, general.constants.pathDir, ReversePath);
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
    element_change = false;

    while w <= general.constants.pathLength && in ~= false
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
        [in, new_Element] = point_in_element(p0, utilities);

        %If the point is outside the local radius, we attempt to find it
        %globablly. The path is projected along its last vector in an
        %attempt to get it to 'land' in another element for the case where
        %its in a small gap between elements.

        if in && new_Element(1).ElementNo ~= Element(1).ElementNo
            [F, Fs1, Fs2] = setInterpFunc(Element,general.constants.pathDir, ReversePath);
        end

        Element = new_Element;
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
    if ReversePath
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
end
function [F, Fs1, Fs2] = setInterpFunc(Element, pathDir, ReversePath)
    %Natural interpolation method is used to form a stress function to then
    %compute the paths.
    % TODO: - How many nodes do i need?
    surr_elements = RunLibrary_surrounding_elemnts(Element, Element);
    nodes = unique([surr_elements(:).nodes]);
    coordx = [nodes(:).xCoordinate]';
    coordy = [nodes(:).yCoordinate]';
    coordz = [nodes(:).zCoordinate]';

    switch general.constants.pathDir
        case 'X'
            stress_tensor = [[nodes(:).xStress]', [nodes(:).xyStress]', [nodes(:).xzStress]'];
        case 'Y'
            stress_tensor = [[nodes(:).yStress]', [nodes(:).xyStress]', [nodes(:).yzStress]'];
        case 'Z'
            stress_tensor = [[nodes(:).zStress]', [nodes(:).xzStress]', [nodes(:).yzStress]'];
    end
    F = scatteredInterpolant(coordx, coordy, coordz, stress_tensor(:,1), 'natural');
    Fs1 = scatteredInterpolant(coordx, coordy, coordz, stress_tensor(:,2), 'natural');
    Fs2 = scatteredInterpolant(coordx, coordy, coordz, stress_tensor(:,3), 'natural');
    %Flips the stress function when the backwards path is being computed.
    if ReversePath
        F=@(x,y, z) -F(x,y, z);
        Fs1=@(x,y,z) -Fs1(x,y,z);
        Fs2=@(x,y,z) -Fs2(x,y,z);
    end
end

function [varargout] = point_in_element(p0, utilities)
    in_test = triIntersect(p0, utilities)
    in = in_test;
    Element = PartArr(1).elements(in_test);
    varargout = {in, Element};
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
            curDir = general.dirs.workingDir + general.dirs.prepPathE...
                     + body + general.files.matExt;
        case "n"
            curDir = general.dirs.workingDir + general.path.nodes;
        case "s"
            curDir = general.dirs.workingDir + general.path.stress;
        case "g"
            curDir = general.dirs.workingDir + general.path.general;
        case "u"
            curDir = general.dirs.workingDir + general.path.utilities;
    end
    if nargin < 3
        loadedVar = load(curDir);
    else
        loadedVar = load(curDir, varToLoad);
        loadedVar = loadedVar.(varToLoad);
    end
end
function [in] = triIntersect(point, utilities)
    %triIntersect - Description
    %
    % Syntax: [in] = triIntersect(point, faceVerticies)
    %
    % Adapted from: http://geomalgorithms.com/a06-_intersect-2.html

    point = point';

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

function [varagout] = getInterpolationData(point, general)
    %getInterpolationData - Description
    %
    % Syntax: [varagout] = getInterpolationData(point, general)

    
end