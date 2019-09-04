

CONSTANTS_FILE = "constants.mat";
delete(CONSTANTS_FILE)

ANSYS.files.ds = "ds.dat";
ANSYS.files.nodalSol = "nodalSolution.txt";

ANSYS.format.nodes = '%u32%f32%f32%f32';
ANSYS.format.elements = repmat('%u32',[1,19]);
ANSYS.format.nodalStress = '%u32%f32%f32%f32%f32%f32%f32';

ANSYS.regex.elementBlockData = "(?:et,)(?<iBody>\d+),(?<elementType>\d+)";
ANSYS.regex.block = "(?<fieldType>n|e)block";
ANSYS.regex.fields = "(?<fields>\d+)|(?<solid>solid)";
ANSYS.regex.elemEnd = "/wb,(?<sectionKey>\w+),(?<context>\w+)";
ANSYS.regex.nodeSolBlock = ["\s+(?<sectionKey>NODE)\s+","(?<stressFields>\<S\w{1,2})\s+"];


general.package.ANSYS = 1;
general.package.STRAND7 = 2;

general.dirs.outPath = "_output_data/"; % Path for outputs
general.dirs.prepPath = "_prep_data/"; % Path for all preprocessed data
general.dirs.prepPathN = "_prep_data/n/"; % Preprocessed node data sub directory
general.dirs.prepPathE = "_prep_data/e/"; % Preprocessed element data sub directory
general.dirs.prepPathS = "_prep_data/s/"; % Preprocessed stress data sub directory
general.dirs.prepPathU = "_prep_data/u/"; % Utilities sub directory
general.dirs.prepPathG = "_prep_data/g/";

general.files.general = "g.mat";
general.files.stress = "s.mat";
general.files.nodes = "n.mat";
general.files.utilities = "u.mat";
general.files.matExt = ".mat";

general.pathSep.OSX = "/";
general.pathSep.PC = "\";

general.varNames.elementIdx = "elementIdx";
general.varNames.connectivity = "connectivity";
general.varNames.elementData = "elementData";
general.varNames.nodeIdx = "nodeIdx";
general.varNames.coords = "coords";
general.varNames.maxNodes = "maxNodes";
general.varNames.nodalStressIdx = "nodalStressIdx";
general.varNames.stress = "stress";
general.varNames.general = "general";
general.varNames.V0 = "V0";
general.varNames.V1 = "V1";
general.varNames.u = "u";
general.varNames.v = "v";
general.varNames.n = "n";
general.varNames.uu = "uu";
general.varNames.uv = "uv";
general.varNames.vv = "vv";
general.varNames.D = "D";
general.varNames.infPoint = "infPoint";

general.constants.infPoint = [1e8 1e8 1e8];

general.path.general = general.dirs.prepPathG + general.files.general;
general.path.nodes = general.dirs.prepPathN + general.files.nodes;
general.path.stress = general.dirs.prepPathS + general.files.stress;
general.path.utilities = general.dirs.prepPathU + general.files.utilities;

vars = {"general", "ANSYS"};

save(CONSTANTS_FILE, "-v7.3", vars{:});

clear
