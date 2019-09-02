

CONSTANTS_FILE = "constants.mat";
delete(CONSTANTS_FILE)

ANSYS.files.ds = "ds.dat";
ANSYS.files.nodalSol = "nodalSolution.txt";

ANSYS.format.nodes = '%u32%f32%f32%f32';
ANSYS.format.elements = repmat('%u32',[1,19]);
ANSYS.format.nodalStress = '%u32%f32%f32%f32%f32%f32%f32';

ANSYS.regex.elementType = "(?:et,)\d*,(?<elementType>\d+)";
ANSYS.regex.block = "(?<fieldType>n|e)block";
ANSYS.regex.fields = "(?<fields>\d+)|(?<solid>solid)";
ANSYS.regex.elemEnd = "/wb,(?<sectionKey>\w+),(?<context>\w+)";
ANSYS.regex.nodeSolBlock = ["\s+(?<sectionKey>NODE)\s+","(?<stressFields>\<S\w{1,2})\s+"];


GENERAL.package.ANSYS = 1;
GENERAL.package.STRAND7 = 2;

GENERAL.dirs.outPath = "_output_data/";
GENERAL.dirs.prepPath = "_prep_data/";
GENERAL.dirs.prepPathN = "_prep_data/n/";
GENERAL.dirs.prepPathE = "_prep_data/e/";
GENERAL.dirs.prepPathS = "_prep_data/s/";
GENERAL.dirs.prepPathGeneral = "_prep_data/g/";

GENERAL.files.general = "globalSimData.mat";
GENERAL.files.stress = "s.mat";
GENERAL.files.nodes = "n.mat";
GENERAL.files.matExt = ".mat";

GENERAL.pathSep.OSX = "/";
GENERAL.pathSep.PC = "\";

GENERAL.varNames.elementIdx = "elementIdx";
GENERAL.varNames.connectivity = "connectivity";
GENERAL.varNames.elementData = "elementData";
GENERAL.varNames.nodeIdx = "nodeIdx";
GENERAL.varNames.coords = "coords";
GENERAL.varNames.maxNodes = "maxNodes";
GENERAL.varNames.nodalStressIdx = "nodalStressIdx";
GENERAL.varNames.stress = "stress";

GENERAL.path.globalData = GENERAL.dirs.prepPathGeneral + GENERAL.files.general;
GENERAL.path.nodes = GENERAL.dirs.prepPathGeneral + GENERAL.files.nodes;
GENERAL.path.stress = GENERAL.dirs.prepPathGeneral + GENERAL.files.stress;

vars = {"GENERAL", "ANSYS"};

save(CONSTANTS_FILE, "-v7.3", vars{:});

clear
