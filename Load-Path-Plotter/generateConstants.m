

CONSTANTS_FILE = "constants.mat"
delete(CONSTANTS_FILE)
% ANSYS.opt = {'MultipleDelimsAsOne',true};
ANSYS.files.ds = "ds.dat";
ANSYS.dirs.outPath = "_output_data/";
ANSYS.dirs.prepPath = "_prep_data/";
ANSYS.dirs.prepPathN = "_prep_data/n/";
ANSYS.dirs.prepPathE = "_prep_data/e/";
ANSYS.format.nodes = "%u32%f32%f32%f32";
ANSYS.format.elements = repmat("%u32",[1,19]);
ANSYS.regex.block = "(?<fieldType>n|e)block";
ANSYS.regex.fields = "(?<fields>\d+)|(?<solid>solid)";


save(CONSTANTS_FILE, "-v7.3", "ANSYS");

% matFileObject = matfile(outputPath, "Writable", true);