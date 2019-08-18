function [output] = preProcess(input)
    % preProcess - This function is the master pre-processing file.
    %
    % Syntax: [output] = preProcess(input)
    %
    % Long description
    fname = strjoin([fpath  path_separator 'ds.dat'],'');

    datafile = fopen(fname);
    % Scans down the ds.dat file until the start of the element definitions.
    trashdata = 'a';
    START_NODES = '/com,*********** Nodes for'
    START_ELEMENTS = '/com,*********** Elements';

    while ~strncmpi(trashdata,startelements, 25)
        trashdata = fgetl(datafile);
    end
    % This loop runs down until it hits the line where the element type is
    % stored. The way the data is stored in the ds.dat file for some
    % elements means that it is necessary to skip a line when reading in
    % data, hence the skip line variable is a booean value. Tet elements
    % haven't been added here as yet.
    elid = 'a';
    numNodes = 0;
    while ~strncmpi(elid, 'eblock', 5)
        elid = fgetl(datafile);
        if strncmpi(elid, 'et', 2)
            [skipLine, numNodes, ~, type] = caseCheck(elid);
        end
    end

end
function [nNodes, nElements] = findSimulationData(resultFileObject)
    
    
end