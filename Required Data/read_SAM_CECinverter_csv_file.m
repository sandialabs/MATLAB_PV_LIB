% Read filename and write .mat file with inverter model
% parameters sets.

filename = 'CEC Inverters.csv';

a = importdata(filename);

txt = a.textdata;
dat = a.data;

hdrs = strsplit(txt{1,1},',');

name=txt(4:end,1);

emptyStruct = struct('name',[] ,  ...
    'Vac', [], 'Pac0', [], 'Pdc0', [], 'Vdc0', [], 'Ps0', [], ...
    'C0', [], 'C1', [], 'C2', [], 'C3', [], 'Pnt', [], ...
    'Vdcmax', [], 'Idcmax', [], 'MPPTLow', [], 'MPPTHi', [], ...
    'LibraryType', [], 'LibraryName', []);

NumberOfEntries = length(name);

InverterDatabase(1:NumberOfEntries) = emptyStruct;
InverterDatabase = InverterDatabase(:);

for i = 1:NumberOfEntries
    InverterDatabase(i).name = name{i};
    InverterDatabase(i).Vac = dat(i,1);
    InverterDatabase(i).Pac0 = dat(i,2);
    InverterDatabase(i).Pdc0 = dat(i,3);
    InverterDatabase(i).Vdc0 = dat(i,4);
    InverterDatabase(i).Ps0 = dat(i,5);
    InverterDatabase(i).C0 = dat(i,6);
    InverterDatabase(i).C1 = dat(i,7);
    InverterDatabase(i).C2 = dat(i,8);
    InverterDatabase(i).C3 = dat(i,9);
    InverterDatabase(i).Pnt = dat(i,10);
    InverterDatabase(i).Vdcmax = dat(i,11);
    InverterDatabase(i).Idcmax = dat(i,12);
    InverterDatabase(i).MPPTLow = dat(i,13);
    InverterDatabase(i).MPPTHi = dat(i,14);
    InverterDatabase(i).LibraryType = 'CEC Inverter';
    InverterDatabase(i).LibraryName = filename;
    
end

% Write the information from the CEC Library to the output variables.
CECInverterDB = InverterDatabase;
InverterNames = [SNLInverterDB.name]';

save 'CECInverterDatabaseSAM2015.6.30.mat' InverterNames CECInverterDB
