% Read filename and write .mat file with inverter model
% parameters sets.

clearvars

filename = 'CEC Modules.csv';

fmt = '%s%s%q%f%f%u%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%s';
a = readtable(filename,'Delimiter',',','Format',fmt,'HeaderLines',3,'ReadVariableNames',false);

emptyStruct = struct('name',[] , 't_noct', [] , 'a_c', [], 'Ns', [], 'i_sc_ref', [], ...
    'v_oc_ref', [], 'i_mp_ref', [], 'v_mp_ref', [], 'alpha_sc', [], 'beta_oc', [], ...
    'a_ref', [], 'IL_ref', [], 'I0_ref', [], 'Rs_ref', [], 'Rsh_ref', [], ...
    'adjust', [], 'gamma_r', [], 'source', [], 'LibraryType', [], 'LibraryName', []);


NumberOfEntries = size(a,1);

CECdatabase(1:NumberOfEntries) = emptyStruct;
CECdatabase = CECdatabase(:);

names = cellstr(a{:,1});
types = cellstr(a{:,22});

for i = 1:NumberOfEntries
    CECdatabase(i).name = names{i};
    CECdatabase(i).t_noct = a{i,4};
    CECdatabase(i).a_c = a{i,5};
    CECdatabase(i).Ns = a{i,6};
    CECdatabase(i).i_sc_ref = a{i,7};
    CECdatabase(i).v_oc_ref = a{i,8};
    CECdatabase(i).i_mp_ref = a{i,9};
    CECdatabase(i).v_mp_ref = a{i,10};
    CECdatabase(i).alpha_sc = a{i,11};
    CECdatabase(i).beta_oc = a{i,12};
    CECdatabase(i).a_ref = a{i,13};
    CECdatabase(i).IL_ref = a{i,14};
    CECdatabase(i).I0_ref = a{i,15};
    CECdatabase(i).Rs_ref = a{i,16};
    CECdatabase(i).Rsh_ref = a{i,17};
    CECdatabase(i).adjust = a{i,18};
    CECdatabase(i).gamma_r = a{i,19};
    CECdatabase(i).source = types{i};
    CECdatabase(i).LibraryType = 'CEC Module';
    CECdatabase(i).LibraryName = filename;
end

% Write the information from the CEC Library to the output variables.
CECModuleDB = CECdatabase;
ModuleNames = names(:);

save 'CECModuleDatabaseSAM2015.6.30.mat' ModuleNames CECModuleDB
