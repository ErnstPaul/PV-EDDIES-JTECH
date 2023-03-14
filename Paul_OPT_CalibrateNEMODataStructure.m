%% Paul_OPT_CalibrateNEMODataStructure.m (version 1.2)
%Author: Paul Ernst
%Date Created: 11/16/2021
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for release
%PE 2/9/2022 - Updated name
%PE 11/16/2021 - Created
%--------------------------------------
%Purpose: Calibrates a data structure for NEMO (or model of choice) on the local file system
%Inputs: Months, Years, input_dir for NEMO (or model of choice)
%Outputs: Large cell array for all NEMO (or model of choice) filenames
%--------------------------------------
function [input_final_NEM] = Paul_OPT_CalibrateNEMODataStructure(input_dir)
%% Months/Years/Indir
monthsC = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09';...
    '10'; '11'; '12'];
yearsC =  ['2016'; '2017'; '2018'; '2019'; '2020'; '2021'];

%% Construction
input_final = cell(31,length(monthsC),length(yearsC));
%Loop through all years/months/days to construct naming cell array to
%access later on when grabbing data from folders. Adjust above indices if
%you add any years to the dataset.
for i = 1:length(yearsC)
    input_dirs(1,i,:) =  [input_dir yearsC(i,:)];
    for j = 1:length(monthsC)
        input_dirs_months(i,j,:) = [transpose(squeeze(input_dirs(1,i,:))), '/', monthsC(j,:)];
        intermediate = dir(input_dirs_months(i,j,:));
        for k = 3:length(intermediate)
            input_final{k-2,j,i} = [transpose(squeeze(input_dirs_months(i,j,:))) '/' intermediate(k).name];
        end
    end
end

%% Final Array
input_final_NEM = input_final;