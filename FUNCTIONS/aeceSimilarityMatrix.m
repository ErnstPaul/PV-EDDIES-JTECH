%% SIMILARITY CALCULATION FUNCTION
function [simTable, simScore] = aeceSimilarityMatrix(labelAntiTest, labelCycloTest, labelAntiKey, labelCycloKey, NantiTest, NcycloTest, NantiKey, NcycloKey, ratioTest, ratioKey, origField)
%Similarity Matrix:
%AE Num Test      - CE Num Test       - AE Num Key        - CE Num Key
%AE True Positive - AE False Positive - AE False Negative - AE True Negative (NAN Masked out)
%CE True Positive - CE False Positive - CE False Negative - CE True Negative (NAN Masked out)
%Number Error     - Spatial Error      - Ratio Error      - Similarity Score

% Run through each type and add the NaN mask of the original land
nanMask = isnan(origField);
%viablePoints = nansum(nansum(~nanMask));
labelAntiTest(nanMask) = NaN; labelCycloTest(nanMask) = NaN;
labelAntiKey(nanMask) = NaN; labelCycloKey(nanMask) = NaN;

%Run through each AE matrix and get our categories
aeTRUEPOS = nansum(nansum((labelAntiTest==1) & (labelAntiKey==1)));
aeTRUENEG = nansum(nansum((labelAntiTest==0) & (labelAntiKey==0)));
aeFALSEPOS = nansum(nansum((labelAntiTest==1) & (labelAntiKey==0)));
aeFALSENEG = nansum(nansum((labelAntiTest==0) & (labelAntiKey==1)));
aeKEYTOTAL = aeTRUEPOS + aeFALSENEG;
aeKEYINCOR = aeTRUENEG + aeFALSEPOS;

%Run through each CE matrix and get our categories
ceTRUEPOS = nansum(nansum((labelCycloTest==1) & (labelCycloKey==1)));
ceTRUENEG = nansum(nansum((labelCycloTest==0) & (labelCycloKey==0)));
ceFALSEPOS = nansum(nansum((labelCycloTest==1) & (labelCycloKey==0)));
ceFALSENEG = nansum(nansum((labelCycloTest==0) & (labelCycloKey==1)));
ceKEYTOTAL = ceTRUEPOS + ceFALSENEG;
ceKEYINCOR = aeTRUENEG + ceFALSEPOS;

%Calculate agg. score: number of eddies (error = [test-key]/key
errorAnti = abs((NantiTest - NantiKey) / (NantiKey));
errorCyclo = abs((NcycloTest - NcycloKey) / (NcycloKey));
errorNumTotal = mean([errorAnti, errorCyclo]);

%Calculate agg. score: area of eddies covered
%First: % error in correctly identified as positive
errorSpatialCorrectTotal = abs(((aeTRUEPOS + ceTRUEPOS) - (aeKEYTOTAL + ceKEYTOTAL)) / (aeKEYTOTAL + ceKEYTOTAL));
%Second: % error in corrected identified as negative
errorSpatialIncorrectTotal = abs(((aeTRUENEG + ceTRUENEG) - (aeKEYINCOR + ceKEYINCOR)) / (aeKEYINCOR + ceKEYINCOR));

%Calculate agg. score: ratios of eddies (shape)
errorRatioTotal = abs((ratioTest - ratioKey)/ratioKey);

%Compile total score out of 100
simScore = (1-mean([errorNumTotal, errorSpatialCorrectTotal, errorSpatialIncorrectTotal, errorRatioTotal]))*100;

%Compile table
simTable = [NantiTest, NcycloTest, NantiKey, NcycloKey; ...
    aeTRUEPOS, aeFALSEPOS, aeFALSENEG, aeTRUENEG; ...
    ceTRUEPOS, ceFALSEPOS, ceFALSENEG, ceTRUENEG; ...
    errorNumTotal, errorSpatialCorrectTotal, errorSpatialIncorrectTotal, errorRatioTotal];
%         FieldsSimilarity = [
%             '(1,1) AE Number Test            '; ...
%             '(1,2) CE Number Test            '; ...
%             '(1,3) AE Number Key             '; ...
%             '(1,4) CE Number Key             '; ...
%             '(2,1) AE True Positive          '; ...
%             '(2,2) AE False Positive         '; ...
%             '(2,3) AE False Negative         '; ...
%             '(2,4) AE True Negative          '; ...
%             '(3,1) CE True Positive          '; ...
%             '(3,2) CE False Positive         '; ...
%             '(3,3) CE False Negative         '; ...
%             '(3,4) CE True Negative          '; ...
%             '(4,1) Number Error              '; ...
%             '(4,2) Spatial Error Correct     '; ...
%             '(4,3) Spatial Error Incorrect   '; ...
%             '(4,4) Ratio Error               '];
end