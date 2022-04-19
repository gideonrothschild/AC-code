


runscript=0;



    Metadatafilename='C:\Users\gid\Dropbox (University of Michigan)\from_box\LabStartup\Code\ACSpeedCodingAA\IMAGING\MetaDataCVGR2_n7.xlsx';
    resultsFile='Z:\DataBackup\gideon\ProcessedData\allResults040520';
    


if runscript
    
    [d v datatxt]=xlsread([Metadatafilename]);
    experiments_tmp=datatxt(:,1);
    stim_tmp=datatxt(:,2);
    animalState_tmp=datatxt(:,3);
    expBatch_tmp=datatxt(:,4);
    videofile=datatxt(:,6);
    locomotionMethod=datatxt(:,10);
    % removing rows that don't have data
    datarows=[];
    for i=1:size(experiments_tmp,1)
        if ~isempty(strfind(experiments_tmp{i},'TSeries'))
            datarows=[datarows i];
        end
        
    end
    
    experiments=experiments_tmp(datarows);
    stim=stim_tmp(datarows);
    animalState=animalState_tmp(datarows);
    expBatch=expBatch_tmp(datarows);
    videofile=videofile(datarows);
    locomotionMethod=locomotionMethod(datarows);
    curExpBatch='';
    experimentsInBatch={};
    stimuliInBatch={};
    videoFilesInBatch={};
    animalStateInBatch={};
    allresults={};
    locomotionMethodInBatch={};
    
    
    for i=1:length(experiments)
        
        % moving to new batch, so processing the previous
        if ~strcmp(expBatch{i},curExpBatch)&i>1
            
            % for when it is just one number, Matlab interprets
            % it as double
            if isnumeric(curExpBatch), curExpBatch=num2str(curExpBatch);end
            if isnan(videoFilesInBatch{1})
             
                perExpResults=analyzeBatchedDataCVNoVideo(curExpBatch,experimentsInBatch,stimuliInBatch,animalStateInBatch,videoFilesInBatch,locomotionMethodInBatch);
                
            else
                
                perExpResults=analyzeBatchedDataCV(curExpBatch,experimentsInBatch,stimuliInBatch,animalStateInBatch,videoFilesInBatch,locomotionMethodInBatch);
            end
            allresults{end+1}=perExpResults;
            save(resultsFile,'allresults')
            
            close all;
            % moving to next batch
            
            curExpBatch=expBatch{i};
            experimentsInBatch={};
            stimuliInBatch={};
            animalStateInBatch={};
            videoFilesInBatch={};
            locomotionMethodInBatch={};
            
            videoFilesInBatch{end+1}=videofile{i};
            experimentsInBatch{end+1}=experiments{i};
            stimuliInBatch{end+1}=stim{i};
            animalStateInBatch{end+1}=animalState{i};
            locomotionMethodInBatch{end+1}=locomotionMethod{i};
            
        else
            experimentsInBatch{end+1}=experiments{i};
            curExpBatch=expBatch{i};
            stimuliInBatch{end+1}=stim{i};
            animalStateInBatch{end+1}=animalState{i};
            videoFilesInBatch{end+1}=videofile{i};
            locomotionMethodInBatch{end+1}=locomotionMethod{i};
            
        end
        
        
    end
    return
else
    load(resultsFile)
end

%%


rng(77)
animalstatetorun='AWK';
plotcellsresponsiveness=0;
minCellsSigCorrsPerSession=4;
toPrintIndividualExps=0;
toCorrectSignificance=1;

cellsresponsivenessFigDir='E:\GideonData\IMAGING\responsiveness\';
cellsresponsivenessbyspeedFigDir='E:\GideonData\IMAGING\responsivenessbyspeed\';
conttracesfolder='E:\GideonData\IMAGING\conttraces\';
NATimagesfolder1='E:\GideonData\IMAGING\NAT1new\';
NATimagesfolder2='E:\GideonData\IMAGING\NAT2\';

spatialfolder='E:\GideonData\IMAGING\Spatial\';
spatialfolder2='E:\GideonData\IMAGING\Spatial2\';


allResponsiveSigCorrs=[];
allResponsiveSigCorrsPerSession=[];
allResponsiveNoiseCorrs=[];

allMeanResponses=[];
allMeanResponsesInclNonSig=[];


alllocomotionEnhancedPairsXCorrs=[];
alllocomotionEnhancedPairsXCorrsImmob=[];
alllocomotionEnhancedPairsXCorrsMove=[];


alllocomotionEnhancedPairsXCorrsShuf=[];
alllocomotionEnhancedPairsXCorrsImmobShuf=[];
alllocomotionEnhancedPairsXCorrsMoveShuf=[];

alllocomotionSuppressedPairsXCorrs=[];
alllocomotionSuppressedPairsXCorrsImmob=[];
alllocomotionSuppressedPairsXCorrsMove=[];

alllocomotionSuppressedPairsXCorrsShuf=[];
alllocomotionSuppressedPairsXCorrsImmobShuf=[];
alllocomotionSuppressedPairsXCorrsMoveShuf=[];

alllocomotionUnmodPairsXCorrs=[];

allacrossClassesXCorrs=[];
allacrossClassesXCorrsMove=[];
allacrossClassesXCorrsImmob=[];


allacrossClassesXCorrsShuf=[];
allacrossClassesXCorrsMoveShuf=[];
allacrossClassesXCorrsImmobShuf=[];


alllocomotionEnhancedNoiseCorrs=[];
alllocomotionEnhancedNoiseCorrsMove=[];
alllocomotionEnhancedNoiseCorrsImmob=[];

alllocomotionSuppressedNoiseCorrs=[];
alllocomotionSuppressedNoiseCorrsMove=[];
alllocomotionSuppressedNoiseCorrsImmob=[];

allacrossClassesNoiseCorrs=[];
allacrossClassesNoiseCorrsMove=[];
allacrossClassesNoiseCorrsImmob=[];


alllocomotionEnhancedNoiseCorrsShuf=[];
alllocomotionEnhancedNoiseCorrsMoveShuf=[];
alllocomotionEnhancedNoiseCorrsImmobShuf=[];

alllocomotionSuppressedNoiseCorrsShuf=[];
alllocomotionSuppressedNoiseCorrsMoveShuf=[];
alllocomotionSuppressedNoiseCorrsImmobShuf=[];

allacrossClassesNoiseCorrsShuf=[];
allacrossClassesNoiseCorrsMoveShuf=[];
allacrossClassesNoiseCorrsImmobShuf=[];


alllocomotionUnmodNoiseCorrs=[];

allPerStimSignificance=[];
allOverallSignificance=[];
allDiscrError=[];
allNumCells=[];
allNumSigCells=[];
allMeanResp=[];
allstdflucts=[];
allstdfluctsSessionMean=[];
allDiscrErrMatchedToSigCorr=[];
sValsAllExpsSpont=[];
sValsAllExpsSpont_onlyResponsive=[];
sValsAllExpsEvoked_onlyResponsive=[];
sValsAllExpsEvoked_SpontonlyResponsive=[];
sValsAllExpsEvoked=[];
sValsAllExpsEvoked_Spont=[];
allmeanspontLocomotion=[];
allmeanspontImmobility=[];
allmeanevokedLocomotion=[];
allmeanevokedImmobility=[];
allmeanEvokedMinusPreLocomotion=[];
allmeanEvokedMinusPreImmobility=[];
participatingEnsemblesDecoding=[];
participatingEnsemblesSpeedPrediction=[];
locomotionEnhancedNoiseCorrs=[];
locomotionSuppressedNoiseCorrs=[];
acrossClassesNoiseCorrs=[];

allmeanspontLocomotion_responsive=[];
allmeanspontImmobility_responsive=[];
allmeanevokedLocomotion_responsive=[];
allmeanevokedImmobility_responsive=[];
allmeanEvokedMinusPreLocomotion_responsive=[];
allmeanEvokedMinusPreImmobility_responsive=[];

allmeanspontLocomotion_responsiveIMMOB=[];
allmeanspontImmobility_responsiveIMMOB=[];
allmeanevokedLocomotion_responsiveIMMOB=[];
allmeanevokedImmobility_responsiveIMMOB=[];
allmeanEvokedMinusPreLocomotion_responsiveIMMOB=[];
allmeanEvokedMinusPreImmobility_responsiveIMMOB=[];

allCellResponsiveness=[];
allCellImmobileResponsiveness=[];
allCellImmobileResponsivenessCONST=[];
difEvokedMinusPreLocomotionVsImmobilityTON3=[];
allmeanEvokedMinusPreLocomotionTON3=[];
allmeanEvokedMinusPreImmobilityTON3=[];

difEvokedMinusPreLocomotionVsImmobilityNAT3=[];
allmeanEvokedMinusPreLocomotionNAT3=[];
allmeanEvokedMinusPreImmobilityNAT3=[];

difEvokedMinusPreLocomotionVsImmobilityNAT=[];
allmeanEvokedMinusPreLocomotionNAT=[];
allmeanEvokedMinusPreImmobilityNAT=[];

acrossStimResponsiveness=[];
sValsENatCell_stim=[];
allmeanevokedLocomotion_natcellstim=[];
allmeanevokedImmobility_natcellstim=[];
allCellResponsiveness_anystim=[];
allCellResponsiveness_perstim=[];
allcell_stimdFFmodulation_nat=[];
allcell_stimSNRmodulation_nat=[];
allcell_acrossStimSNRmodulation_nat=[];
cell_specificmeanevokedLocomotion_natcellstim =[];
cell_specificmeanevokedImmobility_natcellstim =[];
cell_specificmeanSNRLocomotion_natcellstim =[];
cell_specificmeanSNRImmobility_natcellstim =[];
ensemble_ID = [];
ensembleCELL_ID={};
allrelativeDistNeg=[];
allrelativeDistPos=[];
allrelativeDistNeg_cricket=[];
allrelativeDistPos_cricket=[];
allrelativeDistNeg_sparrow=[];
allrelativeDistPos_sparrow=[];
allrelativeDistNeg_scratch=[];
allrelativeDistPos_scratch=[];
allrelativeDistNeg_water=[];
allrelativeDistPos_water=[];

alldffspeedcorrs=[];
alldffspeedcorrsSignificance=[];
alldffspeedcorrsNo0=[];
alldffspeedcorrsSignificanceNo0=[];

allalldff=[];
AllImmobileStimDecodingError=[];
AllLocomotionStimDecodingError=[];
AllLocomotionStimDecodingErrorSHUF=[];
AllLocomotionStimDecodingPVAL=[];

AllPreWinStateDecodingError=[];
AllStimWinStateDecodingError=[];
AllStimWinStateDecodingErrorSHUF=[];
AllImmobileStimDecodingErrorSHUF=[];
AllImmobileStimDecodingPVAL=[];


AllPreWinStateDecodingErrorSHUF=[];
AllPreWinStateDecodingPVAL=[];
AllStimWinStateDecodingErrorSHUF=[];
AllStimWinStateDecodingPVAL=[];


AllRealVsShufStateDecoding=[];
AllRealVsShufStimDecoding=[];



AllStimAndStateDecodingError=[];
AllStimAndStateDecodingErrorShuf=[];

AllStimAndStateDecodingPVAL=[];

allCellIndsForDecoding={};
difSpontLocomotionVsImmobilityAllExps=[];
difSpontLocomotionVsImmobilityAllExps_onlyResponsive=[];
difSpontLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB=[];
difEvokedLocomotionVsImmobilityAllExps=[];
difEvokedLocomotionVsImmobilityAllExps_onlyResponsive=[];

perEnsembleEvokedMinusPreProportions=[];


difEvokedMinusPreLocomotionVsImmobilityAllExps=[];
difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsive=[];
difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB=[];
b2=fir1(21,0.2);
allProportionOfPosSpeedModCellsInEnsemble=[];

ton3folder='E:\GideonData\IMAGING\TON3\'

responsiveinds={};
immobresponsiveinds={};

allinds2={};
stimDetectionErrorImmobSingleCell=[];
stimDetectionErrorMoveSingleCell=[];
stimDetectionPVALImmobSingleCell=[];
stimDetectionPVALMoveSingleCell=[];
stateDetectionErrorSingleCell=[];
stateDetectionPVALSingleCell=[];
stateAndStimDetectionErrorSingleCell=[];
stateAndStimDetectionPVALSingleCell=[];

acrossExperimentsSpeeds=[];
alldffspeedcorrsCONST=[];
alldffspeedcorrsNONE=[];
alldffspeedcorrsRanges=[];
allmeanDuringImmobileNONE=[];
allmeanDuringMoveNONE=[];

allmeanDuringImmobileCONST=[];
allmeanDuringMoveCONST=[];

allSpeedPredictionRs=[];
allSpeedPredictionRsNo0=[];
allSpeedPredictionRsSHUF=[];
allSpeedPredictionRsNo0SHUF=[];
predictionEnsembleSize=[];
decodingCellResponsiveness=[];
allminnumtrials=[];
allnumcellsperensembledecoding=[];
preImmobVar=[];
preMoveVar=[];


allMovePSTHsResponsiveIMMOB=nan(1,15);
allImmobPSTHsResponsiveIMMOB=nan(1,15);
maxLocomotionSuppressedPairsXCorrs=[];
maxLocomotionEnhancedPairsXCorrs=[];
maxWithinClassesXCorrs=[];
maxAcrossClassesXCorrs=[];



for i=1:length(allresults)
    [num2str(i) 'out of ' num2str(length(allresults))]
    numepochs=length(allresults{i});
    NATresps=[];
    BBNresps=[];
    NATrespsLocomotionByStim=[];
    BBNrespsLocomotionByStim=[];
    alldff=[];
    NATrespsMOV=[];
    NATrespsImmobile=[];
    BBNrespsMOV=[];
    BBNrespsImmobile=[];
    allmov=[];
    
    NAT3respsMOV=[];
    NAT3respsImmobile=[];
    NAT3resps=[];
    NAT3respsLocomotionByStim=[];
    alldff3=[];
    allmov3=[];
    
    
    TON3respsMOV=[];
    TON3respsImmobile=[];
    TON3resps=[];
    TON3respsLocomotionByStim=[];
    alldffTON3=[];
    allmovTON3=[];
    
    MIX3respsMOV=[];
    MIX3respsImmobile=[];
    MIX3resps=[];
    MIX3respsLocomotionByStim=[];
    alldffMIX3=[];
    allmovMIX3=[];
    
    allallSpeecCorrVals=[];
    acrossSessionsSpeeds=[];
    acrossSessionsSpeedsNONE=[];
    acrossSessionsSpeedsCONST=[];
    epochdffspeedcorrs=[];
    epochdffspeedcorrsSignificance=[];
    epochdffspeedcorrsNo0=[];
    epochdffspeedcorrsSignificanceNo0=[];
    epochdffspeedcorrsNONE=[];
    epochdffspeedcorrsCONST=[];
    alldffNONE=[];
    allmovNONE=[];
    alldffCONST=[];
    allmovCONST=[];
    
    for j=1:numepochs
        if strcmp(allresults{1,i}{1,j}.state,animalstatetorun)
            if strcmp(allresults{1,i}{1,j}.stim,'MIX3')
                MIX3respsMOV=cat(3,MIX3respsMOV, allresults{1,i}{1,j}.expRespsMove);
                MIX3respsImmobile=cat(3,MIX3respsImmobile, allresults{1,i}{1,j}.expRespsImmobile);
                MIX3resps=cat(3,MIX3resps,allresults{1,i}{1,j}.resps);
                MIX3respsLocomotionByStim=cat(3,MIX3respsLocomotionByStim,allresults{1,i}{1,j}.respsLocomotionByStim);
                
                alldffMIX3=[alldff3 allresults{1,i}{1,j}.dff];
                allmovMIX3=[allmov3 allresults{1,i}{1,j}.movementTrace];
                
            end
            
            if strcmp(allresults{1,i}{1,j}.stim,'TON3')
                TON3respsMOV=cat(3,TON3respsMOV, allresults{1,i}{1,j}.expRespsMove);
                TON3respsImmobile=cat(3,TON3respsImmobile, allresults{1,i}{1,j}.expRespsImmobile);
                TON3resps=cat(3,TON3resps,allresults{1,i}{1,j}.resps);
                TON3respsLocomotionByStim=cat(3,TON3respsLocomotionByStim,allresults{1,i}{1,j}.respsLocomotionByStim);
                
                alldffTON3=[alldff3 allresults{1,i}{1,j}.dff];
                allmovTON3=[allmov3 allresults{1,i}{1,j}.movementTrace];
                
            end
            
            if strcmp(allresults{1,i}{1,j}.stim,'NAT3')
                NAT3respsMOV=cat(3,NAT3respsMOV, allresults{1,i}{1,j}.expRespsMove);
                NAT3respsImmobile=cat(3,NAT3respsImmobile, allresults{1,i}{1,j}.expRespsImmobile);
                NAT3resps=cat(3,NAT3resps,allresults{1,i}{1,j}.resps);
                NAT3respsLocomotionByStim=cat(3,NAT3respsLocomotionByStim,allresults{1,i}{1,j}.respsLocomotionByStim);
                
                alldff3=[alldff3 allresults{1,i}{1,j}.dff];
                allmov3=[allmov3 allresults{1,i}{1,j}.movementTrace];
                
            end
            
            if strcmp(allresults{1,i}{1,j}.stim,'NAT')|strcmp(allresults{1,i}{1,j}.stim,'NAT2')
                NATrespsMOV=cat(3,NATrespsMOV, allresults{1,i}{1,j}.expRespsMove);
                NATrespsImmobile=cat(3,NATrespsImmobile, allresults{1,i}{1,j}.expRespsImmobile);
                NATresps=cat(3,NATresps,allresults{1,i}{1,j}.resps);
                NATrespsLocomotionByStim=cat(3,NATrespsLocomotionByStim,allresults{1,i}{1,j}.respsLocomotionByStim);
                
                alldff=[alldff allresults{1,i}{1,j}.dff];
                allmov=[allmov allresults{1,i}{1,j}.movementTrace];
                
            end
            if strcmp(allresults{1,i}{1,j}.stim,'BBN')|strcmp(allresults{1,i}{1,j}.stim,'BBN2')
                
                
                if size(allresults{1,i}{1,j}.resps,4)==size(BBNresps,4)+1
                    allresults{1,i}{1,j}.resps=allresults{1,i}{1,j}.resps(:,:,:,1:end-1);
                    allresults{1,i}{1,j}.respsLocomotionByStim=allresults{1,i}{1,j}.respsLocomotionByStim(:,:,:,1:end-1);
                    
                    allresults{1,i}{1,j}.expRespsMove=allresults{1,i}{1,j}.expRespsMove(:,:,:,1:end-1);
                    allresults{1,i}{1,j}.expRespsImmobile=allresults{1,i}{1,j}.expRespsImmobile(:,:,:,1:end-1);
                end
                BBNresps=cat(3,BBNresps,allresults{1,i}{1,j}.resps);
                BBNrespsLocomotionByStim=cat(3,BBNrespsLocomotionByStim,allresults{1,i}{1,j}.respsLocomotionByStim);
                
                BBNrespsMOV=cat(3,BBNrespsMOV, allresults{1,i}{1,j}.expRespsMove);
                BBNrespsImmobile=cat(3,BBNrespsImmobile, allresults{1,i}{1,j}.expRespsImmobile);
                
                alldff=[alldff allresults{1,i}{1,j}.dff];
                allmov=[allmov allresults{1,i}{1,j}.movementTrace];
                
            end
            if strcmp(allresults{1,i}{1,j}.stim,'NONE')
                
                
                alldffNONE=[alldffNONE allresults{1,i}{1,j}.dff];
                allmovNONE=[allmovNONE allresults{1,i}{1,j}.movementTrace];
                acrossSessionsSpeedsNONE=[acrossSessionsSpeedsNONE allresults{1,i}{1,j}.locomotionPerFrame'];
                
            end
            if strcmp(allresults{1,i}{1,j}.stim,'CONST')
                alldffCONST=[alldffCONST allresults{1,i}{1,j}.dff];
                allmovCONST=[allmovCONST allresults{1,i}{1,j}.movementTrace];
                acrossSessionsSpeedsCONST=[acrossSessionsSpeedsCONST allresults{1,i}{1,j}.locomotionPerFrame'];
                
            end
            if ~isempty(allresults{1,i}{1,j}.locomotionPerFrame)&~strcmp(allresults{1,i}{1,j}.stim,'CONST')&~strcmp(allresults{1,i}{1,j}.stim,'NONE')
                
                acrossSessionsSpeeds=[acrossSessionsSpeeds allresults{1,i}{1,j}.locomotionPerFrame'];
                allallSpeecCorrVals=[allallSpeecCorrVals allresults{1,i}{1,j}.allSpeecCorrVals(:,1)];
            else
                acrossSessionsSpeeds=[acrossSessionsSpeeds];
                allallSpeecCorrVals=[allallSpeecCorrVals];
            end
        end
    end
    
    
    %++++++++++++++++++++++++++BBN ANALYSIS+++++++++++++++++++++++++++++
    
    stimwinlength=4;
    prestimwin=floor((size(BBNresps,4)/2)-stimwinlength:(size(BBNresps,4)/2)-1);
    stimwin=floor((size(BBNresps,4)/2)+1:(size(BBNresps,4)/2)+stimwinlength);
    
    
    numcells=size(BBNrespsMOV,1);
    difSpontLocomotionVsImmobility=[];
    difEvokedLocomotionVsImmobility=[];
    difEvokedMinusPreLocomotionVsImmobility=[];
    sValsS_onlyResponsive=[];
    cellResponsiveness=[];
    cellResponsivenessImmobile=[];
    ensembleNormSoundEvokedResps=[];
    ensembleNormSoundEvokedRespsMove=[];
    ensembleNormSoundEvokedRespsImmob=[];
    
    minnumtrialsperstate=8;
    for c=1:numcells
        
        % Responses across locomotion and immobility
        curCellMOVandImmobility=squeeze(BBNresps(c,1,:,:));
        notAllNANrowsMOVandImmobility=find(sum(~isnan(curCellMOVandImmobility)')~=0);
        curCellMOVandImmobilitynoNAN=curCellMOVandImmobility(notAllNANrowsMOVandImmobility,:);
        
        curCellRespBySpeed=squeeze(BBNrespsLocomotionByStim(c,1,:,:));
        curCellRespBySpeed=nanmean(curCellRespBySpeed(notAllNANrowsMOVandImmobility,:),2);
        [sortedSpeeds sortedSpeedsInds]=sort(curCellRespBySpeed);
        
        curCellEvokedResps=nanmean(curCellMOVandImmobilitynoNAN(:,stimwin),2);
        curCellPreResps=nanmean(curCellMOVandImmobilitynoNAN(:,prestimwin),2);
        
        goodinds=find(~isnan(curCellPreResps)&~isnan(curCellEvokedResps));
        curCellPreResps=curCellPreResps(goodinds);
        curCellEvokedResps=curCellEvokedResps(goodinds);
        
        % DETERMINING SOUND RESPONSIVENESS
        [curCellResponsivenessH, curCellResponsivenessP]=ttest(curCellPreResps,curCellEvokedResps,'Tail','left');
        cellResponsiveness=[cellResponsiveness;curCellResponsivenessH];
        
        ensembleNormSoundEvokedResps=[ensembleNormSoundEvokedResps;(curCellEvokedResps-curCellPreResps)'];
        
        if plotcellsresponsiveness
            figure;
            subplot(2,1,1);
            errorbar(nanmean(curCellMOVandImmobilitynoNAN),nanstd(curCellMOVandImmobilitynoNAN)/sqrt(size(curCellMOVandImmobilitynoNAN,1)),'k','linewidth',2);
            axis([1 13 min(0,min(nanmean(curCellMOVandImmobilitynoNAN))) max(1,1.2*max(nanmean(curCellMOVandImmobilitynoNAN)))]);
            title(['Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c) ' P= ' num2str(curCellResponsivenessP)]);
            subplot(2,1,2);
            imagesc(curCellMOVandImmobilitynoNAN)
            savecellresponsiveness=1;
            if savecellresponsiveness
                saveas(gcf,[cellsresponsivenessFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'jpg');
                saveas(gcf,[cellsresponsivenessFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'fig');
                print([cellsresponsivenessFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'-dpdf','-bestfit')
                
            end
            close all
        end
        plotSpeedSorting=0;
        if plotSpeedSorting
            colormap copper
            figure('Position',[50,50,400, 800]);
            subplot(3,2,1)
            imagesc(curCellMOVandImmobilitynoNAN)
            xlabel('Time (frames)')
            ylabel('Trial #');
            colormap copper
            
            subplot(3,2,2)
            plot(curCellRespBySpeed);
            axis([1 length(curCellRespBySpeed) 0 60])
            camroll(-90);
            ylabel('Speed')
            
            
            subplot(3,2,3)
            imagesc(curCellMOVandImmobilitynoNAN(sortedSpeedsInds,:))
            xlabel('Time (frames)')
            ylabel('Trial #');
            colormap copper
            
            subplot(3,2,4)
            plot(sortedSpeeds);
            
            axis([1 length(curCellRespBySpeed) 0 60])
            camroll(-90);
            ylabel('Speed')
            
            subplot(3,2,5)
            plot(nanmean(curCellMOVandImmobilitynoNAN(curCellRespBySpeed<4,:)),'g','linewidth',2)
            hold on
            try
                shadedErrorBar(1:size(curCellMOVandImmobilitynoNAN,2),nanmean(curCellMOVandImmobilitynoNAN(curCellRespBySpeed<4,:)),nanstd(curCellMOVandImmobilitynoNAN(curCellRespBySpeed<4,:))/sqrt(size(curCellMOVandImmobilitynoNAN(curCellRespBySpeed<4,:),1)),'lineprops','g')
            catch
            end
            plot(nanmean(curCellMOVandImmobilitynoNAN(curCellRespBySpeed>=4,:)),'r','linewidth',2)
            try
                shadedErrorBar(1:size(curCellMOVandImmobilitynoNAN,2),nanmean(curCellMOVandImmobilitynoNAN(curCellRespBySpeed>=4,:)),nanstd(curCellMOVandImmobilitynoNAN(curCellRespBySpeed>=4,:))/sqrt(size(curCellMOVandImmobilitynoNAN(curCellRespBySpeed>=4,:),1)),'lineprops','r')
            catch
            end
            xlim([1 size(curCellMOVandImmobilitynoNAN,2)]);
            curY=ylim;
            ylim([min(0,curY(1)),max(1,curY(2))])
            xlabel('Time (frames)')
            ylabel('dFF')
            savecellresponsivenessbyspeed=1;
            if savecellresponsivenessbyspeed
                saveas(gcf,[cellsresponsivenessbyspeedFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'jpg');
                saveas(gcf,[cellsresponsivenessbyspeedFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'fig');
                print([cellsresponsivenessbyspeedFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'-dpdf','-bestfit')
            end
            close all
        end
        
        
        
        % LOCOMOTION RESPONSES
        curCellMOV=squeeze(BBNrespsMOV(c,1,:,:));
        notAllNANrowsMOV=find(sum(~isnan(curCellMOV)')~=0);
        curCellMOVnoNAN=curCellMOV(notAllNANrowsMOV,:);
        
        % IMMOBILITY RESPONSES
        curCellImmobile=squeeze(BBNrespsImmobile(c,1,:,:));
        notAllNANrowsImmobile=find(sum(~isnan(curCellImmobile)')~=0);
        curCellImmobilenoNAN=curCellImmobile(notAllNANrowsImmobile,:);
        
        
        % EVOKED RESPONSES DURING LOCOMOTION
        curCellMove_evoked =nanmean(curCellMOVnoNAN(:,stimwin),2);
        
        % EVOKED RESPONSES DURING IMMOBILITY
        curCellImmobile_evoked=nanmean(curCellImmobilenoNAN(:,stimwin),2);
        
        % PRE-STIM ACTIVITY DURING LOCOMOTION
        curCellMove_pre =nanmean(curCellMOVnoNAN(:,prestimwin),2);
        
        % PRE-STIM ACTIVITY DURING IMMOBILITY
        curCellImmobile_pre=nanmean(curCellImmobilenoNAN(:,prestimwin),2);
        
        % EVOKED MINUS PRE-STIM DURING LOCOMOTION
        evokedMinusPre_Move=curCellMove_evoked-curCellMove_pre;
        
        % EVOKED MINUS PRE-STIM DURING IMMOBILITY
        evokedMinusPre_Immobile=curCellImmobile_evoked-curCellImmobile_pre;
        
        [curCellResponsivenessImmobileH, curCellResponsivenessImmobileP]=ttest(curCellImmobile_pre,curCellImmobile_evoked,'Tail','left');
        cellResponsivenessImmobile=[cellResponsivenessImmobile;curCellResponsivenessImmobileH];
        
        ensembleNormSoundEvokedRespsMove=[ensembleNormSoundEvokedRespsMove;evokedMinusPre_Move'];
        ensembleNormSoundEvokedRespsImmob=[ensembleNormSoundEvokedRespsImmob;evokedMinusPre_Immobile'];
        
        
        if length(evokedMinusPre_Move)>minnumtrialsperstate & length(evokedMinusPre_Immobile)>minnumtrialsperstate
            
            
          
            [h, p]=ttest2(curCellMove_pre,curCellImmobile_pre);
            sigModulation=0;
            if p<0.05
                if mean(curCellMove_pre)>mean(curCellImmobile_pre)
                    sigModulation=1;
                else
                    sigModulation=-1;
                end
            end
            
            difSpontLocomotionVsImmobility=[difSpontLocomotionVsImmobility;sigModulation];
            allmeanspontLocomotion=[allmeanspontLocomotion;mean(curCellMove_pre)];
            allmeanspontImmobility=[allmeanspontImmobility;mean(curCellImmobile_pre)];
            difSpontLocomotionVsImmobilityAllExps=[difSpontLocomotionVsImmobilityAllExps;sigModulation];
            
            if curCellResponsivenessH
                allmeanspontLocomotion_responsive=[allmeanspontLocomotion_responsive;mean(curCellMove_pre)];
                allmeanspontImmobility_responsive=[allmeanspontImmobility_responsive;mean(curCellImmobile_pre)];
                difSpontLocomotionVsImmobilityAllExps_onlyResponsive=[difSpontLocomotionVsImmobilityAllExps_onlyResponsive; sigModulation];
                
                if sum(isnan(allmeanspontImmobility_responsive))>0
                    keyboard
                end
            end
            if curCellResponsivenessImmobileH
                allmeanevokedLocomotion_responsiveIMMOB=[allmeanevokedLocomotion_responsiveIMMOB;mean(curCellMove_evoked)];
                allmeanevokedImmobility_responsiveIMMOB =[allmeanevokedImmobility_responsiveIMMOB;mean(curCellImmobile_evoked)];
                
                allmeanspontLocomotion_responsiveIMMOB=[allmeanspontLocomotion_responsiveIMMOB;mean(curCellMove_pre)];
                allmeanspontImmobility_responsiveIMMOB=[allmeanspontImmobility_responsiveIMMOB;mean(curCellImmobile_pre)];
                difSpontLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB=[difSpontLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB; sigModulation];
                
            end
            
      
            [h, p]=ttest2(curCellMove_evoked,curCellImmobile_evoked);
            sigModulation=0;
            if p<0.05
                if mean(curCellMove_evoked)>mean(curCellImmobile_evoked)
                    sigModulation=1;
                else
                    sigModulation=-1;
                end
                
            end
            
            difEvokedLocomotionVsImmobility=[difEvokedLocomotionVsImmobility;sigModulation];
            allmeanevokedLocomotion=[allmeanevokedLocomotion;mean(curCellMove_evoked)];
            allmeanevokedImmobility=[allmeanevokedImmobility;mean(curCellImmobile_evoked)];
            difEvokedLocomotionVsImmobilityAllExps=[difEvokedLocomotionVsImmobilityAllExps;sigModulation];
            
            if curCellResponsivenessH
                allmeanevokedLocomotion_responsive=[allmeanevokedLocomotion_responsive;mean(curCellMove_evoked)];
                allmeanevokedImmobility_responsive=[allmeanevokedImmobility_responsive;mean(curCellImmobile_evoked)];
                difEvokedLocomotionVsImmobilityAllExps_onlyResponsive=[difEvokedLocomotionVsImmobilityAllExps_onlyResponsive; sigModulation];
                
            end
            
         
            [h, p]=ttest2(evokedMinusPre_Move,evokedMinusPre_Immobile);
            sigModulation=0;
            if p<0.05
                if mean(evokedMinusPre_Move)>mean(evokedMinusPre_Immobile)
                    sigModulation=1;
                else
                    sigModulation=-1;
                end
                
            end
            
            difEvokedMinusPreLocomotionVsImmobility=[difEvokedMinusPreLocomotionVsImmobility;sigModulation];
            allmeanEvokedMinusPreLocomotion=[allmeanEvokedMinusPreLocomotion;mean(evokedMinusPre_Move)];
            allmeanEvokedMinusPreImmobility=[allmeanEvokedMinusPreImmobility;mean(evokedMinusPre_Immobile)];
            
            difEvokedMinusPreLocomotionVsImmobilityAllExps=[difEvokedMinusPreLocomotionVsImmobilityAllExps;sigModulation];
            curcellstr={[allresults{1,i}{1,1}.batch ' Cell ' num2str(c)]};
            allinds2{end+1}=curcellstr;
            %
            if curCellResponsivenessH
                allmeanEvokedMinusPreLocomotion_responsive=[ allmeanEvokedMinusPreLocomotion_responsive;mean(evokedMinusPre_Move)];
                allmeanEvokedMinusPreImmobility_responsive=[allmeanEvokedMinusPreImmobility_responsive;mean(evokedMinusPre_Immobile)];
                difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsive=[difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsive; sigModulation];
                curcellstr={[allresults{1,i}{1,1}.batch ' Cell ' num2str(c)]};
                responsiveinds{end+1}=curcellstr;
                
            end
            if curCellResponsivenessImmobileH
                allmeanEvokedMinusPreLocomotion_responsiveIMMOB=[ allmeanEvokedMinusPreLocomotion_responsiveIMMOB;mean(evokedMinusPre_Move)];
                allmeanEvokedMinusPreImmobility_responsiveIMMOB=[allmeanEvokedMinusPreImmobility_responsiveIMMOB;mean(evokedMinusPre_Immobile)];
                difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB=[difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB; sigModulation];
                immobresponsiveinds{end+1}=curcellstr;
                movPSTH=nanmean(curCellMOVnoNAN);
                immobPSTH=nanmean(curCellImmobilenoNAN);
                allMovePSTHsResponsiveIMMOB(end+1,1:length(movPSTH))=movPSTH;
                allImmobPSTHsResponsiveIMMOB(end+1,1:length(immobPSTH))=immobPSTH;
                
            end
            %==========
            cellsresponsivenessMovVsImmobFigDir='E:\GideonData\IMAGING\responsivenessMoveVsImmob3\';
            plotcellsresponsivenessMovVsImmob=0;
            if plotcellsresponsivenessMovVsImmob
                if curCellResponsivenessH | curCellResponsivenessImmobileH
                    figure('Position', [50, 50, 400, 500])
                    %                    figure('Position', [100, 100, 400, 600])
                    % subplot(3,1,1)
                    shadedErrorBar(1:size(curCellImmobilenoNAN,2),nanmean(curCellImmobilenoNAN),nanstd(curCellImmobilenoNAN)/sqrt(size(curCellImmobilenoNAN,1)),'lineprops','r')
                    hold on
                    shadedErrorBar(1:size(curCellMOVnoNAN,2),nanmean(curCellMOVnoNAN),nanstd(curCellMOVnoNAN)/sqrt(size(curCellMOVnoNAN,1)),'lineprops','g')
                    ht=title(['Cell ' allresults{1,i}{1,1}.batch '-' num2str(c) ' Resp= ' num2str(curCellResponsivenessH) ' ImResp = ' num2str(curCellResponsivenessImmobileH) ' S/N: ' num2str(sigModulation)]);
                    set(ht,'FontSize',8)
                    xlim([1 13])
                    set(gca,'FontSize',16)
                    set(ht,'FontSize',8)
                    
                    
                    savecellresponsivenessMovVsImmob=1;
                    if savecellresponsivenessMovVsImmob
                        saveas(gcf,[cellsresponsivenessMovVsImmobFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'jpg');
                        % saveas(gcf,[cellsresponsivenessMovVsImmobFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'fig');
                        print([cellsresponsivenessMovVsImmobFigDir 'Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c)],'-dpdf','-bestfit')
                        
                    end
                    close all
                end
            end
            %==========
        end
    end % end of cell-loop of current ensemble
    
    
    % Examining cross-correlations between cells that show locomotion-induced increase or decrease
    % in sound-evoked responses
    
    if ~isempty(difEvokedMinusPreLocomotionVsImmobility)
        
        locomotionEnhancedCells=find(difEvokedMinusPreLocomotionVsImmobility==1);
        locomotionSuppressedCells=find(difEvokedMinusPreLocomotionVsImmobility==-1);
        locomotionUnmodCells=find(difEvokedMinusPreLocomotionVsImmobility==0);
        
        locomotionEnhancedPairs=[];
        locomotionSuppressedPairs=[];
        locomotionUnmodPairs=[];
        acrossClassPairs=[];
        
        % xcorr between locomotion-enhanced pairs
        if length(locomotionEnhancedCells)>1
            locomotionEnhancedPairs=nchoosek(locomotionEnhancedCells,2);
            locomotionEnhancedPairsXCorrs=[];
            locomotionEnhancedPairsXCorrsImmob=[];
            locomotionEnhancedPairsXCorrsMove=[];
            
            locomotionEnhancedPairsXCorrsShuf=[];
            locomotionEnhancedPairsXCorrsImmobShuf=[];
            locomotionEnhancedPairsXCorrsMoveShuf=[];
            
            for w=1:size(locomotionEnhancedPairs,1)
                xc=xcorr(alldff(locomotionEnhancedPairs(w,1),:),alldff(locomotionEnhancedPairs(w,2),:),20,'coeff');
                locomotionEnhancedPairsXCorrs=[locomotionEnhancedPairsXCorrs;xc];
        
                xc_immob=xcorr(alldff(locomotionEnhancedPairs(w,1),allmov==0),alldff(locomotionEnhancedPairs(w,2),allmov==0),20,'coeff');
                locomotionEnhancedPairsXCorrsImmob=[locomotionEnhancedPairsXCorrsImmob;xc_immob];
                xc_move=xcorr(alldff(locomotionEnhancedPairs(w,1),allmov==1),alldff(locomotionEnhancedPairs(w,2),allmov==1),20,'coeff');
                locomotionEnhancedPairsXCorrsMove=[locomotionEnhancedPairsXCorrsMove;xc_move];
                
                
                % shuffling
                f1=alldff(locomotionEnhancedPairs(w,1),:);
                f2=alldff(locomotionEnhancedPairs(w,2),:);
                f1_immob=alldff(locomotionEnhancedPairs(w,1),allmov==0);
                f2_immob=alldff(locomotionEnhancedPairs(w,2),allmov==0);
                f1_move=alldff(locomotionEnhancedPairs(w,1),allmov==1);
                f2_move=alldff(locomotionEnhancedPairs(w,2),allmov==1);
                
                numcircshufs=100;
                allxcs=[];
                allxcs_immob=[];
                allxcs_move=[];
                for kk=1:100
                    circshiftamount=round(rand(1)*size(alldff,2));
                    f1s=circshift(f1,circshiftamount);
                    xcs1=xcorr(f1s,f2,20,'coeff');
                    allxcs=[allxcs; xcs1];
                    
                    f1s_immob=circshift(f1_immob,circshiftamount);
                    xcs_immob1=xcorr(f1s_immob,f2_immob,20,'coeff');
                    allxcs_immob=[allxcs_immob; xcs_immob1];
                    
                    f1s_move=circshift(f1_move,circshiftamount);
                    xcs_move1=xcorr(f1s_move,f2_move,20,'coeff');
                    allxcs_move=[allxcs_move; xcs_move1];
                end
                
                xcs=mean(allxcs);
                xcs_immob=mean(allxcs_immob);
                xcs_move=mean(allxcs_move);
                
                locomotionEnhancedPairsXCorrsShuf=[locomotionEnhancedPairsXCorrsShuf;xcs];
                locomotionEnhancedPairsXCorrsImmobShuf=[locomotionEnhancedPairsXCorrsImmobShuf;xcs_immob];
                locomotionEnhancedPairsXCorrsMoveShuf=[locomotionEnhancedPairsXCorrsMoveShuf;xcs_move];
                
            end
            alllocomotionEnhancedPairsXCorrs=[alllocomotionEnhancedPairsXCorrs;locomotionEnhancedPairsXCorrs];
            alllocomotionEnhancedPairsXCorrsImmob=[alllocomotionEnhancedPairsXCorrsImmob;locomotionEnhancedPairsXCorrsImmob];
            alllocomotionEnhancedPairsXCorrsMove=[alllocomotionEnhancedPairsXCorrsMove;locomotionEnhancedPairsXCorrsMove];
            
            alllocomotionEnhancedPairsXCorrsShuf=[alllocomotionEnhancedPairsXCorrsShuf;locomotionEnhancedPairsXCorrsShuf];
            alllocomotionEnhancedPairsXCorrsImmobShuf=[alllocomotionEnhancedPairsXCorrsImmobShuf;locomotionEnhancedPairsXCorrsImmobShuf];
            alllocomotionEnhancedPairsXCorrsMoveShuf=[alllocomotionEnhancedPairsXCorrsMoveShuf;locomotionEnhancedPairsXCorrsMoveShuf];
            
            
        end
        
        % xcorr between locomotion-suppressed pairs
        if length(locomotionSuppressedCells)>1
            locomotionSuppressedPairs=nchoosek(locomotionSuppressedCells,2);
            locomotionSuppressedPairsXCorrs=[];
            locomotionSuppressedPairsXCorrsImmob=[];
            locomotionSuppressedPairsXCorrsMove=[];
            
            locomotionSuppressedPairsXCorrsShuf=[];
            locomotionSuppressedPairsXCorrsImmobShuf=[];
            locomotionSuppressedPairsXCorrsMoveShuf=[];
            
            for w=1:size(locomotionSuppressedPairs,1)
                xc=xcorr(alldff(locomotionSuppressedPairs(w,1),:),alldff(locomotionSuppressedPairs(w,2),:),20,'coeff');
                locomotionSuppressedPairsXCorrs=[locomotionSuppressedPairsXCorrs;xc];
           
                xc_immob=xcorr(alldff(locomotionSuppressedPairs(w,1),allmov==0),alldff(locomotionSuppressedPairs(w,2),allmov==0),20,'coeff');
                locomotionSuppressedPairsXCorrsImmob=[locomotionSuppressedPairsXCorrsImmob;xc_immob];
                xc_move=xcorr(alldff(locomotionSuppressedPairs(w,1),allmov==1),alldff(locomotionSuppressedPairs(w,2),allmov==1),20,'coeff');
                locomotionSuppressedPairsXCorrsMove=[locomotionSuppressedPairsXCorrsMove;xc_move];
                
                f1=alldff(locomotionSuppressedPairs(w,1),:);
                f2=alldff(locomotionSuppressedPairs(w,2),:);
                f1_immob=alldff(locomotionSuppressedPairs(w,1),allmov==0);
                f2_immob=alldff(locomotionSuppressedPairs(w,2),allmov==0);
                f1_move=alldff(locomotionSuppressedPairs(w,1),allmov==1);
                f2_move=alldff(locomotionSuppressedPairs(w,2),allmov==1);
                
                numcircshufs=100;
                allxcs=[];
                allxcs_immob=[];
                allxcs_move=[];
                for kk=1:100
                    circshiftamount=round(rand(1)*size(alldff,2));
                    f1s=circshift(f1,circshiftamount);
                    xcs1=xcorr(f1s,f2,20,'coeff');
                    allxcs=[allxcs; xcs1];
                    
                    f1s_immob=circshift(f1_immob,circshiftamount);
                    xcs_immob1=xcorr(f1s_immob,f2_immob,20,'coeff');
                    allxcs_immob=[allxcs_immob; xcs_immob1];
                    
                    f1s_move=circshift(f1_move,circshiftamount);
                    xcs_move1=xcorr(f1s_move,f2_move,20,'coeff');
                    allxcs_move=[allxcs_move; xcs_move1];
                end
                
                xcs=mean(allxcs);
                xcs_immob=mean(allxcs_immob);
                xcs_move=mean(allxcs_move);
                
                locomotionSuppressedPairsXCorrsShuf=[locomotionSuppressedPairsXCorrsShuf;xcs];
                locomotionSuppressedPairsXCorrsImmobShuf=[locomotionSuppressedPairsXCorrsImmobShuf;xcs_immob];
                locomotionSuppressedPairsXCorrsMoveShuf=[locomotionSuppressedPairsXCorrsMoveShuf;xcs_move];
                
            end
            alllocomotionSuppressedPairsXCorrs=[alllocomotionSuppressedPairsXCorrs;locomotionSuppressedPairsXCorrs];
            alllocomotionSuppressedPairsXCorrsImmob=[alllocomotionSuppressedPairsXCorrsImmob;locomotionSuppressedPairsXCorrsImmob];
            alllocomotionSuppressedPairsXCorrsMove=[alllocomotionSuppressedPairsXCorrsMove;locomotionSuppressedPairsXCorrsMove];
            
            alllocomotionSuppressedPairsXCorrsShuf=[alllocomotionSuppressedPairsXCorrsShuf;locomotionSuppressedPairsXCorrsShuf];
            alllocomotionSuppressedPairsXCorrsImmobShuf=[alllocomotionSuppressedPairsXCorrsImmobShuf;locomotionSuppressedPairsXCorrsImmobShuf];
            alllocomotionSuppressedPairsXCorrsMoveShuf=[alllocomotionSuppressedPairsXCorrsMoveShuf;locomotionSuppressedPairsXCorrsMoveShuf];
            
        end
        
        % xcorr between locomotion-unmod pairs
        if length(locomotionUnmodCells)>1
            locomotionUnmodPairs=nchoosek(locomotionUnmodCells,2);
            locomotionUnmodPairsXCorrs=[];
            for w=1:size(locomotionUnmodPairs,1)
                xc=xcorr(alldff(locomotionUnmodPairs(w,1),:),alldff(locomotionUnmodPairs(w,2),:),20,'coeff');
                locomotionUnmodPairsXCorrs=[locomotionUnmodPairsXCorrs;xc];
            end
            alllocomotionUnmodPairsXCorrs=[alllocomotionUnmodPairsXCorrs;locomotionUnmodPairsXCorrs];
            
        end
        
        % xcorr between pairs where one is locomotion-enhanced and one is
        % locomotion-suppressed
        if length(locomotionEnhancedCells)>1 & length(locomotionSuppressedCells)>1
 
            [A B]=meshgrid(locomotionEnhancedCells,locomotionSuppressedCells);
            c=cat(2,A',B');
            acrossClassPairs=reshape(c,[],2);
            
            acrossClassesXCorrs=[];
            acrossClassesXCorrsImmob=[];
            acrossClassesXCorrsMove=[];
            
            acrossClassesXCorrsShuf=[];
            acrossClassesXCorrsImmobShuf=[];
            acrossClassesXCorrsMoveShuf=[];
            
            for w=1:size(acrossClassPairs,1)
                xc=xcorr(alldff(acrossClassPairs(w,1),:),alldff(acrossClassPairs(w,2),:),20,'coeff');
                acrossClassesXCorrs=[acrossClassesXCorrs;xc];
              
                xc_immob=xcorr(alldff(acrossClassPairs(w,1),allmov==0),alldff(acrossClassPairs(w,2),allmov==0),20,'coeff');
                acrossClassesXCorrsImmob=[acrossClassesXCorrsImmob;xc_immob];
                xc_move=xcorr(alldff(acrossClassPairs(w,1),allmov==1),alldff(acrossClassPairs(w,2),allmov==1),20,'coeff');
                acrossClassesXCorrsMove=[acrossClassesXCorrsMove;xc_move];
                
                
             
                f1=alldff(acrossClassPairs(w,1),:);
                f2=alldff(acrossClassPairs(w,2),:);
                f1_immob=alldff(acrossClassPairs(w,1),allmov==0);
                f2_immob=alldff(acrossClassPairs(w,2),allmov==0);
                f1_move=alldff(acrossClassPairs(w,1),allmov==1);
                f2_move=alldff(acrossClassPairs(w,2),allmov==1);
                
                numcircshufs=100;
                allxcs=[];
                allxcs_immob=[];
                allxcs_move=[];
                for kk=1:100
                    circshiftamount=round(rand(1)*size(alldff,2));
                    f1s=circshift(f1,circshiftamount);
                    xcs1=xcorr(f1s,f2,20,'coeff');
                    allxcs=[allxcs; xcs1];
                    
                    f1s_immob=circshift(f1_immob,circshiftamount);
                    xcs_immob1=xcorr(f1s_immob,f2_immob,20,'coeff');
                    allxcs_immob=[allxcs_immob; xcs_immob1];
                    
                    f1s_move=circshift(f1_move,circshiftamount);
                    xcs_move1=xcorr(f1s_move,f2_move,20,'coeff');
                    allxcs_move=[allxcs_move; xcs_move1];
                end
                
                xcs=mean(allxcs);
                xcs_immob=mean(allxcs_immob);
                xcs_move=mean(allxcs_move);
                
                acrossClassesXCorrsShuf=[acrossClassesXCorrsShuf;xcs];
                acrossClassesXCorrsImmobShuf=[acrossClassesXCorrsImmobShuf;xcs_immob];
                acrossClassesXCorrsMoveShuf=[acrossClassesXCorrsShuf;xcs_move];
                
            end
            
            allacrossClassesXCorrs=[allacrossClassesXCorrs;acrossClassesXCorrs];
            allacrossClassesXCorrsImmob=[allacrossClassesXCorrsImmob;acrossClassesXCorrsImmob];
            allacrossClassesXCorrsMove=[allacrossClassesXCorrsMove;acrossClassesXCorrsMove];
            
            allacrossClassesXCorrsShuf=[allacrossClassesXCorrsShuf;acrossClassesXCorrsShuf];
            allacrossClassesXCorrsImmobShuf=[allacrossClassesXCorrsImmobShuf;acrossClassesXCorrsImmobShuf];
            allacrossClassesXCorrsMoveShuf=[allacrossClassesXCorrsMoveShuf;acrossClassesXCorrsMoveShuf];
            
        end
        
        
        % Examining noise-correlations between cells that show locomotion-induced increase or decrease
        % in sound-evoked responses
        [noiseCorrs ppp]=corrcoef(ensembleNormSoundEvokedResps');
        [noiseCorrsMove pppMove]=corrcoef(ensembleNormSoundEvokedRespsMove');
        [noiseCorrsImmob pppImmob]=corrcoef(ensembleNormSoundEvokedRespsImmob');
        
    
        numshufsnoisecorr=100;
        allnoiseCorrsShuf=[];
        allnoiseCorrsMoveShuf=[];
        allnoiseCorrsImmobShuf=[];
        
        for zz=1:numshufsnoisecorr
            ensembleNormSoundEvokedRespsS=[];
            for kk=1:size(ensembleNormSoundEvokedResps,1)
                ensembleNormSoundEvokedRespsS(kk,:)=circshift(ensembleNormSoundEvokedResps(kk,:),round(rand(1)*size(ensembleNormSoundEvokedResps,2)));
            end
            [noiseCorrsShuf1 pppShuf1]=corrcoef(ensembleNormSoundEvokedRespsS');
            allnoiseCorrsShuf(zz,:,:)=noiseCorrsShuf1;
            
            ensembleNormSoundEvokedRespsMoveS=[];
            for kk=1:size(ensembleNormSoundEvokedRespsMove,1)
                ensembleNormSoundEvokedRespsMoveS(kk,:)=circshift(ensembleNormSoundEvokedRespsMove(kk,:),round(rand(1)*size(ensembleNormSoundEvokedRespsMove,2)));
            end
            [noiseCorrsMoveShuf1 pppMoveShuf1]=corrcoef(ensembleNormSoundEvokedRespsMoveS');
            allnoiseCorrsMoveShuf(zz,:,:)=noiseCorrsMoveShuf1;
            
            
            ensembleNormSoundEvokedRespsImmobS=[];
            for kk=1:size(ensembleNormSoundEvokedRespsImmob,1)
                ensembleNormSoundEvokedRespsImmobS(kk,:)=circshift(ensembleNormSoundEvokedRespsImmob(kk,:),round(rand(1)*size(ensembleNormSoundEvokedRespsImmob,2)));
            end
          [noiseCorrsImmobShuf1 pppImmobShuf1]=corrcoef(ensembleNormSoundEvokedRespsImmobS');
            allnoiseCorrsImmobShuf(zz,:,:)=noiseCorrsImmobShuf1;
            
            
        end
        noiseCorrsShuf=squeeze(mean(allnoiseCorrsShuf,1));
        noiseCorrsMoveShuf=squeeze(mean(allnoiseCorrsMoveShuf,1));
        noiseCorrsImmobShuf=squeeze(mean(allnoiseCorrsImmobShuf,1));
        
        
        locomotionEnhancedNoiseCorrs=[];
        locomotionEnhancedNoiseCorrsMove=[];
        locomotionEnhancedNoiseCorrsImmob=[];
        
        locomotionSuppressedNoiseCorrs=[];
        locomotionSuppressedNoiseCorrsMove=[];
        locomotionSuppressedNoiseCorrsImmob=[];
        
        acrossClassesNoiseCorrs=[];
        acrossClassesNoiseCorrsMove=[];
        acrossClassesNoiseCorrsImmob=[];
        
        
        locomotionEnhancedNoiseCorrsShuf=[];
        locomotionEnhancedNoiseCorrsMoveShuf=[];
        locomotionEnhancedNoiseCorrsImmobShuf=[];
        
        locomotionSuppressedNoiseCorrsShuf=[];
        locomotionSuppressedNoiseCorrsMoveShuf=[];
        locomotionSuppressedNoiseCorrsImmobShuf=[];
        
        acrossClassesNoiseCorrsShuf=[];
        acrossClassesNoiseCorrsMoveShuf=[];
        acrossClassesNoiseCorrsImmobShuf=[];
        
        locomotionUnmodNoiseCorrs=[];
        for w=1:size(locomotionEnhancedPairs,1)
            locomotionEnhancedNoiseCorrs=[locomotionEnhancedNoiseCorrs noiseCorrs(locomotionEnhancedPairs(w,1),locomotionEnhancedPairs(w,2))];
            locomotionEnhancedNoiseCorrsMove=[locomotionEnhancedNoiseCorrsMove noiseCorrsMove(locomotionEnhancedPairs(w,1),locomotionEnhancedPairs(w,2))];
            locomotionEnhancedNoiseCorrsImmob=[locomotionEnhancedNoiseCorrsImmob noiseCorrsImmob(locomotionEnhancedPairs(w,1),locomotionEnhancedPairs(w,2))];
            
            locomotionEnhancedNoiseCorrsShuf=[locomotionEnhancedNoiseCorrsShuf noiseCorrsShuf(locomotionEnhancedPairs(w,1),locomotionEnhancedPairs(w,2))];
            locomotionEnhancedNoiseCorrsMoveShuf=[locomotionEnhancedNoiseCorrsMoveShuf noiseCorrsMoveShuf(locomotionEnhancedPairs(w,1),locomotionEnhancedPairs(w,2))];
            locomotionEnhancedNoiseCorrsImmobShuf=[locomotionEnhancedNoiseCorrsImmobShuf noiseCorrsImmobShuf(locomotionEnhancedPairs(w,1),locomotionEnhancedPairs(w,2))];
            
        end
        
        for w=1:size(locomotionSuppressedPairs,1)
            locomotionSuppressedNoiseCorrs=[locomotionSuppressedNoiseCorrs noiseCorrs(locomotionSuppressedPairs(w,1),locomotionSuppressedPairs(w,2))];
            locomotionSuppressedNoiseCorrsMove=[locomotionSuppressedNoiseCorrsMove noiseCorrsMove(locomotionSuppressedPairs(w,1),locomotionSuppressedPairs(w,2))];
            locomotionSuppressedNoiseCorrsImmob=[locomotionSuppressedNoiseCorrsImmob noiseCorrsImmob(locomotionSuppressedPairs(w,1),locomotionSuppressedPairs(w,2))];
            
            locomotionSuppressedNoiseCorrsShuf=[locomotionSuppressedNoiseCorrsShuf noiseCorrsShuf(locomotionSuppressedPairs(w,1),locomotionSuppressedPairs(w,2))];
            locomotionSuppressedNoiseCorrsMoveShuf=[locomotionSuppressedNoiseCorrsMoveShuf noiseCorrsMoveShuf(locomotionSuppressedPairs(w,1),locomotionSuppressedPairs(w,2))];
            locomotionSuppressedNoiseCorrsImmobShuf=[locomotionSuppressedNoiseCorrsImmobShuf noiseCorrsImmobShuf(locomotionSuppressedPairs(w,1),locomotionSuppressedPairs(w,2))];
            
        end
        
        for w=1:size(acrossClassPairs,1)
            acrossClassesNoiseCorrs=[acrossClassesNoiseCorrs noiseCorrs(acrossClassPairs(w,1),acrossClassPairs(w,2))];
            acrossClassesNoiseCorrsMove=[acrossClassesNoiseCorrsMove noiseCorrsMove(acrossClassPairs(w,1),acrossClassPairs(w,2))];
            acrossClassesNoiseCorrsImmob=[acrossClassesNoiseCorrsImmob noiseCorrsImmob(acrossClassPairs(w,1),acrossClassPairs(w,2))];
            
            acrossClassesNoiseCorrsShuf=[acrossClassesNoiseCorrsShuf noiseCorrsShuf(acrossClassPairs(w,1),acrossClassPairs(w,2))];
            acrossClassesNoiseCorrsMoveShuf=[acrossClassesNoiseCorrsMoveShuf noiseCorrsMoveShuf(acrossClassPairs(w,1),acrossClassPairs(w,2))];
            acrossClassesNoiseCorrsImmobShuf=[acrossClassesNoiseCorrsImmobShuf noiseCorrsImmobShuf(acrossClassPairs(w,1),acrossClassPairs(w,2))];
            
        end
        for w=1:size(locomotionUnmodPairs,1)
            locomotionUnmodNoiseCorrs=[locomotionUnmodNoiseCorrs noiseCorrs(locomotionUnmodPairs(w,1),locomotionUnmodPairs(w,2))];
        end
        
        alllocomotionEnhancedNoiseCorrs=[alllocomotionEnhancedNoiseCorrs locomotionEnhancedNoiseCorrs];
        alllocomotionEnhancedNoiseCorrsMove=[alllocomotionEnhancedNoiseCorrsMove locomotionEnhancedNoiseCorrsMove];
        alllocomotionEnhancedNoiseCorrsImmob=[alllocomotionEnhancedNoiseCorrsImmob locomotionEnhancedNoiseCorrsImmob];
        
        alllocomotionSuppressedNoiseCorrs=[alllocomotionSuppressedNoiseCorrs locomotionSuppressedNoiseCorrs];
        alllocomotionSuppressedNoiseCorrsMove=[alllocomotionSuppressedNoiseCorrsMove locomotionSuppressedNoiseCorrsMove];
        alllocomotionSuppressedNoiseCorrsImmob=[alllocomotionSuppressedNoiseCorrsImmob locomotionSuppressedNoiseCorrsImmob];
        
        allacrossClassesNoiseCorrs=[allacrossClassesNoiseCorrs acrossClassesNoiseCorrs];
        allacrossClassesNoiseCorrsMove=[allacrossClassesNoiseCorrsMove acrossClassesNoiseCorrsMove];
        allacrossClassesNoiseCorrsImmob=[allacrossClassesNoiseCorrsImmob acrossClassesNoiseCorrsImmob];
        
        alllocomotionUnmodNoiseCorrs=[alllocomotionUnmodNoiseCorrs locomotionUnmodNoiseCorrs];
        
        alllocomotionEnhancedNoiseCorrsShuf=[alllocomotionEnhancedNoiseCorrsShuf locomotionEnhancedNoiseCorrsShuf];
        alllocomotionEnhancedNoiseCorrsMoveShuf=[alllocomotionEnhancedNoiseCorrsMoveShuf locomotionEnhancedNoiseCorrsMoveShuf];
        alllocomotionEnhancedNoiseCorrsImmobShuf=[alllocomotionEnhancedNoiseCorrsImmobShuf locomotionEnhancedNoiseCorrsImmobShuf];
        
        alllocomotionSuppressedNoiseCorrsShuf=[alllocomotionSuppressedNoiseCorrsShuf locomotionSuppressedNoiseCorrsShuf];
        alllocomotionSuppressedNoiseCorrsMoveShuf=[alllocomotionSuppressedNoiseCorrsMoveShuf locomotionSuppressedNoiseCorrsMoveShuf];
        alllocomotionSuppressedNoiseCorrsImmobShuf=[alllocomotionSuppressedNoiseCorrsImmobShuf locomotionSuppressedNoiseCorrsImmobShuf];
        
        allacrossClassesNoiseCorrsShuf=[allacrossClassesNoiseCorrsShuf acrossClassesNoiseCorrsShuf];
        allacrossClassesNoiseCorrsMoveShuf=[allacrossClassesNoiseCorrsMoveShuf acrossClassesNoiseCorrsMoveShuf];
        allacrossClassesNoiseCorrsImmobShuf=[allacrossClassesNoiseCorrsImmobShuf acrossClassesNoiseCorrsImmobShuf];
        
        % example plotting pairs
        plotexample=0;
        if plotexample
            figure;
            subplot(1,3,1)
            plot(ensembleNormSoundEvokedResps(3,:),ensembleNormSoundEvokedResps(20,:),'ko');
            [m n]=polyfit(ensembleNormSoundEvokedResps(3,:),ensembleNormSoundEvokedResps(20,:),1);
            hold on
            plot([-1 3],polyval(m,[-1 3]),'r--')
            axis([-1 3 -1 3])
            title(['Cell 3 and 20, R= ' num2str(noiseCorrs(3,20))])
            
            subplot(1,3,2)
            plot(ensembleNormSoundEvokedResps(1,:),ensembleNormSoundEvokedResps(22,:),'ko');
            [m n]=polyfit(ensembleNormSoundEvokedResps(1,:),ensembleNormSoundEvokedResps(22,:),1);
            hold on
            plot([-1 3],polyval(m,[-1 3]),'r--')
            axis([-.5 .5 -.5 .5])
            title(['Cell 1 and 22, R= ' num2str(noiseCorrs(1,22))])
            
            subplot(1,3,3)
            plot(ensembleNormSoundEvokedResps(3,:),ensembleNormSoundEvokedResps(22,:),'ko');
            [m n]=polyfit(ensembleNormSoundEvokedResps(3,:),ensembleNormSoundEvokedResps(22,:),1);
            hold on
            plot([-1 3],polyval(m,[-1 3]),'r--')
            axis([-.5 3 -.5 3])
            
            title(['Cell 3 and 22, R= ' num2str(noiseCorrs(3,22))])
        end
        
        
        
        
    end
    
    
    
    
    allCellResponsiveness=[allCellResponsiveness;cellResponsiveness];
    allCellImmobileResponsiveness=[allCellImmobileResponsiveness;cellResponsivenessImmobile];
    
    
    % calculating the correlation between each session-concatenated dff and
    % session-concatenated speed
    if ~isempty(acrossSessionsSpeeds)
        acrossExperimentsSpeeds=[acrossExperimentsSpeeds acrossSessionsSpeeds];
        
        acrossSessionsSpeeds=smooth(acrossSessionsSpeeds,6);
        
        for kp=1:size(alldff,1)
            [r p]=corrcoef(alldff(kp,:),acrossSessionsSpeeds);
            epochdffspeedcorrs=[epochdffspeedcorrs r(1,2)];
            epochdffspeedcorrsSignificance=[epochdffspeedcorrsSignificance p(1,2)];
            
            [r2 p2]=corrcoef(alldff(kp,acrossSessionsSpeeds>0),acrossSessionsSpeeds(acrossSessionsSpeeds>0));
            epochdffspeedcorrsNo0=[epochdffspeedcorrsNo0 r2(1,2)];
            epochdffspeedcorrsSignificanceNo0=[epochdffspeedcorrsSignificanceNo0 p2(1,2)];
            
            
        end
        alldffspeedcorrs=[alldffspeedcorrs epochdffspeedcorrs];
        alldffspeedcorrsSignificance=[alldffspeedcorrsSignificance epochdffspeedcorrsSignificance];
        alldffspeedcorrsRanges=[alldffspeedcorrsRanges range(epochdffspeedcorrs)];
        
        alldffspeedcorrsNo0=[alldffspeedcorrsNo0 epochdffspeedcorrsNo0];
        alldffspeedcorrsSignificanceNo0=[alldffspeedcorrsSignificanceNo0 epochdffspeedcorrsSignificanceNo0];
        
        %  SPEED PREDICTION
        dologspeeds=1;
     
        numrepeats=200;
        if dologspeeds
            [realR, realRNo0, shufR, shufRNo0]=decodeSpeedFromDFF(alldff,log(acrossSessionsSpeeds+1),numrepeats,i);
        else
            [realR, realRNo0, shufR, shufRNo0]=decodeSpeedFromDFF(alldff,acrossSessionsSpeeds,numrepeats,i);
            
        end
        participatingEnsemblesSpeedPrediction=[participatingEnsemblesSpeedPrediction 1];
        allSpeedPredictionRs=[allSpeedPredictionRs realR];
        allSpeedPredictionRsNo0=[allSpeedPredictionRsNo0 realRNo0];
        allSpeedPredictionRsSHUF=[allSpeedPredictionRsSHUF shufR];
        allSpeedPredictionRsNo0SHUF=[allSpeedPredictionRsNo0SHUF shufRNo0];
        predictionEnsembleSize=[predictionEnsembleSize size(alldff,1)];
        
        
        % comparing NONE to CONST
        curnumcells=size(alldffNONE,1);
        if ~isempty(alldffNONE) & ~isempty(alldffCONST) & sum(allmovNONE)>0
            allCellImmobileResponsivenessCONST=[allCellImmobileResponsivenessCONST; cellResponsivenessImmobile];
            figure;
            subplot(1,2,1);
            plot(find(allmovNONE),0.5,'ko')
            h=area(allmovNONE*curnumcells*1.1);
            h.FaceColor = [1 1 1]*0.8;
            h.LineStyle= 'none';
            hold on
            for pp=1:curnumcells,plot(alldffNONE(pp,:)+pp);hold on;end
            title('No background sound')
            subplot(1,2,2);
            h=area(allmovCONST*curnumcells*1.1);
            h.FaceColor = [1 1 1]*0.8;
            h.LineStyle= 'none';
            hold on
            plot(find(allmovCONST),0.5,'ko')
            
            for pp=1:curnumcells,plot(alldffCONST(pp,:)+pp);hold on;end
            title('Background noise')
            
            for kp=1:size(alldffNONE,1)
                [rNONE pNONE]=corrcoef(alldffNONE(kp,:),acrossSessionsSpeedsNONE);
                epochdffspeedcorrsNONE=[epochdffspeedcorrsNONE rNONE(1,2)];
                meanDuringImmobileNONE=mean(alldffNONE(kp,acrossSessionsSpeedsNONE==0));
                meanDuringMoveNONE=mean(alldffNONE(kp,acrossSessionsSpeedsNONE>0));
                allmeanDuringImmobileNONE=[allmeanDuringImmobileNONE meanDuringImmobileNONE];
                allmeanDuringMoveNONE=[allmeanDuringMoveNONE meanDuringMoveNONE];
                
                [rCONST pCONST]=corrcoef(alldffCONST(kp,:),acrossSessionsSpeedsCONST);
                if rCONST>0.4 | rCONST <-0.2
                    % keyboard
                end
                epochdffspeedcorrsCONST=[epochdffspeedcorrsCONST rCONST(1,2)];
                meanDuringImmobileCONST=mean(alldffCONST(kp,acrossSessionsSpeedsCONST==0));
                meanDuringMoveCONST=mean(alldffCONST(kp,acrossSessionsSpeedsCONST>0));
                allmeanDuringImmobileCONST=[allmeanDuringImmobileCONST meanDuringImmobileCONST];
                allmeanDuringMoveCONST=[allmeanDuringMoveCONST meanDuringMoveCONST];
                
            end
            figure;plot(epochdffspeedcorrsNONE,epochdffspeedcorrsCONST,'ko');hold on;plot([-1 1],[-1 1]);xlabel('No background sound');ylabel('Background sound')
            alldffspeedcorrsCONST=[alldffspeedcorrsCONST epochdffspeedcorrsCONST];
            alldffspeedcorrsNONE=[alldffspeedcorrsNONE epochdffspeedcorrsNONE];
            
        end
    else
        participatingEnsemblesSpeedPrediction=[participatingEnsemblesSpeedPrediction 0];
        
    end
    
    

    curCellCenters=allresults{i}{1}.cellCenters;
    alldists=dist(curCellCenters');
    
    alldistsunique=alldists(tril(ones(numcells),-1)==1);
    meanDistAllPairs=mean(alldistsunique);
    positivelyLocomotionModulatedCells=find(difSpontLocomotionVsImmobility==1);
    
    nonPositivelyLocomotionModulatedCells=find(difSpontLocomotionVsImmobility<1);
    relativeDistPos=nan;
    relativeDistNeg=nan;
    if length(nonPositivelyLocomotionModulatedCells)>1
        
        negativedists=dist(curCellCenters(nonPositivelyLocomotionModulatedCells,:)');
        negativedistsunique=negativedists(tril(ones(length(nonPositivelyLocomotionModulatedCells)),-1)==1);
        relativeDistNeg=mean(negativedistsunique)/mean(alldistsunique);
        allrelativeDistNeg=[allrelativeDistNeg; relativeDistNeg];
        
        
        
    end
    if length(positivelyLocomotionModulatedCells)>1
        positivedists=dist(curCellCenters(positivelyLocomotionModulatedCells,:)');
        positivedistsunique=positivedists(tril(ones(length(positivelyLocomotionModulatedCells)),-1)==1);
        relativeDistPos=mean(positivedistsunique)/mean(alldistsunique);
        allrelativeDistPos=[allrelativeDistPos; relativeDistPos];
        
    end
    
    allProportionOfPosSpeedModCellsInEnsemble=[allProportionOfPosSpeedModCellsInEnsemble mean(difSpontLocomotionVsImmobility==1)];
    
    
    
 
    
    prestimwin_ton3=floor((size(TON3resps,4)/2)-stimwinlength:(size(TON3resps,4)/2)-1);
    stimwin_ton3=floor((size(TON3resps,4)/2)+1:(size(TON3resps,4)/2)+stimwinlength);
    numTONE3stims=size(TON3respsMOV,2);
    numcells3=size(TON3respsMOV,1);
    for c=1:numcells3
        
        for ss=1:numTONE3stims
            % LOCOMOTION RESPONSES
            curCellMOV=squeeze(TON3respsMOV(c,ss,:,:));
            notAllNANrowsMOV=find(sum(~isnan(curCellMOV)')~=0);
            curCellMOVnoNAN=curCellMOV(notAllNANrowsMOV,:);
            
            % IMMOBILITY RESPONSES
            curCellImmobile=squeeze(TON3respsImmobile(c,ss,:,:));
            notAllNANrowsImmobile=find(sum(~isnan(curCellImmobile)')~=0);
            curCellImmobilenoNAN=curCellImmobile(notAllNANrowsImmobile,:);
            
            
            % EVOKED RESPONSES DURING LOCOMOTION
            curCellMove_evoked =nanmean(curCellMOVnoNAN(:,stimwin_ton3),2);
            
            % EVOKED RESPONSES DURING IMMOBILITY
            curCellImmobile_evoked=nanmean(curCellImmobilenoNAN(:,stimwin_ton3),2);
            
            % PRE-STIM ACTIVITY DURING LOCOMOTION
            curCellMove_pre =nanmean(curCellMOVnoNAN(:,prestimwin_ton3),2);
            
            % PRE-STIM ACTIVITY DURING IMMOBILITY
            curCellImmobile_pre=nanmean(curCellImmobilenoNAN(:,prestimwin_ton3),2);
            
            % EVOKED MINUS PRE-STIM DURING LOCOMOTION
            evokedMinusPre_Move=curCellMove_evoked-curCellMove_pre;
            
            % EVOKED MINUS PRE-STIM DURING IMMOBILITY
            evokedMinusPre_Immobile=curCellImmobile_evoked-curCellImmobile_pre;
            
            % responsiveness to this stim in immobility
            [curCellResponsiveness3, curCellResponsiveness3P]=ttest([curCellImmobile_pre],[curCellImmobile_evoked],'Tail','left');
            
       
            if curCellResponsiveness3==1
                % DIFFERENCE IN EVOKED MINUS PRE DURING LOCOMOTION VS IMMOBILITY
                % significantly higher during locomotion = 1, immobility = -1, no difference = 0
                [h, p]=ttest2(evokedMinusPre_Move,evokedMinusPre_Immobile);
                sigModulation=0;
                if p<0.05
                    if mean(evokedMinusPre_Move)>mean(evokedMinusPre_Immobile)
                        sigModulation=1;
                    else
                        sigModulation=-1;
                    end
                    
                end
                
                difEvokedMinusPreLocomotionVsImmobilityTON3=[difEvokedMinusPreLocomotionVsImmobilityTON3;sigModulation];
                allmeanEvokedMinusPreLocomotionTON3=[allmeanEvokedMinusPreLocomotionTON3;mean(evokedMinusPre_Move)];
                allmeanEvokedMinusPreImmobilityTON3=[allmeanEvokedMinusPreImmobilityTON3;mean(evokedMinusPre_Immobile)];
                toplotton3=0;
                if toplotton3
                    figure('Position', [50, 50, 300, 800])
                    subplot(3,1,1);
                    errorbar(nanmean(curCellImmobilenoNAN),nanstd(curCellImmobilenoNAN)/sqrt(size(curCellImmobilenoNAN,1)),'r', 'Linewidth', 1);
                    hold on
                    errorbar(nanmean(curCellMOVnoNAN),nanstd(curCellMOVnoNAN)/sqrt(size(curCellMOVnoNAN,1)),'g', 'Linewidth', 1);
                    xlim([1 12])
                    title([allresults{1,i}{1,1}.batch ' cell ' num2str(c) ' stim ' num2str(ss) ' mod= ' num2str(sigModulation)])
                    y=ylim;
                    ylim([y(1) max(y(2),.5)])
                    
                    subplot(3,1,2);
                    imagesc(curCellImmobilenoNAN);
                    xlim([1 12])
                    title('Immobile')
                    c1=caxis;
                    
                    subplot(3,1,3);
                    imagesc(curCellMOVnoNAN);
                    xlim([1 12])
                    title('Locomotion')
                    c2=caxis;
                    
                    subplot(3,1,2), caxis([min(c1(1),c2(1)) max(c1(2),c2(2))])
                    subplot(3,1,3), caxis([min(c1(1),c2(1)) max(c1(2),c2(2))])
                    
                    saveas(gcf,[ton3folder  allresults{1,i}{1,1}.batch ' cell ' num2str(c) ' stim ' num2str(ss)],'jpg');
                    saveas(gcf,[ton3folder  allresults{1,i}{1,1}.batch ' cell ' num2str(c) ' stim ' num2str(ss)],'fig');
                    print([ton3folder allresults{1,i}{1,1}.batch ' cell ' num2str(c) ' stim ' num2str(ss)],'-dpdf','-bestfit');
                    close all
                end
                
            end
        end
        
    end
    
    
    

    
    prestimwin_NAT3=floor((size(NAT3resps,4)/2)-stimwinlength:(size(NAT3resps,4)/2)-1);
    stimwin_NAT3=floor((size(NAT3resps,4)/2)+1:(size(NAT3resps,4)/2)+stimwinlength);
    numNAT3stims=size(NAT3respsMOV,2);
    numcells3=size(NAT3respsMOV,1);
    for c=1:numcells3
        
        for ss=1:numNAT3stims
            % LOCOMOTION RESPONSES
            curCellMOV=squeeze(NAT3respsMOV(c,ss,:,:));
            notAllNANrowsMOV=find(sum(~isnan(curCellMOV)')~=0);
            curCellMOVnoNAN=curCellMOV(notAllNANrowsMOV,:);
            
            % IMMOBILITY RESPONSES
            curCellImmobile=squeeze(NAT3respsImmobile(c,ss,:,:));
            notAllNANrowsImmobile=find(sum(~isnan(curCellImmobile)')~=0);
            curCellImmobilenoNAN=curCellImmobile(notAllNANrowsImmobile,:);
            
            
            % EVOKED RESPONSES DURING LOCOMOTION
            curCellMove_evoked =nanmean(curCellMOVnoNAN(:,stimwin_NAT3),2);
            
            % EVOKED RESPONSES DURING IMMOBILITY
            curCellImmobile_evoked=nanmean(curCellImmobilenoNAN(:,stimwin_NAT3),2);
            
            % PRE-STIM ACTIVITY DURING LOCOMOTION
            curCellMove_pre =nanmean(curCellMOVnoNAN(:,prestimwin_NAT3),2);
            
            % PRE-STIM ACTIVITY DURING IMMOBILITY
            curCellImmobile_pre=nanmean(curCellImmobilenoNAN(:,prestimwin_NAT3),2);
            
            % EVOKED MINUS PRE-STIM DURING LOCOMOTION
            evokedMinusPre_Move=curCellMove_evoked-curCellMove_pre;
            
            % EVOKED MINUS PRE-STIM DURING IMMOBILITY
            evokedMinusPre_Immobile=curCellImmobile_evoked-curCellImmobile_pre;
            
            % responsiveness to this stim across states
            [curCellResponsiveness3, curCellResponsiveness3P]=ttest([curCellImmobile_pre; curCellMove_pre],[curCellImmobile_evoked;curCellMove_evoked],'Tail','left');
            %   cellResponsiveness=[cellResponsiveness;curCellResponsivenessH];
            
            if curCellResponsiveness3==1
                % DIFFERENCE IN EVOKED MINUS PRE DURING LOCOMOTION VS IMMOBILITY
                % significantly higher during locomotion = 1, immobility = -1, no difference = 0
                [h, p]=ttest2(evokedMinusPre_Move,evokedMinusPre_Immobile);
                sigModulation=0;
                if p<0.05
                    if mean(evokedMinusPre_Move)>mean(evokedMinusPre_Immobile)
                        sigModulation=1;
                    else
                        sigModulation=-1;
                    end
                    
                end
                
                difEvokedMinusPreLocomotionVsImmobilityNAT3=[difEvokedMinusPreLocomotionVsImmobilityNAT3;sigModulation];
                allmeanEvokedMinusPreLocomotionNAT3=[allmeanEvokedMinusPreLocomotionNAT3;mean(evokedMinusPre_Move)];
                allmeanEvokedMinusPreImmobilityNAT3=[allmeanEvokedMinusPreImmobilityNAT3;mean(evokedMinusPre_Immobile)];
                
                
                
                % difEvokedMinusPreLocomotionVsImmobilityAllExps=[difEvokedMinusPreLocomotionVsImmobilityAllExps;sigModulation];
                
            end
        end
        
    end
    
    

    minnumtrialsNAT=8;
    if size(NATresps,2)==5
        NATresps(:,2,:,:)=[];
        NATrespsMOV(:,2,:,:)=[];
        NATrespsImmobile(:,2,:,:)=[];
    end
    
    prestimwin_nat=floor((size(NATresps,4)/2)-stimwinlength:(size(NATresps,4)/2)-1);
    stimwin_nat=floor((size(NATresps,4)/2)+1:(size(NATresps,4)/2)+stimwinlength);
    numstims_nat=size(NATresps,2);
    numcells_nat=size(NATresps,1);
    
    
    for c=1:numcells_nat
        
        for ss=1:numstims_nat
            % LOCOMOTION RESPONSES
            curCellMOV=squeeze(NATrespsMOV(c,ss,:,:));
            notAllNANrowsMOV=find(sum(~isnan(curCellMOV)')~=0);
            curCellMOVnoNAN=curCellMOV(notAllNANrowsMOV,:);
            
            % IMMOBILITY RESPONSES
            curCellImmobile=squeeze(NATrespsImmobile(c,ss,:,:));
            notAllNANrowsImmobile=find(sum(~isnan(curCellImmobile)')~=0);
            curCellImmobilenoNAN=curCellImmobile(notAllNANrowsImmobile,:);
            
            numMoveTrials=size(curCellMOVnoNAN,1);
            numImmobileTrials=size(curCellImmobilenoNAN,1);
            % only gathering if there are enough trials in each state
            if numMoveTrials>minnumtrialsNAT & numImmobileTrials>minnumtrialsNAT
                
                % EVOKED RESPONSES DURING LOCOMOTION
                curCellMove_evoked =nanmean(curCellMOVnoNAN(:,stimwin_nat),2);
                
                % EVOKED RESPONSES DURING IMMOBILITY
                curCellImmobile_evoked=nanmean(curCellImmobilenoNAN(:,stimwin_nat),2);
                
                % PRE-STIM ACTIVITY DURING LOCOMOTION
                curCellMove_pre =nanmean(curCellMOVnoNAN(:,prestimwin_nat),2);
                
                % PRE-STIM ACTIVITY DURING IMMOBILITY
                curCellImmobile_pre=nanmean(curCellImmobilenoNAN(:,prestimwin_nat),2);
                
                % EVOKED MINUS PRE-STIM DURING LOCOMOTION
                evokedMinusPre_Move=curCellMove_evoked-curCellMove_pre;
                
                % EVOKED MINUS PRE-STIM DURING IMMOBILITY
                evokedMinusPre_Immobile=curCellImmobile_evoked-curCellImmobile_pre;
                
                % responsiveness to this stim across states
                [curCellResponsiveness3, curCellResponsiveness3P]=ttest([curCellImmobile_pre; curCellMove_pre],[curCellImmobile_evoked;curCellMove_evoked],'Tail','left');
                
                % responsiveness to this stim in immobility
                
                [curCellResponsiveness3IMMOB, curCellResponsiveness3PIMMOB]=ttest([curCellImmobile_pre],[curCellImmobile_evoked],'Tail','left');
                
                
                if curCellResponsiveness3IMMOB==1
                    % DIFFERENCE IN EVOKED MINUS PRE DURING LOCOMOTION VS IMMOBILITY
                    % significantly higher during locomotion = 1, immobility = -1, no difference = 0
                    [h, p]=ttest2(evokedMinusPre_Move,evokedMinusPre_Immobile);
                    sigModulation=0;
                    if p<0.05
                        if mean(evokedMinusPre_Move)>mean(evokedMinusPre_Immobile)
                            sigModulation=1;
                        else
                            sigModulation=-1;
                        end
                        
                    end
                    
                    difEvokedMinusPreLocomotionVsImmobilityNAT=[difEvokedMinusPreLocomotionVsImmobilityNAT;sigModulation];
                    allmeanEvokedMinusPreLocomotionNAT=[allmeanEvokedMinusPreLocomotionNAT;mean(evokedMinusPre_Move)];
                    allmeanEvokedMinusPreImmobilityNAT=[allmeanEvokedMinusPreImmobilityNAT;mean(evokedMinusPre_Immobile)];
                    
                    
                end
            end
            
        end
        
    end
    
    
    
    numcells_nat=size(NATrespsMOV,1);
    prestimwin_nat=floor((size(NATresps,4)/2)-stimwinlength:(size(NATresps,4)/2)-1);
    stimwin_nat=floor((size(NATresps,4)/2)+1:(size(NATresps,4)/2)+stimwinlength);
    allCellResponsiveness_nat=[];
    allCellResponsiveness_perstim_ensemble=[];
    % DETERMINE SOUND RESPONSIVENESS
    for c=1:numcells_nat
        cellResponsiveness_nat=[];
        cellresponsive_nat_anystim=[];
        if size(NATresps,2)==5
            NATresps(:,2,:,:)=[];
        end
        acrossStimPre=[];
        acrossStimDur=[];
        for s=1:size(NATresps,2)
            curCellEvoked_natstim=nanmean(squeeze(NATresps(c,s,:,stimwin_nat))');
            curCellPreEvoked_natstim=nanmean(squeeze(NATresps(c,s,:,prestimwin_nat))');
            
            goodinds_nat=find(~isnan(curCellPreEvoked_natstim)&~isnan(curCellEvoked_natstim));
            
            curCellEvoked_natstim=curCellEvoked_natstim(goodinds_nat);
            curCellPreEvoked_natstim=curCellPreEvoked_natstim(goodinds_nat);
            acrossStimPre=[acrossStimPre curCellPreEvoked_natstim];
            acrossStimDur=[acrossStimDur curCellEvoked_natstim];
            
            [h, p]=ttest(curCellPreEvoked_natstim,curCellEvoked_natstim,'Tail','left');
            cellResponsiveness_nat=[cellResponsiveness_nat h];
        end
        [hAcrossStim, pAcrossStim]=ttest(acrossStimPre,acrossStimDur,'Tail','left');
        
        acrossStimResponsiveness=[acrossStimResponsiveness hAcrossStim];
        if any(cellResponsiveness_nat(:)==1)
            cellresponsive_nat_anystim = 1;
        else
            cellresponsive_nat_anystim=0;
        end
        allCellResponsiveness_anystim=[allCellResponsiveness_anystim;cellresponsive_nat_anystim];
        allCellResponsiveness_perstim_ensemble=[allCellResponsiveness_perstim_ensemble; cellResponsiveness_nat];
        
        plotcellsresponsivenessNAT=0;
        if plotcellsresponsivenessNAT
            
            figure('position',[30 30 200 700]);
            numstim1=size(NATresps,2);
            for s=1:numstim1
                curCellEvoked_natstim1=nanmean(squeeze(NATresps(c,s,:,stimwin_nat))');
                curCellPreEvoked_natstim1=nanmean(squeeze(NATresps(c,s,:,prestimwin_nat))');
                
                goodinds_nat=find(~isnan(curCellPreEvoked_natstim1)&~isnan(curCellEvoked_natstim1));
                
                curCellEvoked_natstim1=curCellEvoked_natstim1(goodinds_nat);
                curCellPreEvoked_natstim1=curCellPreEvoked_natstim1(goodinds_nat);
                [h1, p1]=ttest(curCellPreEvoked_natstim1,curCellEvoked_natstim1,'Tail','left');
                subplot(numstim1,1,s)
                plot(squeeze(NATresps(c,s,:,:))','k');
                hold on;
                plot(nanmean(squeeze(NATresps(c,s,:,:))),'r','linewidth',2);
                
                title(['P= ' num2str(p1)])
                
            end
            
            
            close all
        end
        
        
        
    end
    allCellResponsiveness_perstim=[allCellResponsiveness_perstim;allCellResponsiveness_perstim_ensemble];
    
    
    % Determine stimulus-specific locomotion modulation
    allcell_stimdFFmodulation_nat_ensemble=[];
    allcell_stimSNRmodulation_nat_ensemble=[];
    minnumtrialsNAT=4;
    for c=1:numcells_nat
        cell_stimdFFmodulation_nat=[];
        cell_stimSNRmodulation_nat=[];
        allmeanevokedLocomotion_natcellstim=[];
        allmeanevokedImmobility_natcellstim =[];
        allmeanSNRLocomotion_natcellstim=[];
        allmeanSNRImmobility_natcellstim=[];
        if size(NATrespsMOV,2)==5
            NATrespsMOV(:,2,:,:)=[];
            NATrespsImmobile(:,2,:,:)=[];
        end
        %For each stimulus within this cell
        acrossStimSNRsMove=[];
        acrossStimSNRsImmobile=[];
        
        for s=1:size(NATrespsMOV,2)
            %take data from Mov and Immobile vector
            curCell_stimEvokedMove=squeeze(NATrespsMOV(c,s,:,:));
            curCell_stimEvokedImmobile=squeeze(NATrespsImmobile(c,s,:,:));
            %Evoked Move & Immobile analysis get movement from both stim
            %windows of interest, then get SNR (Evoked-pre), subtract out NaNs
            %from all three at the end
            curCell_stimEvokedMove_during =nanmean(curCell_stimEvokedMove(:,stimwin)')';
            curCell_stimEvokedMove_pre =nanmean(curCell_stimEvokedMove(:,prestimwin)')';
            
            curCell_stimEvokedMove_SNR =(curCell_stimEvokedMove_during -curCell_stimEvokedMove_pre);
            
            curCell_stimEvokedMove_pre =curCell_stimEvokedMove_pre(~isnan(curCell_stimEvokedMove_pre));
            curCell_stimEvokedMove_during =curCell_stimEvokedMove_during(~isnan(curCell_stimEvokedMove_during));
            curCell_stimEvokedMove_SNR = curCell_stimEvokedMove_SNR(~isnan(curCell_stimEvokedMove_SNR));
            
            curCell_stimEvokedImmobile_during=nanmean(curCell_stimEvokedImmobile(:,stimwin)')';
            curCell_stimEvokedImmobile_pre =nanmean(curCell_stimEvokedImmobile(:,prestimwin)')';
            
            curCell_stimEvokedImmobile_SNR =(curCell_stimEvokedImmobile_during -curCell_stimEvokedImmobile_pre);
            
            curCell_stimEvokedImmobile_pre =curCell_stimEvokedImmobile_pre(~isnan(curCell_stimEvokedImmobile_pre));
            curCell_stimEvokedImmobile_during=curCell_stimEvokedImmobile_during(~isnan(curCell_stimEvokedImmobile_during));
            curCell_stimEvokedImmobile_SNR=curCell_stimEvokedImmobile_SNR(~isnan(curCell_stimEvokedImmobile_SNR));
            
            
            
            acrossStimSNRsMove=[acrossStimSNRsMove; curCell_stimEvokedMove_SNR];
            acrossStimSNRsImmobile=[acrossStimSNRsImmobile; curCell_stimEvokedImmobile_SNR];
            %Test difference in dff between locomotion and immobility during
            %evoked window
            [h11, p11]=ttest2(curCell_stimEvokedMove_during,curCell_stimEvokedImmobile_during);
            
            dFFsigModulation=0;
            if p11<0.05
                if mean(curCell_stimEvokedMove_during)>mean(curCell_stimEvokedImmobile_during)
                    dFFsigModulation=1;
                else
                    dFFsigModulation=-1;
                end
            end
            
            
            [h22, p22]=ttest2(curCell_stimEvokedMove_SNR,curCell_stimEvokedImmobile_SNR);
            
            SNRsigModulation=0;
            if p22<0.05
                if mean(curCell_stimEvokedMove_SNR)>mean(curCell_stimEvokedImmobile_SNR)
                    SNRsigModulation=1;
                else
                    SNRsigModulation=-1;
                end
            end
            
            numMoveTrials=length(curCell_stimEvokedMove_during);
            numImmobileTrials=length(curCell_stimEvokedImmobile_during);
            if numMoveTrials>minnumtrialsNAT & numImmobileTrials>minnumtrialsNAT
              
                cell_stimdFFmodulation_nat = [cell_stimdFFmodulation_nat dFFsigModulation];
                
                cell_stimSNRmodulation_nat = [cell_stimSNRmodulation_nat SNRsigModulation];
                allmeanevokedLocomotion_natcellstim=[allmeanevokedLocomotion_natcellstim mean(curCell_stimEvokedMove_during)];
                allmeanevokedImmobility_natcellstim=[allmeanevokedImmobility_natcellstim mean(curCell_stimEvokedImmobile_during)];
                allmeanSNRLocomotion_natcellstim=[allmeanSNRLocomotion_natcellstim mean(curCell_stimEvokedMove_SNR)];
                allmeanSNRImmobility_natcellstim=[allmeanSNRImmobility_natcellstim mean(curCell_stimEvokedImmobile_SNR)];
            else
                cell_stimdFFmodulation_nat = [cell_stimdFFmodulation_nat nan];
                
                cell_stimSNRmodulation_nat = [cell_stimSNRmodulation_nat nan];
                allmeanevokedLocomotion_natcellstim=[allmeanevokedLocomotion_natcellstim nan];
                allmeanevokedImmobility_natcellstim=[allmeanevokedImmobility_natcellstim nan];
                allmeanSNRLocomotion_natcellstim=[allmeanSNRLocomotion_natcellstim nan];
                allmeanSNRImmobility_natcellstim=[allmeanSNRImmobility_natcellstim nan];
                
            end
        end
        
        
        
        allcell_stimdFFmodulation_nat_ensemble=[allcell_stimdFFmodulation_nat_ensemble;cell_stimdFFmodulation_nat];
        allcell_stimSNRmodulation_nat_ensemble=[allcell_stimSNRmodulation_nat_ensemble;cell_stimSNRmodulation_nat];
      
        ensemble_ID = [ensemble_ID; i];
        ensembleCELL_ID{end+1}=[allresults{1,i}{1,1}.batch '-' num2str(c)];
        
        allcell_stimdFFmodulation_nat=[allcell_stimdFFmodulation_nat;cell_stimdFFmodulation_nat];
        allcell_stimSNRmodulation_nat=[allcell_stimSNRmodulation_nat;cell_stimSNRmodulation_nat];
        cell_specificmeanevokedLocomotion_natcellstim =[cell_specificmeanevokedLocomotion_natcellstim; allmeanevokedLocomotion_natcellstim];
        cell_specificmeanevokedImmobility_natcellstim =[cell_specificmeanevokedImmobility_natcellstim; allmeanevokedImmobility_natcellstim];
        cell_specificmeanSNRLocomotion_natcellstim =[cell_specificmeanSNRLocomotion_natcellstim; allmeanSNRLocomotion_natcellstim];
        cell_specificmeanSNRImmobility_natcellstim =[cell_specificmeanSNRImmobility_natcellstim; allmeanSNRImmobility_natcellstim];
        
      
        [hAcrossStimSNR, pAcrossStimSNR]=ttest2(acrossStimSNRsMove,acrossStimSNRsImmobile);
        acrossStimSNRsigModulation=0;
        if pAcrossStimSNR<0.05
            if mean(acrossStimSNRsMove)>mean(acrossStimSNRsImmobile)
                acrossStimSNRsigModulation=1;
            else
                acrossStimSNRsigModulation=-1;
            end
        end
        allcell_acrossStimSNRmodulation_nat = [allcell_acrossStimSNRmodulation_nat acrossStimSNRsigModulation];
        
  
    end
    dodecoding=1;
    if dodecoding
        rng(77)
        numshufs1=200;
      
        respAmpsMovePre=[];
        respAmpsMoveDuring=[];
        respAmpsImmobilePre=[];
        respAmpsImmobileDuring=[];
        
        for c=1:numcells
            
            curCellMove=squeeze(BBNrespsMOV(c,1,:,:));
            curCellImmobile=squeeze(BBNrespsImmobile(c,1,:,:));
            
            curCellMove_pre=nanmean(curCellMove(:,prestimwin)')';
            curCellMove_pre=curCellMove_pre(~isnan(curCellMove_pre));
            
            curCellImmobile_pre=nanmean(curCellImmobile(:,prestimwin)')';
            curCellImmobile_pre=curCellImmobile_pre(~isnan(curCellImmobile_pre));
            
            curCellMove_during =nanmean(curCellMove(:,stimwin)')';
            curCellMove_during =curCellMove_during(~isnan(curCellMove_during));
            
            curCellImmobile_during=nanmean(curCellImmobile(:,stimwin)')';
            curCellImmobile_during=curCellImmobile_during(~isnan(curCellImmobile_during));
            
            respAmpsMovePre=[respAmpsMovePre curCellMove_pre ];
            respAmpsMoveDuring=[respAmpsMoveDuring curCellMove_during];
            respAmpsImmobilePre=[respAmpsImmobilePre curCellImmobile_pre];
            respAmpsImmobileDuring=[respAmpsImmobileDuring curCellImmobile_during];
            
        end
        
        
        minnumtrials=min(size(respAmpsMovePre,1),size(respAmpsImmobilePre,1));
        
        if minnumtrials>12&numcells>10
            allminnumtrials=[allminnumtrials minnumtrials];
            allnumcellsperensembledecoding=[allnumcellsperensembledecoding;numcells];
            participatingEnsemblesDecoding=[participatingEnsemblesDecoding 1];
            respAmpsMatImmobile=([respAmpsImmobilePre(1:minnumtrials,:); respAmpsImmobileDuring(1:minnumtrials,:)]);
            respAmpsMatMove=([respAmpsMovePre(1:minnumtrials,:); respAmpsMoveDuring(1:minnumtrials,:)]);
            respAmpsMatState=([respAmpsImmobileDuring(1:minnumtrials,:); respAmpsMoveDuring(1:minnumtrials,:)]);
            respAmpsALLMat=[respAmpsImmobilePre(1:minnumtrials,:); respAmpsImmobileDuring(1:minnumtrials,:);respAmpsMovePre(1:minnumtrials,:); respAmpsMoveDuring(1:minnumtrials,:)];
            tagsAll=[zeros(1,minnumtrials) ones(1,minnumtrials) 2*ones(1,minnumtrials) 3*ones(1,minnumtrials)];
            
            
            
            %how well can network detect a stim in immobility
            tags1=[zeros(1,minnumtrials) ones(1,minnumtrials)];
            dosinglecells=0;
            if dosinglecells
            for c1=1:numcells
                preImmobVar=[preImmobVar var(respAmpsImmobilePre(1:minnumtrials,c1))];
                preMoveVar=[preMoveVar var(respAmpsMovePre(1:minnumtrials,c1))];
                decodingCellResponsiveness=[decodingCellResponsiveness cellResponsiveness(c1)];
                
                MdlImmob = fitcsvm(respAmpsMatImmobile(:,c1),tags1);
                cvmodelImmob = crossval(MdlImmob);
                LImmob = kfoldLoss(cvmodelImmob);
                
                MdlMove = fitcsvm(respAmpsMatMove(:,c1),tags1);
                cvmodelMove = crossval(MdlMove);
                LMove = kfoldLoss(cvmodelMove);
                
                MdlState = fitcsvm(respAmpsMatState(:,c1),tags1);
                cvmodelState = crossval(MdlState);
                LState = kfoldLoss(cvmodelState);
                
                MdlStateAndStim = fitcsvm(respAmpsALLMat(:,c1),tagsAll);
                cvmodelStateAndStim = crossval(MdlStateAndStim);
                LStateAndStim = kfoldLoss(cvmodelStateAndStim);
                
                curCellImmobSHUF=[];
                curCellMoveSHUF=[];
                curCellStateSHUF=[];
                curCellStateAndStimSHUF=[];
                for xx1=1:numshufs1
                    MdlImmobSHUF = fitcsvm(respAmpsMatImmobile(:,c1),tags1(randperm(length(tags1))));
                    cvmodelImmobSHUF = crossval(MdlImmobSHUF);
                    LImmobSHUF = kfoldLoss(cvmodelImmobSHUF);
                    curCellImmobSHUF=[curCellImmobSHUF LImmobSHUF];
                    
                    MdlMoveSHUF = fitcsvm(respAmpsMatMove(:,c1),tags1(randperm(length(tags1))));
                    cvmodelMoveSHUF = crossval(MdlMoveSHUF);
                    LMoveSHUF = kfoldLoss(cvmodelMoveSHUF);
                    curCellMoveSHUF=[curCellMoveSHUF LMoveSHUF];
                    
                    MdlStateSHUF = fitcsvm(respAmpsMatState(:,c1),tags1(randperm(length(tags1))));
                    cvmodelStateSHUF = crossval(MdlStateSHUF);
                    LStateSHUF = kfoldLoss(cvmodelStateSHUF);
                    curCellStateSHUF=[curCellStateSHUF LStateSHUF];
                    
                    MdlStateAndStimSHUF = fitcsvm(respAmpsALLMat(:,c1),tagsAll(randperm(length(tagsAll))));
                    cvmodelStateAndStimSHUF = crossval(MdlStateAndStimSHUF);
                    LStateAndStimSHUF = kfoldLoss(cvmodelStateAndStimSHUF);
                    curCellStateAndStimSHUF=[curCellStateAndStimSHUF LStateAndStimSHUF];
                    
                end
                
                stimDetectionErrorImmobSingleCell=[stimDetectionErrorImmobSingleCell LImmob];
                stimDetectionErrorMoveSingleCell=[stimDetectionErrorMoveSingleCell LMove];
                stateDetectionErrorSingleCell=[stateDetectionErrorSingleCell LState];
                stateAndStimDetectionErrorSingleCell=[stateAndStimDetectionErrorSingleCell LStateAndStim];
                
                
                stimDetectionPVALImmobSingleCell=[stimDetectionPVALImmobSingleCell mean(LImmob>curCellImmobSHUF)];
                stimDetectionPVALMoveSingleCell=[stimDetectionPVALMoveSingleCell mean(LMove>curCellMoveSHUF)];
                stateDetectionPVALSingleCell=[stateDetectionPVALSingleCell mean(LState>curCellStateSHUF)];
                stateAndStimDetectionPVALSingleCell=[stateAndStimDetectionPVALSingleCell mean(LStateAndStim>curCellStateAndStimSHUF)];
                
                ['Batch ' allresults{1,i}{1,1}.batch ' Cell ' num2str(c1)]
                allCellIndsForDecoding{end+1}= [allresults{1,i}{1,1}.batch ' Cell ' num2str(c1)];
                ['Immob: ' num2str(LImmob) ' Move: ' num2str(LMove) ' State during: ' num2str(LState) ' StateAndStim: ' num2str(LStateAndStim)]
             
            end
            
            end
            
            % ENSEMBLE DETECTION
            
            % stim detection in immobility
            MdlImmobEnsemble = fitcsvm(respAmpsMatImmobile,tags1);
            cvmodelImmobEnsemble = crossval(MdlImmobEnsemble);
            immobileStimDecodingError = kfoldLoss(cvmodelImmobEnsemble);
            
            % stim detection in locomotion
            MdlMoveEnsemble = fitcsvm(respAmpsMatMove,tags1);
            cvmodelMoveEnsemble = crossval(MdlMoveEnsemble);
            locomotionStimDecodingError = kfoldLoss(cvmodelMoveEnsemble);
            
            % state detection from spont
            preAmpsMatImmobile=respAmpsImmobilePre(1:minnumtrials,:);
            preAmpsMatMove=respAmpsMovePre(1:minnumtrials,:);
            preAll=[preAmpsMatImmobile;preAmpsMatMove];
            
            MdlStateEnsembleFromPre = fitcsvm(preAll,tags1);
            cvmodelStateEnsembleFromPre = crossval(MdlStateEnsembleFromPre);
            prewinStateDecodingError = kfoldLoss(cvmodelStateEnsembleFromPre);
            
            
            % State detection from stim
            respAmpsMatImmobileDuring=respAmpsImmobileDuring(1:minnumtrials,:);
            respAmpsMatMoveDuring=respAmpsMoveDuring(1:minnumtrials,:);
            respAll=[respAmpsMatImmobileDuring;respAmpsMatMoveDuring];
            
            MdlStateEnsembleFromResp = fitcsvm(respAll,tags1);
            cvmodelStateEnsembleFromResp = crossval(MdlStateEnsembleFromResp);
            stimwinStateDecodingError = kfoldLoss(cvmodelStateEnsembleFromResp);
            
            % state and stim detection
            
            stateAndStim=[respAmpsImmobilePre(1:minnumtrials,:); respAmpsImmobileDuring(1:minnumtrials,:);respAmpsMovePre(1:minnumtrials,:); respAmpsMoveDuring(1:minnumtrials,:)];
            tagsstateAndStim=[zeros(1,minnumtrials) ones(1,minnumtrials) 2*ones(1,minnumtrials) 3*ones(1,minnumtrials)];
            MdlStateAndStim = fitcdiscr(stateAndStim,tagsstateAndStim);
            cvmodelStateAndStim  = crossval(MdlStateAndStim);
            stimAndStateDecodingError = kfoldLoss(cvmodelStateAndStim);
            
            
            
            shufvalsImmob=[];
            shufvalsMove=[];
            shufStateVals=[];
            shufStateValsD=[];
            shufStimAndStateDecodingErrors=[];
            for i2=1:numshufs1
                
                % shuffled stim detection in immobility
                Mdl = fitcsvm(respAmpsMatImmobile,tags1(randperm(length(tags1))));
                cvmodel = crossval(Mdl);
                L = kfoldLoss(cvmodel);
                shufvalsImmob=[shufvalsImmob L];
                
                % shuffled stim detection in locomotion
                Mdl = fitcsvm(respAmpsMatMove,tags1(randperm(length(tags1))));
                cvmodel = crossval(Mdl);
                L = kfoldLoss(cvmodel);
                shufvalsMove=[shufvalsMove L];
                
                % shuffled state detection from spont
                Mdl = fitcsvm(preAll,tags1(randperm(length(tags1))));
                cvmodel = crossval(Mdl);
                L = kfoldLoss(cvmodel);
                shufStateVals=[shufStateVals L];
                
                % shuffled state detection from stim
                Mdl = fitcsvm(respAll,tags1(randperm(length(tags1))));
                cvmodel = crossval(Mdl);
                L = kfoldLoss(cvmodel);
                shufStateValsD=[shufStateValsD L];
                
                % shuffled state and stim detection
                MdlStateAndStimS = fitcdiscr(stateAndStim,tagsstateAndStim(randperm(length(tagsstateAndStim))));
                cvmodelStateAndStimS  = crossval(MdlStateAndStimS);
                stimAndStateDecodingErrorS = kfoldLoss(cvmodelStateAndStimS);
                shufStimAndStateDecodingErrors=[shufStimAndStateDecodingErrors stimAndStateDecodingErrorS];
                
                
            end
            
            
            
            AllImmobileStimDecodingError=[AllImmobileStimDecodingError immobileStimDecodingError];
            AllImmobileStimDecodingErrorSHUF=[AllImmobileStimDecodingErrorSHUF mean(shufvalsImmob)];
            curPVal=mean(immobileStimDecodingError>shufvalsImmob);
            AllImmobileStimDecodingPVAL=[AllImmobileStimDecodingPVAL curPVal];
            
            
            AllLocomotionStimDecodingError=[AllLocomotionStimDecodingError locomotionStimDecodingError];
            AllLocomotionStimDecodingErrorSHUF=[AllLocomotionStimDecodingErrorSHUF  mean(shufvalsMove)];
            curPVal=mean(locomotionStimDecodingError>shufvalsMove);
            AllLocomotionStimDecodingPVAL=[AllLocomotionStimDecodingPVAL curPVal];
            
            
            AllPreWinStateDecodingError=[AllPreWinStateDecodingError prewinStateDecodingError];
            AllPreWinStateDecodingErrorSHUF=[AllPreWinStateDecodingErrorSHUF mean(shufStateVals)];
            curPVal=mean(prewinStateDecodingError>shufStateVals);
            AllPreWinStateDecodingPVAL=[AllPreWinStateDecodingPVAL curPVal];
            
            
            AllStimWinStateDecodingError=[AllStimWinStateDecodingError stimwinStateDecodingError];
            AllStimWinStateDecodingErrorSHUF=[AllStimWinStateDecodingErrorSHUF mean(shufStateValsD)];
            curPVal=mean(stimwinStateDecodingError>shufStateValsD);
            AllStimWinStateDecodingPVAL=[AllStimWinStateDecodingPVAL curPVal];
            
            
            AllStimAndStateDecodingError=[AllStimAndStateDecodingError stimAndStateDecodingError];
            AllStimAndStateDecodingErrorShuf=[AllStimAndStateDecodingErrorShuf mean(shufStimAndStateDecodingErrors)];
            curPVal=mean(stimAndStateDecodingError>shufStimAndStateDecodingErrors);
            AllStimAndStateDecodingPVAL=[AllStimAndStateDecodingPVAL curPVal];
            
            
        else
            participatingEnsemblesDecoding=[participatingEnsemblesDecoding 0];
            
        end
    end
    
    
    curXpixels=allresults{1,i}{1,1}.cellPixels_X;
    curYpixels=allresults{1,i}{1,1}.cellPixels_Y;
    
    %Spatial organization of locomotion modulation
    animStr=allresults{1,i}{1,1}.animal;
    expBatch=allresults{1,i}{1,1}.batch;
    plotspatial=0;
    if plotspatial
        
        %all cell locations
        figure('Position', [50, 250, 1400, 300])
        subplot(1,4,1)
        imagesc(allresults{1,i}{1,1}.meanImg);colormap(gray)
        set(gca,'FontSize',6)
        title(['Animal ' animStr ' batch ' expBatch])
        %all cell ROIs
        subplot(1,4,2)
        im=zeros(256,256);
        for j=1:numcells
            im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= 1;
        end
        imagesc(im);colormap(gray)
        
        %ROIS sig mod by locomotion
        
        if length(difSpontLocomotionVsImmobility)>0
            subplot(1,4,3)
            im=zeros(256,256);
            for j=1:numcells
                if difSpontLocomotionVsImmobility(j) == 1
                    im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= 1;
                elseif difSpontLocomotionVsImmobility(j)==-1
                    im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= -1;
                else
                    im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= 0;
                    
                end
            end
            imagesc(im);colormap(gray)
            title(['Loc-mod. of spont, relDPos= ' num2str(round(relativeDistPos,2)) ' relDNeg= ' num2str(round(relativeDistNeg,2)) ' ProbPos= ' num2str(round(allProportionOfPosSpeedModCellsInEnsemble(end),2))])
            
            subplot(1,4,4)
            im=zeros(256,256);
            for j=1:numcells
                if difEvokedMinusPreLocomotionVsImmobility(j) == 1
                    im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= 1;
                elseif difEvokedMinusPreLocomotionVsImmobility(j)==-1
                    im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= -1;
                else
                    im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= 0;
                    
                end
            end
            imagesc(im);colormap(gray)
            title('Loc-mod. of snd S/N')
            
            
            plotcellnums=0;
            if plotcellnums
                
                for ii=1:numcells
                    meanx=allresults{1,i}{1,1}.cellCenters(ii,1);
                    meany=allresults{1,i}{1,1}.cellCenters(ii,2);
                    subplot(1,4,1)
                    text(meanx-5,meany,num2str(ii),'color','magenta','fontweight','bold','fontsize',6);
                    set(gca,'FontSize',6)
                    
                    subplot(1,4,2)
                    text(meanx-5,meany,num2str(ii),'color','magenta','fontweight','bold','fontsize',6);
                    set(gca,'FontSize',6)
                    
                    subplot(1,4,3)
                    text(meanx-5,meany,num2str(ii),'color','magenta','fontweight','bold','fontsize',6);
                    set(gca,'FontSize',6)
                    
                    subplot(1,4,4)
                    text(meanx-5,meany,num2str(ii),'color','magenta','fontweight','bold','fontsize',6);
                    set(gca,'FontSize',6)
                    
                    
                    
                end
            end
            animStr=allresults{1,i}{1,1}.animal;
            expBatch=allresults{1,i}{1,1}.batch;
            dff=allresults{1,i}{1,1}.dff;
            
            cc=[1 0 0; 1 1 1;0 1 0];
            caxis([-1 1])
            colormap(cc)
        end
        savespatial=0;
        if savespatial
            saveas(gcf,[spatialfolder  allresults{1,i}{1,1}.batch 'color'],'jpg');
            saveas(gcf,[spatialfolder  allresults{1,i}{1,1}.batch 'color'],'fig');
            print([spatialfolder allresults{1,i}{1,1}.batch 'color'],'-dpdf','-bestfit');
        end
        
        if ~isempty(epochdffspeedcorrs)
            
            
            figure;
            
            im=-100*ones(256,256);
            cc=colormap('jet');
            cc=[1 1 1;cc];
            
            for j=1:numcells
                
                im(sub2ind(size(im),curYpixels{j},curXpixels{j}))= epochdffspeedcorrs(j);
                
            end
            imagesc(im);
            
            caxis([-0.5 0.5])
            colormap(cc)
            colorbar
            title(['dff-speed correlations ' allresults{1,i}{1,1}.batch] )
            savespatial2=1;
            if savespatial2
                saveas(gcf,[spatialfolder2  allresults{1,i}{1,1}.batch 'color'],'jpg');
                saveas(gcf,[spatialfolder2  allresults{1,i}{1,1}.batch 'color'],'fig');
                print([spatialfolder2 allresults{1,i}{1,1}.batch 'color'],'-dpdf','-bestfit');
            end
            
        end
    end
    %
    printconttraces=0;
    if printconttraces
        
        figure('Position', [800, 250, 900, 500])
        bar(find(allmov==1),numcells*1.2*ones(size(find(allmov==1))),1,'FaceColor',[1 1 1]*0.8,'EdgeColor',[1 1 1]*0.8)
        hold on;
        for ii=1:numcells,plot(alldff(ii,:)+ii,'linewidth',1);hold on;end
        set(gca,'FontSize',6)
        
        saveas(gcf,[conttracesfolder  allresults{1,i}{1,1}.batch],'jpg');
        saveas(gcf,[conttracesfolder  allresults{1,i}{1,1}.batch],'fig');
        print([conttracesfolder allresults{1,i}{1,1}.batch],'-dpdf','-bestfit');
    end
    
    
    
    %==========
    
    plotNATensemble1=0;
    if plotNATensemble1
        numNATfigs=ceil(numcells/10);
        numstim=size(NATresps,2);
        numcolumns=numstim;
        curnumcells=0;
        c1=0;
        
        for c=1:numcells
            respacrossstim=[];
            for jjj=1:numstim
                respacrossstim=[respacrossstim allCellResponsiveness_perstim_ensemble(c,jjj)];
            end
            if sum(respacrossstim)==0
                continue
            end
            curnumcells=curnumcells+1;
            if curnumcells==1
                figure('position',[30 30 1000 700]);
            end
            for s=1:numstim
                subplot(10,numcolumns,numstim*(curnumcells-1)+s);
                currespMOV=squeeze(NATrespsMOV(c,s,:,:));
                currespImmob=squeeze(NATrespsImmobile(c,s,:,:));
                shadedErrorBar(1:size(currespMOV,2),nanmean(currespMOV),nanstd(currespMOV)/sqrt(size(currespMOV,1)),'lineprops','r')
                
                hold on;
                shadedErrorBar(1:size(currespImmob,2),nanmean(currespImmob),nanstd(currespImmob)/sqrt(size(currespImmob,1)),'lineprops','g')
                
                axis tight
                if allcell_stimSNRmodulation_nat_ensemble(c,s)==1
                    title(['Cell=' num2str(c) ' DFF= ' num2str(allcell_stimdFFmodulation_nat_ensemble(c,s)) ' SNR= ' num2str(allcell_stimSNRmodulation_nat_ensemble(c,s)) ' Resp= ' num2str(allCellResponsiveness_perstim_ensemble(c,s))],'color','r');
                elseif allcell_stimSNRmodulation_nat_ensemble(c,s)==-1
                    title(['Cell=' num2str(c) ' DFF= ' num2str(allcell_stimdFFmodulation_nat_ensemble(c,s)) ' SNR= ' num2str(allcell_stimSNRmodulation_nat_ensemble(c,s)) ' Resp= ' num2str(allCellResponsiveness_perstim_ensemble(c,s))],'color','g');
                else
                    title(['Cell=' num2str(c) ' DFF= ' num2str(allcell_stimdFFmodulation_nat_ensemble(c,s)) ' SNR= ' num2str(allcell_stimSNRmodulation_nat_ensemble(c,s)) ' Resp= ' num2str(allCellResponsiveness_perstim_ensemble(c,s))],'color','k');
                    
                end
                
                set(gca,'FontSize',6)
                %   end
            end
            
            sameYforall=1
            if sameYforall
                maxY=0.3;
                minY=-0.25;
                for ss1=1:numstim
                    subplot(10,numcolumns,numstim*(curnumcells-1)+ss1);
                    curY=ylim;
                    if curY(2)>maxY, maxY=curY(2);end
                    if curY(1)<minY, minY=curY(1);end
                end
                for ss1=1:numstim
                    subplot(10,numcolumns,numstim*(curnumcells-1)+ss1);
                    ylim([minY maxY])
                    xlim([1 12])
                end
            end
            if curnumcells==10 | c==numcells
                
                saveas(gcf,[NATimagesfolder1 'Natstimresps'  allresults{1,i}{1,1}.batch '-' num2str(c1) 'b'],'jpg');
                saveas(gcf,[NATimagesfolder1 'Natstimresps' allresults{1,i}{1,1}.batch '-' num2str(c1) 'b'],'fig');
                print([NATimagesfolder1 'Natstimresps'  allresults{1,i}{1,1}.batch '-' num2str(c1) 'b'],'-dpdf','-fillpage');
                c1=c1+1
                curnumcells=0;
            end
            
        end
        
    end
  
    plotNATensemble2=0;
    if plotNATensemble2
        numNATfigs=ceil(numcells/10);
        numstim=size(NATresps,2);
        numcolumns=numstim;
        curnumcells=0;
        xpt = 1:13;
        for c=1:numcells
            curnumcells=curnumcells+1;
            if curnumcells==1
                figure('position',[60 60 1000 700]);
            end
            for s=1:numstim
                xpt1=find(~isnan(squeeze(nanmean(NATrespsMOV(1,s,:,:)))));
                xpt2=find(~isnan(squeeze(nanmean(NATrespsImmobile(1,s,:,:)))));
                xptcommon=intersect(xpt1,xpt2)';
                subplot(10,numcolumns,numstim*(curnumcells-1)+s);
                patch([xptcommon fliplr(xptcommon)],[nanmean(squeeze(NATrespsMOV(c,s,:,xptcommon)))-nanstd(squeeze(NATrespsMOV(c,s,:,xptcommon)))/sqrt(size(NATrespsMOV,3)), fliplr([nanmean(squeeze(NATrespsMOV(c,s,:,xptcommon)))+ nanstd(squeeze(NATrespsMOV(c,s,:,xptcommon)))/sqrt(size(NATrespsMOV,3))])],'r','LineStyle','none');
                hold on; patch([xptcommon fliplr(xptcommon)],[nanmean(squeeze(NATrespsImmobile(c,s,:,xptcommon)))-nanstd(squeeze(NATrespsImmobile(c,s,:,xptcommon)))/sqrt(size(NATrespsImmobile,3)), fliplr([nanmean(squeeze(NATrespsImmobile(c,s,:,xptcommon)))+ nanstd(squeeze(NATrespsImmobile(c,s,:,xptcommon)))/sqrt(size(NATrespsImmobile,3))])],'g','LineStyle','none');
                hold on; plot(xptcommon,nanmean(squeeze(NATrespsMOV(c,s,:,xptcommon))),'k');
                hold on; plot(xptcommon,nanmean(squeeze(NATrespsImmobile(c,s,:,xptcommon))),'k');
                axis tight
                % ylim([-0.1 1]);
                title(['Cell=' num2str(c) ' DFF= ' num2str(allcell_stimdFFmodulation_nat_ensemble(c,s)) ' SNR= ' num2str(allcell_stimSNRmodulation_nat_ensemble(c,s))]);
                set(gca,'FontSize',6)
            end
            if curnumcells==10
                
                curnumcells=0;
                saveas(gcf,[NATimagesfolder2 'ANatstimresps' allresults{1,i}{1,1}.batch num2str(c)],'jpg');
                saveas(gcf,[NATimagesfolder2 'ANatstimresps' allresults{1,i}{1,1}.batch num2str(c)],'fig');
                print([NATimagesfolder2 'ANatstimresps' allresults{1,i}{1,1}.batch num2str(c)],'-dpdf','-fillpage');
                
            end
            
            
        end
        close all
    end
    
    
    close all
end

%% num cells per ensemble in decoding
mean(allnumcellsperensembledecoding)
std(allnumcellsperensembledecoding)/sqrt(length(allnumcellsperensembledecoding))



%%
load dataAfterNewImagingRun041422




%% Scatter of evoked-minus-pre in locomotion vs. immobility of IMMOBILITY-SOUND-RESPONSIVE CELLS


figure('Position',[50,50,800,320]);


subplot(1,2,1)
plot(allmeanEvokedMinusPreImmobility_responsiveIMMOB(difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB==0), allmeanEvokedMinusPreLocomotion_responsiveIMMOB(difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB==0),'ko','markerfacecolor','b')
hold on
plot(allmeanEvokedMinusPreImmobility_responsiveIMMOB(difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB == -1), allmeanEvokedMinusPreLocomotion_responsiveIMMOB(difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB == -1),'ko','markerfacecolor','r')
plot(allmeanEvokedMinusPreImmobility_responsiveIMMOB(difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB == 1), allmeanEvokedMinusPreLocomotion_responsiveIMMOB(difEvokedMinusPreLocomotionVsImmobilityAllExps_onlyResponsiveIMMOB == 1),'ko','markerfacecolor','g')


plot([-1 4.5],[-1 4.5],'r--')
axis([-1 4.5 -1 4.5])
xlabel('Evoked-minus-pre in Immobility')
ylabel('Evoked-minus-pre in Locomotion')
title('evoked-minus-pre in locomotion vs. immobility of SOUND-RESPONSIVE CELLS')

subplot(1,2,2)
d4=allmeanEvokedMinusPreLocomotion_responsiveIMMOB-allmeanEvokedMinusPreImmobility_responsiveIMMOB;
boxplot([d4],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.3 0.3])


set(gca,'xticklabel',{'Locomotion minus immobility'})

p2=signrank(allmeanEvokedMinusPreImmobility_responsiveIMMOB,allmeanEvokedMinusPreLocomotion_responsiveIMMOB)
title(['P= ' num2str(p2)])

%% nat (supplemental)

figure('Position',[50,50,800,320]);

subplot(1,2,1)

plot(allmeanEvokedMinusPreImmobilityNAT,allmeanEvokedMinusPreLocomotionNAT,'ko');hold on;plot([-1 7],[-1 7],'r--')
hold on
plot(allmeanEvokedMinusPreImmobilityNAT(difEvokedMinusPreLocomotionVsImmobilityNAT==0),allmeanEvokedMinusPreLocomotionNAT(difEvokedMinusPreLocomotionVsImmobilityNAT==0),'ko','markerfacecolor','b')
plot(allmeanEvokedMinusPreImmobilityNAT(difEvokedMinusPreLocomotionVsImmobilityNAT==1),allmeanEvokedMinusPreLocomotionNAT(difEvokedMinusPreLocomotionVsImmobilityNAT==1),'ko','markerfacecolor','g')
plot(allmeanEvokedMinusPreImmobilityNAT(difEvokedMinusPreLocomotionVsImmobilityNAT==-1),allmeanEvokedMinusPreLocomotionNAT(difEvokedMinusPreLocomotionVsImmobilityNAT==-1),'ko','markerfacecolor','r')
xlabel('Evoked-minus-pre in Immobility')
ylabel('Evoked-minus-pre in Locomotion')
axis([-.5 6 -.5 6])
p2=signrank(allmeanEvokedMinusPreImmobilityNAT,allmeanEvokedMinusPreLocomotionNAT)
subplot(1,2,2)
d4=allmeanEvokedMinusPreLocomotionNAT-allmeanEvokedMinusPreImmobilityNAT;
boxplot([d4],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.3 0.3])


%% new tones
figure('Position',[50,50,800,320]);

subplot(1,2,1)

plot(allmeanEvokedMinusPreImmobilityTON3,allmeanEvokedMinusPreLocomotionTON3,'ko');hold on;plot([-0.5 4],[-0.5 4],'r--')
hold on
plot(allmeanEvokedMinusPreImmobilityTON3(difEvokedMinusPreLocomotionVsImmobilityTON3==0),allmeanEvokedMinusPreLocomotionTON3(difEvokedMinusPreLocomotionVsImmobilityTON3==0),'ko','markerfacecolor','b')
plot(allmeanEvokedMinusPreImmobilityTON3(difEvokedMinusPreLocomotionVsImmobilityTON3==1),allmeanEvokedMinusPreLocomotionTON3(difEvokedMinusPreLocomotionVsImmobilityTON3==1),'ko','markerfacecolor','g')
plot(allmeanEvokedMinusPreImmobilityTON3(difEvokedMinusPreLocomotionVsImmobilityTON3==-1),allmeanEvokedMinusPreLocomotionTON3(difEvokedMinusPreLocomotionVsImmobilityTON3==-1),'ko','markerfacecolor','r')
xlabel('Evoked-minus-pre in Immobility')
ylabel('Evoked-minus-pre in Locomotion')
p2=signrank(allmeanEvokedMinusPreImmobilityTON3,allmeanEvokedMinusPreLocomotionTON3)
axis([-.2 4 -.2 4])
subplot(1,2,2)
d2=allmeanEvokedMinusPreLocomotionTON3-allmeanEvokedMinusPreImmobilityTON3;
boxplot([d2],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.3 0.3])







%% state-dependent noise correlations, boxplot
groupinds=[ 1*ones(1,length(alllocomotionEnhancedNoiseCorrs)) 2*ones(1,length(alllocomotionSuppressedNoiseCorrs)) 3*ones(1,length(allacrossClassesNoiseCorrs)) ];

figure;
subplot(1,3,1);
boxplot([ alllocomotionEnhancedNoiseCorrs alllocomotionSuppressedNoiseCorrs allacrossClassesNoiseCorrs],groupinds,'Symbol','')
ylabel('Noise correlation')
hold on;
plot([0 5],[0 0],'r--')
ylim([-0.7 0.9])
title('ACROSS STATES')

subplot(1,3,2)
groupindsMove=[0*ones(1,length(alllocomotionEnhancedNoiseCorrsMove)) 1*ones(1,length(alllocomotionSuppressedNoiseCorrsMove)) 2*ones(1,length(allacrossClassesNoiseCorrsMove)) ];
boxplot([alllocomotionEnhancedNoiseCorrsMove alllocomotionSuppressedNoiseCorrsMove allacrossClassesNoiseCorrsMove ],groupindsMove,'Symbol','')
ylabel('Noise correlation')
hold on;
plot([0 5],[0 0],'r--')
ylim([-0.7 0.9])
title('MOVEMENT')


subplot(1,3,3)
groupindsImmob=[3*ones(1,length(alllocomotionEnhancedNoiseCorrsImmob)) 4*ones(1,length(alllocomotionSuppressedNoiseCorrsImmob))  5*ones(1,length(allacrossClassesNoiseCorrsImmob))];
boxplot([ alllocomotionEnhancedNoiseCorrsImmob alllocomotionSuppressedNoiseCorrsImmob  allacrossClassesNoiseCorrsImmob],groupindsImmob,'Symbol','')
ylabel('Noise correlation')
hold on;
plot([0 5],[0 0],'r--')
ylim([-0.7 0.9])
title('IMMOBILITY')

%% state-dependent noise correlations, anova

groupindsAll=[ ones(1,length(alllocomotionEnhancedNoiseCorrs)), 2*ones(1,length(alllocomotionSuppressedNoiseCorrs)), 3*ones(1,length(allacrossClassesNoiseCorrs))];
dataall=[ alllocomotionEnhancedNoiseCorrs alllocomotionSuppressedNoiseCorrs allacrossClassesNoiseCorrs];
[p1 tbl1 stats1]=anova1(dataall,groupindsAll);
multcompare(stats1)

groupindsMove=[0*ones(1,length(alllocomotionEnhancedNoiseCorrsMove)) 1*ones(1,length(alllocomotionSuppressedNoiseCorrsMove)) 2*ones(1,length(allacrossClassesNoiseCorrsMove)) ];
dataMove=[alllocomotionEnhancedNoiseCorrsMove alllocomotionSuppressedNoiseCorrsMove allacrossClassesNoiseCorrsMove ];
[p2 tbl2 stats2]=anova1(dataMove,groupindsMove);
multcompare(stats2)

groupindsImmob=[3*ones(1,length(alllocomotionEnhancedNoiseCorrsImmob)) 4*ones(1,length(alllocomotionSuppressedNoiseCorrsImmob))  5*ones(1,length(allacrossClassesNoiseCorrsImmob))];
dataImmob=[ alllocomotionEnhancedNoiseCorrsImmob alllocomotionSuppressedNoiseCorrsImmob  allacrossClassesNoiseCorrsImmob];
[p3 tbl3 stats3]=anova1(dataImmob,groupindsImmob);
multcompare(stats3)


%% SHUFFLED state-dependent noise correlations, boxplot
groupinds=[ 1*ones(1,length(alllocomotionEnhancedNoiseCorrsShuf)) 2*ones(1,length(alllocomotionSuppressedNoiseCorrsShuf)) 3*ones(1,length(allacrossClassesNoiseCorrsShuf)) ];

figure;
subplot(1,3,1);
boxplot([ alllocomotionEnhancedNoiseCorrsShuf alllocomotionSuppressedNoiseCorrsShuf allacrossClassesNoiseCorrsShuf],groupinds,'Symbol','')
ylabel('Noise correlation')
hold on;
plot([0 5],[0 0],'r--')
ylim([-0.7 0.9])
title('ACROSS STATES')

subplot(1,3,2)
groupindsMove=[0*ones(1,length(alllocomotionEnhancedNoiseCorrsMoveShuf)) 1*ones(1,length(alllocomotionSuppressedNoiseCorrsMoveShuf)) 2*ones(1,length(allacrossClassesNoiseCorrsMoveShuf)) ];
boxplot([alllocomotionEnhancedNoiseCorrsMoveShuf alllocomotionSuppressedNoiseCorrsMoveShuf allacrossClassesNoiseCorrsMoveShuf ],groupindsMove,'Symbol','')
ylabel('Noise correlation')
hold on;
plot([0 5],[0 0],'r--')
ylim([-0.7 0.9])
title('MOVEMENT')


subplot(1,3,3)
groupindsImmob=[3*ones(1,length(alllocomotionEnhancedNoiseCorrsImmobShuf)) 4*ones(1,length(alllocomotionSuppressedNoiseCorrsImmobShuf))  5*ones(1,length(allacrossClassesNoiseCorrsImmobShuf))];
boxplot([ alllocomotionEnhancedNoiseCorrsImmobShuf alllocomotionSuppressedNoiseCorrsImmobShuf  allacrossClassesNoiseCorrsImmobShuf],groupindsImmob,'Symbol','')
ylabel('Noise correlation')
hold on;
plot([0 5],[0 0],'r--')
ylim([-0.7 0.9])
title('IMMOBILITY')

%% state-dependent cross-corr

figure;
subplot(1,3,1)
shadedErrorBar(1:size(alllocomotionEnhancedPairsXCorrs,2),nanmean(alllocomotionEnhancedPairsXCorrs),nanstd(alllocomotionEnhancedPairsXCorrs)/sqrt(size(alllocomotionEnhancedPairsXCorrs,1)),'lineprops','g')
hold on
shadedErrorBar(1:size(alllocomotionSuppressedPairsXCorrs,2),nanmean(alllocomotionSuppressedPairsXCorrs),nanstd(alllocomotionSuppressedPairsXCorrs)/sqrt(size(alllocomotionSuppressedPairsXCorrs,1)),'lineprops','r')
shadedErrorBar(1:size(allacrossClassesXCorrs,2),nanmean(allacrossClassesXCorrs),nanstd(allacrossClassesXCorrs)/sqrt(size(allacrossClassesXCorrs,1)),'lineprops','k')

ylabel('Cross corr')
axis([12 30 0 0.28])
title('ACROSS STATES')


subplot(1,3,2)
shadedErrorBar(1:size(alllocomotionEnhancedPairsXCorrsMove,2),nanmean(alllocomotionEnhancedPairsXCorrsMove),nanstd(alllocomotionEnhancedPairsXCorrsMove)/sqrt(size(alllocomotionEnhancedPairsXCorrsMove,1)),'lineprops','g')
hold on
shadedErrorBar(1:size(alllocomotionSuppressedPairsXCorrsMove,2),nanmean(alllocomotionSuppressedPairsXCorrsMove),nanstd(alllocomotionSuppressedPairsXCorrsMove)/sqrt(size(alllocomotionSuppressedPairsXCorrsMove,1)),'lineprops','r')
shadedErrorBar(1:size(allacrossClassesXCorrsMove,2),nanmean(allacrossClassesXCorrsMove),nanstd(allacrossClassesXCorrsMove)/sqrt(size(allacrossClassesXCorrsMove,1)),'lineprops','k')

ylabel('Cross corr')
axis([12 30 0 0.28])
title('MOVEMENT')

subplot(1,3,3)
shadedErrorBar(1:size(alllocomotionEnhancedPairsXCorrsImmob,2),nanmean(alllocomotionEnhancedPairsXCorrsImmob),nanstd(alllocomotionEnhancedPairsXCorrsImmob)/sqrt(size(alllocomotionEnhancedPairsXCorrsImmob,1)),'lineprops','g')
hold on
shadedErrorBar(1:size(alllocomotionSuppressedPairsXCorrsImmob,2),nanmean(alllocomotionSuppressedPairsXCorrsImmob),nanstd(alllocomotionSuppressedPairsXCorrsImmob)/sqrt(size(alllocomotionSuppressedPairsXCorrsImmob,1)),'lineprops','r')
shadedErrorBar(1:size(allacrossClassesXCorrsImmob,2),nanmean(allacrossClassesXCorrsImmob),nanstd(allacrossClassesXCorrsImmob)/sqrt(size(allacrossClassesXCorrsImmob,1)),'lineprops','k')

ylabel('Cross corr')
axis([12 30 0 0.28])
title('IMMOBILITY')
%% ANOVA

% ALL
sumvalsacross=sum(allacrossClassesXCorrs(:,19:23)');
sumvalssuppressed=sum(alllocomotionSuppressedPairsXCorrs(:,19:23)');
sumvalsenhanced=sum(alllocomotionEnhancedPairsXCorrs(:,19:23)');

valsAllCC=[ sumvalsenhanced sumvalssuppressed sumvalsacross];
tagsAllCC=[ ones(1,length(sumvalsenhanced)), 2*ones(1,length(sumvalssuppressed)), 3*ones(1,length(sumvalsacross))];
[pALLCC tblALLCC statsALLCC]=anova1(valsAllCC,tagsAllCC);
multcompare(statsALLCC)


% Move
sumvalsacrossMove=sum(allacrossClassesXCorrsMove(:,19:23)');
sumvalssuppressedMove=sum(alllocomotionSuppressedPairsXCorrsMove(:,19:23)');
sumvalsenhancedMove=sum(alllocomotionEnhancedPairsXCorrsMove(:,19:23)');

valsMove=[ sumvalsenhancedMove sumvalssuppressedMove sumvalsacrossMove];
tagsMove=[ones(1,length(sumvalsenhancedMove)), 2*ones(1,length(sumvalssuppressedMove)), 3*ones(1,length(sumvalsacrossMove))];
[pMove1 tblMove1 statsMove1]=anova1(valsMove,tagsMove);
multcompare(statsMove1)

%Immobility
sumvalsacrossImmob=sum(allacrossClassesXCorrsImmob(:,19:23)');
sumvalssuppressedImmob=sum(alllocomotionSuppressedPairsXCorrsImmob(:,19:23)');
sumvalsenhancedImmob=sum(alllocomotionEnhancedPairsXCorrsImmob(:,19:23)');

valsImmob=[ sumvalsenhancedImmob sumvalssuppressedImmob sumvalsacrossImmob];
tagsImmob=[ones(1,length(sumvalsenhancedImmob)), 2*ones(1,length(sumvalssuppressedImmob)), 3*ones(1,length(sumvalsacrossImmob))];
[pImmob1 tblImmob1 statsImmob1]=anova1(valsImmob,tagsImmob);
multcompare(statsImmob1)

%% Population-average PSTHs across states
figure;
shadedErrorBar(1:size(allImmobPSTHsResponsiveIMMOB,2),nanmean(allImmobPSTHsResponsiveIMMOB),nanstd(allImmobPSTHsResponsiveIMMOB)/sqrt(size(allImmobPSTHsResponsiveIMMOB,1)),'lineprops','r')
hold on
shadedErrorBar(1:size(allMovePSTHsResponsiveIMMOB,2),nanmean(allMovePSTHsResponsiveIMMOB),nanstd(allMovePSTHsResponsiveIMMOB)/sqrt(size(allMovePSTHsResponsiveIMMOB,1)),'lineprops','g')
xlim([1 13])
%% SPONT boxplots, move minus immob
figure('Position',[50,50,300,600]);
subplot(1,2,2)
boxplot([allmeanspontLocomotion-allmeanspontImmobility],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.4 0.4])
ylabel('move minus immob')
title('all cells')
p2=signrank(allmeanspontLocomotion,allmeanspontImmobility)

subplot(1,2,1)
boxplot([allmeanspontLocomotion_responsiveIMMOB-allmeanspontImmobility_responsiveIMMOB],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.4 0.4])
ylabel('move minus immob')
title('significant cells')
p3=signrank(allmeanspontLocomotion_responsiveIMMOB,allmeanspontImmobility_responsiveIMMOB)



%% locomotion influence on pre-stim window against locomotion influence on stim window
%
dPre=allmeanspontLocomotion_responsiveIMMOB-allmeanspontImmobility_responsiveIMMOB;
dDuring=allmeanevokedLocomotion_responsiveIMMOB-allmeanevokedImmobility_responsiveIMMOB;
figure;

plot(dPre,dDuring,'ko','markerfacecolor','m','markersize',14)
[rr pp]=corrcoef(dPre,dDuring);
pp1=polyfit(dPre,dDuring,1);
hold on
xlabel('Locomotion influence on ongoing (prestim) activity');
ylabel('Locomotion influence on evoked (stim) activity');

plot([-2:2],polyval(pp1,[-2:2]),'r--')
%axis([-1.5 2 -1.5 2])
plot([-2 2],[-2 2],'k')
p5=signrank(dPre,dDuring)
title(['R= ' num2str(rr(1,2)) ' Pcorr= ' num2str(pp(1,2)) ' Pdiff= ' num2str(p5)])

figure('Position',[50,50,300,600]);
boxplot([dPre-dDuring],'Symbol','');
ylim([-0.3 0.3])
hold on
plot([-1 3],[0 0],'r--')


%% CONST vs NONE

alldffmat=[allmeanDuringImmobileCONST' allmeanDuringMoveCONST'];
alldffmatResp=[allmeanDuringImmobileCONST(allCellImmobileResponsivenessCONST==1)' allmeanDuringMoveCONST(allCellImmobileResponsivenessCONST==1)'];

figure('Position', [50, 50, 300, 500])
subplot(2,1,1)
bar(mean(alldffmat))
hold on
errorbar(mean(alldffmat),std(alldffmat)/sqrt(size(alldffmat,1)),'.k','linewidth',2)
ylabel('Mean dff')
set(gca,'xticklabel',{'BG sound, immobile','BG sound, moving'})
title('all cells')

subplot(2,1,2)
bar(mean(alldffmatResp))
hold on
errorbar(mean(alldffmatResp),std(alldffmatResp)/sqrt(size(alldffmatResp,1)),'.k','linewidth',2)
ylabel('Mean dff')
set(gca,'xticklabel',{'BG sound, immobile','BG sound, moving'})
title('responsive cells')
%%
signrank(allmeanDuringImmobileCONST, allmeanDuringMoveCONST)

signrank(allmeanDuringImmobileCONST(allCellImmobileResponsivenessCONST==1), allmeanDuringMoveCONST(allCellImmobileResponsivenessCONST==1))

%% boxplots, locomotion minus immob
figure('Position',[50,50,300,600]);
subplot(1,2,1)
boxplot([allmeanDuringMoveCONST-allmeanDuringImmobileCONST],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.2 0.2])
ylabel('move minus immobile')
title('all cells')

subplot(1,2,2)
boxplot([allmeanDuringMoveCONST(allCellImmobileResponsivenessCONST==1)-allmeanDuringImmobileCONST(allCellImmobileResponsivenessCONST==1)],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-0.2 0.2])
ylabel('move minus immobile')
title('significant cells')

signrank(allmeanDuringMoveCONST-allmeanDuringImmobileCONST,0)

signrank(allmeanDuringMoveCONST(allCellImmobileResponsivenessCONST==1)-allmeanDuringImmobileCONST(allCellImmobileResponsivenessCONST==1),0)




%% dff-speed correlation
figure('Position',[50,50,600,300]);
histogram(alldffspeedcorrs,[-1:0.1:1],'Normalization','probability');
xlabel('dff-speed correlation');ylabel('cell count')
xlim([-1 1])
skewness(alldffspeedcorrs)
%% comparison to electro
load spikingspeedcorr
figure('Position',[50,50,600,300]);
histogram(alldffspeedcorrs,[-1:0.1:1],'Normalization','probability');
hold on
histogram(spikingspeedcorr,[-1:0.1:1],'Normalization','probability');
xlabel('dff-speed correlation');ylabel('cell count')
xlim([-1 1])
%% Proportion of cells whose spontaneous activity is correlated with locomotion speed
ratePosMod=mean(alldffspeedcorrsSignificance<0.05 & alldffspeedcorrs>0);
ratePosNeg=mean(alldffspeedcorrsSignificance<0.05 & alldffspeedcorrs<0);
ratePodUn=mean(alldffspeedcorrsSignificance>0.05);

figure;
bar([ratePosMod ratePosNeg ratePodUn])
hold on;plot([0 4],[0.05 0.05],'r--')
ylim([0 1])



%% ensemble-level ranges of dff-speed correlations
rng(22)
figure;
bar(mean(alldffspeedcorrsRanges))
hold on
errorbar(mean(alldffspeedcorrsRanges), std(alldffspeedcorrsRanges)/sqrt(length(alldffspeedcorrsRanges)),'.k','linewidth',2);
%ylim([0 1])
ylabel('Range of dff-speed correlations')
plot(1+0.1*randn(1,length(alldffspeedcorrsRanges)),alldffspeedcorrsRanges,'ko','markerfacecolor','m')
%plot(alldffspeedcorrsRanges

%% ENSEMBLE ANALYSIS
notrun=1;
if notrun
    load dataAfterNewImagingRun041422

    
end
%%
AllImmobileStimDecodingGain=1-AllImmobileStimDecodingError;
AllLocomotionStimDecodingGain=1-AllLocomotionStimDecodingError;
figure;
plot(AllImmobileStimDecodingGain,AllLocomotionStimDecodingGain,'ko');
hold on;
plot([0,1],[0,1])
onlyMove=find(AllImmobileStimDecodingPVAL>0.05&AllLocomotionStimDecodingPVAL<=0.05);
onlyImmob=find(AllImmobileStimDecodingPVAL<=0.05&AllLocomotionStimDecodingPVAL>0.05);
bothsig=find(AllImmobileStimDecodingPVAL<=0.05&AllLocomotionStimDecodingPVAL<=0.05);

plot(AllImmobileStimDecodingGain(onlyMove),AllLocomotionStimDecodingGain(onlyMove),'ko','markerfacecolor','r');
plot(AllImmobileStimDecodingGain(onlyImmob),AllLocomotionStimDecodingGain(onlyImmob),'ko','markerfacecolor','g');
plot(AllImmobileStimDecodingGain(bothsig),AllLocomotionStimDecodingGain(bothsig),'go','markersize',8,'markerfacecolor','r','linewidth',2);

xlabel('Detection GAIN during immobility')
ylabel('Detection GAIN during locomotion')
title('Ensemble stim detection')




%% SPEED PREDICTION
allspeedpredictiondata=[allSpeedPredictionRs; allSpeedPredictionRsSHUF; allSpeedPredictionRsNo0; allSpeedPredictionRsNo0SHUF];
figure('Position',[50,50,300,600]);
bar(mean(allspeedpredictiondata'))
hold on
jitterMat=[ones(1,length(allspeedpredictiondata)); 2*ones(1,length(allspeedpredictiondata)); 3*ones(1,length(allspeedpredictiondata)); 4*ones(1,length(allspeedpredictiondata))];
jitterMat=jitterMat+0.1*randn(size(jitterMat));
plot(jitterMat,allspeedpredictiondata,'ko','markerfacecolor','k')

errorbar(mean(allspeedpredictiondata'), std(allspeedpredictiondata')/sqrt(size(allspeedpredictiondata,2)),'r.','linewidth',2)
ylim([-1 1])
ylabel('Speed prediction performance')

set(gca,'xticklabel',{'Real','Shuffled'})
%%
figure;
plot(predictionEnsembleSize,allSpeedPredictionRs,'ko','markerfacecolor','m')
ylim([0 1])
xlabel('Ensemble size')
ylabel('Speed prediction')
pp=polyfit(predictionEnsembleSize,allSpeedPredictionRs,1);
hold on

plot([0:100],polyval(pp,[0:100]),'r')
[rrr ppp]=corrcoef(predictionEnsembleSize,allSpeedPredictionRs)



%% State decoding (from pre Stim Window) against stim detection in locomotion
AllPreWinStateDecodingGain=1-AllPreWinStateDecodingError;

meanState=mean(AllPreWinStateDecodingGain);
stdState=std(AllPreWinStateDecodingGain);%;/sqrt(length((AllPreWinStateDecodingGain)));

meanStim=mean(AllLocomotionStimDecodingGain);
stdStim=std(AllLocomotionStimDecodingGain);%;/sqrt(length((AllLocomotionStimDecodingGain)));

figure('Position',[50,50,1100,300]);

subplot(1,3,1)
plot([0 1],[0.5 0.5],'r--')
hold on
plot([0.5 0.5],[0 1],'r--')

plot(AllPreWinStateDecodingGain,AllLocomotionStimDecodingGain,'ko');axis([0 1 0 1])
hold on
plot(AllPreWinStateDecodingGain(AllPreWinStateDecodingPVAL<0.05),AllLocomotionStimDecodingGain(AllPreWinStateDecodingPVAL<0.05),'ko','markerfacecolor',[0 0.8 1])
plot(AllPreWinStateDecodingGain(AllPreWinStateDecodingPVAL<0.05&AllLocomotionStimDecodingPVAL<0.05),AllLocomotionStimDecodingGain(AllPreWinStateDecodingPVAL<0.05&AllLocomotionStimDecodingPVAL<0.05),'o','color',[1 0.2 0],'linewidth',2)

hold on
plot([meanState-stdState meanState+stdState],[meanStim meanStim],'k','linewidth',3)
plot([meanState meanState],[meanStim-stdStim meanStim+stdStim],'k','linewidth',3)

xlabel('State decoding from pre stim window')
ylabel('Stim dection in locomotion')

signrank(AllPreWinStateDecodingGain,0.5)
signrank(AllLocomotionStimDecodingGain,0.5)

%% ENSEMBLE Decoding of State AND Stim
AllStimAndStateDecodingGain=(1-AllStimAndStateDecodingError);
figure('Position',[50,50,1100,300]);

subplot(1,3,1)
bar(2,mean(AllStimAndStateDecodingGain));
hold on
errorbar(2,mean(AllStimAndStateDecodingGain),std(AllStimAndStateDecodingGain)/sqrt(length(AllStimAndStateDecodingGain)),'k','linewidth',2);

hold on
rng(7)
rr=0.1*randn(1,length(AllStimAndStateDecodingGain));
plot(2+rr,AllStimAndStateDecodingGain,'ko');hold on;
plot(2+rr(AllStimAndStateDecodingPVAL<0.05),AllStimAndStateDecodingGain(AllStimAndStateDecodingPVAL<0.05),'ko','markerface','m');
plot([1 3],[0.25 0.25],'r--')
ylim([0 1])

signrank(AllStimAndStateDecodingGain,0.25)






