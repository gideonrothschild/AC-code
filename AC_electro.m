
gatherandsave=0;
combineddatafile='E:\GideonData\combinedSpeedModData.mat';
load(combineddatafile);
%% cell numbers
rng(77)
allindsX=[];
celltypes=[];
for xx=1:length(combinedSBTrackData)
    if any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKmat'))
        
        allindsX=[allindsX;combinedSBTrackData{xx}.cellind];
        celltypes=[celltypes combinedSBTrackData{xx}.celltype];
    end
end
totalnumcellsontrack=length(celltypes);
totalnumEXCcellsontrack=sum(celltypes==1);
%% Determining sound responsiveness on track
speedthresh=4;

tracksoundresponsivenessNEW=[];
tracksoundresponsivenessNEWIMMOB=[];

for xx=1:length(combinedSBTrackData)
    
    if any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKmat'))
        
        
        allraster=combinedSBTrackData{xx}.TRACKmat;
        curprerespamps=sum(allraster(:,101:550),2);
        currespamps=sum(allraster(:,551:1000),2);
        
        [hResp pResp]=ttest(currespamps,curprerespamps,'tail','right');
        combinedSBTrackData{xx}.TRACK_NEW_SOUNDRESP_P=pResp;
        tracksoundresponsivenessNEW=[tracksoundresponsivenessNEW pResp<0.05];
        
        sp=combinedSBTrackData{xx}.TRACKspeeds;
        immobraster=allraster(sp<speedthresh,:);
        curprerespampsIMMOB=sum(immobraster(:,101:550),2);
        currespampsIMMOB=sum(immobraster(:,551:1000),2);
        [hRespIMMOB pRespIMMOB]=ttest(currespampsIMMOB,curprerespampsIMMOB,'tail','right');
        combinedSBTrackData{xx}.TRACK_NEW_SOUNDRESP_P_IMMOB=pRespIMMOB;
        tracksoundresponsivenessNEWIMMOB=[tracksoundresponsivenessNEWIMMOB pRespIMMOB<0.05];
      
        g=1;
    end
end

numOfEXCCellsResponsiveToTargetSound=sum(tracksoundresponsivenessNEW==1&celltypes==1);
rateOfEXCCellsResponsiveToTargetSound=numOfEXCCellsResponsiveToTargetSound/totalnumEXCcellsontrack;

numOfEXCCellsResponsiveToTargetSoundIMMOB=sum(tracksoundresponsivenessNEWIMMOB==1&celltypes==1);
rateOfEXCCellsResponsiveToTargetSound=numOfEXCCellsResponsiveToTargetSoundIMMOB/totalnumEXCcellsontrack;
%% Sound-evoked response in immobility vs. locomotion

speedthresh=4;
minnumtrials=10;
allcurinds1=[];
soundstart=551;
prestart=51;

soundtimewin=450;
soundend=soundstart+soundtimewin-1;
preend=prestart+soundtimewin-1;
allrespspeedcorR=[];
allrespspeedcorP=[];
allsoundresponsive=[];
allMove=[];
allImmob=[];
allsigs=[];
allsigsP=[];
allMoveUN=[];
allImmobUN=[];
allsigsUN=[];
allsigsPUN=[];
allMovePRE=[];
allImmobPRE=[];
allsigsPRE=[];
allsigsPPRE=[];
allresponsiveness=[];
allresponsivenessIMMOB=[];

allinds1=[];
for xx=1:length(combinedSBTrackData)
    if any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKmat'))
        if combinedSBTrackData{xx}.celltype==1
            
            
            
            currespmat=combinedSBTrackData{xx}.TRACKmat;
            trackaudresp=nansum(currespmat(:,soundstart:soundend),2)-nansum(currespmat(:,prestart:preend),2);
            trackaudrespUN=nansum(currespmat(:,soundstart:soundend),2);
            trackaudrespPRE=nansum(currespmat(:,prestart:preend),2);
            
            sp=combinedSBTrackData{xx}.TRACKspeeds;
            
            immobResps=trackaudresp(sp<speedthresh);
            moveResps=trackaudresp(sp>=speedthresh);
            [h p]=ttest2(immobResps,moveResps);
            
            immobRespsUN=trackaudrespUN(sp<speedthresh);
            moveRespsUN=trackaudrespUN(sp>=speedthresh);
            [hUN pUN]=ttest2(immobRespsUN,moveRespsUN);
            
            immobRespsPRE=trackaudrespPRE(sp<speedthresh);
            moveRespsPRE=trackaudrespPRE(sp>=speedthresh);
            [hPRE pPRE]=ttest2(immobRespsPRE,moveRespsPRE);
            
            if length(moveResps)>minnumtrials & length(immobResps)>minnumtrials
                
                if combinedSBTrackData{xx}.TRACK_NEW_SOUNDRESP_P<0.05
                    allresponsiveness=[allresponsiveness 1];
                else
                    allresponsiveness=[allresponsiveness 0];
                end
                if combinedSBTrackData{xx}.TRACK_NEW_SOUNDRESP_P_IMMOB<0.05
                    allresponsivenessIMMOB=[allresponsivenessIMMOB 1];
                else
                    allresponsivenessIMMOB=[allresponsivenessIMMOB 0];
                end
                
                
                allMove=[allMove nanmean(moveResps)];
                allImmob=[allImmob nanmean(immobResps)];
                allsigs=[allsigs h];
                allsigsP=[allsigsP p];
                allcurinds1=[allcurinds1; combinedSBTrackData{xx}.cellind];
                
                if h==1 & nanmean(immobResps)>nanmean(moveResps)
                    combinedSBTrackData{xx}.pref='immob';
                elseif h==1 & nanmean(immobResps)<nanmean(moveResps)
                    combinedSBTrackData{xx}.pref='move';
                elseif h==0
                    combinedSBTrackData{xx}.pref='no pref';
                end
                
                
                allMoveUN=[allMoveUN nanmean(moveRespsUN)];
                allImmobUN=[allImmobUN nanmean(immobRespsUN)];
                allsigsUN=[allsigsUN hUN];
                allsigsPUN=[allsigsPUN pUN];
                
                
                allMovePRE=[allMovePRE nanmean(moveRespsPRE)];
                allImmobPRE=[allImmobPRE nanmean(immobRespsPRE)];
                allsigsPRE=[allsigsPRE hPRE];
                allsigsPPRE=[allsigsPPRE pPRE];
                
            end
        end
    end
end


%% NORMALIZED sound-evoked responses in immobility vs locomotion (Excitatory immobility-sound-responsive cells only)

allImmob_responsive=allImmob(allresponsivenessIMMOB==1);
allMove_responsive=allMove(allresponsivenessIMMOB==1);
allsigs_responsive=allsigs(allresponsivenessIMMOB==1);

figure('Position',[50,50,1100,300]);

subplot(1,3,1)
plot(allImmob_responsive,allMove_responsive,'ko');hold on;
plot(allImmob_responsive(allsigs_responsive==0),allMove_responsive(allsigs_responsive==0),'ko','markerfacecolor','b');hold on;
plot(allImmob_responsive(allsigs_responsive==1&allImmob_responsive<allMove_responsive),allMove_responsive(allsigs_responsive==1&allImmob_responsive<allMove_responsive),'ko','markerfacecolor','g');hold on;
plot(allImmob_responsive(allsigs_responsive==1&allImmob_responsive>allMove_responsive),allMove_responsive(allsigs_responsive==1&allImmob_responsive>allMove_responsive),'ko','markerfacecolor','r');hold on;

plot([-2 10],[-2 10],'k--')
xlabel('Sound-evoked response magnitude during immobility')
ylabel('Sound-evoked response magnitude during locomotion')
title('Excitatory sound-responsive cells only')
axis([-2 10 -2 10])

subplot(1,3,2)
b=bar([mean(allImmob_responsive) mean(allMove_responsive)]);
b.FaceColor = 'flat';
b.CData(1,:)=[1 0 0];
b.CData(2,:)=[0 1 0];

hold on
errorbar([mean(allImmob_responsive) mean(allMove_responsive)],[std(allImmob_responsive)/sqrt(length(allImmob_responsive)) std(allMove_responsive)/sqrt(length(allMove_responsive))],'k.','linewidth',2);
ylabel('Sound-evoked response (spikes)')
set(gca,'xticklabel',{'Immobility','Locomotion'})
p2=signrank(allImmob_responsive,allMove_responsive)
title(['P= ' num2str(p2)])

subplot(1,3,3)
dd=allMove_responsive-allImmob_responsive;
boxplot([dd],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylim([-3 3])




%% gather all psths of immobility-responsive neurons, immobility and locomotion


bb=fir1(51,0.1);
curxaxis=-549:550;
b=fir1(51,0.05);
f1=ones(3,10)/30;
allimmobilityPSTHs=[];
alllocomotionPSTHs=[];
tmpinds=[];
for xx=1:length(combinedSBTrackData)
    
    if any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKmat'))
        if combinedSBTrackData{xx}.celltype==1 & combinedSBTrackData{xx}.TRACK_NEW_SOUNDRESP_P_IMMOB<0.05
            
            
            sp=combinedSBTrackData{xx}.TRACKspeeds;
            %[sps spi]=sort(sp);
            immobilitytrials=find(sp<=4);
            locomotiontrials=find(sp>4);
            if length(immobilitytrials)>minnumtrials & length(locomotiontrials)>minnumtrials
                allraster=combinedSBTrackData{xx}.TRACKmat;
                immobilityraster=allraster(immobilitytrials,:);
                locomotionraster=allraster(locomotiontrials,:);
                allimmobilityPSTHs=[allimmobilityPSTHs; filtfilt(b,1,nanmean(immobilityraster))];
                alllocomotionPSTHs=[alllocomotionPSTHs; filtfilt(b,1,nanmean(locomotionraster))];
                tmpinds=[tmpinds;combinedSBTrackData{xx}.cellind];
            end
            % end
        end
    end
end
%
curxaxis=-549:550;

figure;
shadedErrorBar(curxaxis,nanmean(allimmobilityPSTHs),nanstd(allimmobilityPSTHs)/sqrt(size(allimmobilityPSTHs,1)),'lineprops','r')
hold on
shadedErrorBar(curxaxis,nanmean(alllocomotionPSTHs),nanstd(alllocomotionPSTHs)/sqrt(size(alllocomotionPSTHs,1)),'lineprops','g')
xlim([-300 500])
%% SPONT (PRE-Sound) in immobility vs locomotion (Excitatory immobility-sound-responsive cells only)

allImmobPRE_responsive=allImmobPRE(allresponsivenessIMMOB==1);
allMovePRE_responsive=allMovePRE(allresponsivenessIMMOB==1);
allsigsPRE_responsive=allsigsPRE(allresponsivenessIMMOB==1);

figure('Position',[50,50,1100,300]);
subplot(1,3,1)
plot(allImmobPRE_responsive,allMovePRE_responsive,'ko');hold on;
plot(allImmobPRE_responsive(allsigsPRE_responsive==0),allMovePRE_responsive(allsigsPRE_responsive==0),'ko','markerfacecolor','b');hold on;
plot(allImmobPRE_responsive(allsigsPRE_responsive==1&allImmobPRE_responsive<allMovePRE_responsive),allMovePRE_responsive(allsigsPRE_responsive==1&allImmobPRE_responsive<allMovePRE_responsive),'ko','markerfacecolor','g');hold on;
plot(allImmobPRE_responsive(allsigsPRE_responsive==1&allImmobPRE_responsive>allMovePRE_responsive),allMovePRE_responsive(allsigsPRE_responsive==1&allImmobPRE_responsive>allMovePRE_responsive),'ko','markerfacecolor','r');hold on;

plot([-1 8],[-1 8],'k--')
axis([-1 8 -1 8])
xlabel('PRE response during immobility')
ylabel('PRE response during locomotion')
title('Excitatory sound-responsive cells only')
subplot(1,3,2)

b=bar([mean(allImmobPRE_responsive) mean(allMovePRE_responsive)]);
b.FaceColor = 'flat';
b.CData(1,:)=[1 0 0];
b.CData(2,:)=[0 1 0];

hold on
errorbar([mean(allImmobPRE_responsive) mean(allMovePRE_responsive)],[std(allImmobPRE_responsive)/sqrt(length(allImmobPRE_responsive)) std(allMovePRE_responsive)/sqrt(length(allMovePRE_responsive))],'k.','linewidth',2);
ylabel('PREnormalized soPREd-evoked response (spikes)')
set(gca,'xticklabel',{'Immobility','Locomotion'})
p2=signrank(allImmobPRE_responsive,allMovePRE_responsive);
title(['P= ' num2str(p2)])
subplot(1,3,3)
boxplot([allMovePRE_responsive-allImmobPRE_responsive],'symbol','');
hold on
plot([-1 3],[0 0],'r--')

%% %% Change in UNNORMALIZED evoked response against change in spont
allImmobUN_responsive=allImmobUN(allresponsivenessIMMOB==1);
allMoveUN_responsive=allMoveUN(allresponsivenessIMMOB==1);
allsigsUN_responsive=allsigsUN(allresponsivenessIMMOB==1);
figure;

plot(allMovePRE_responsive-allImmobPRE_responsive,allMoveUN_responsive-allImmobUN_responsive,'ko','markerfacecolor','m')
xlabel('Locomotion-modulation  of spont')
ylabel('Locomotion-modulation of UNNORMALIZED sound-evoked')
[rr pp]=corrcoef(allMovePRE_responsive-allImmobPRE_responsive,allMoveUN_responsive-allImmobUN_responsive);
title(['R= ' num2str(rr(1,2)) ' P= ' num2str(pp(1,2))])
pp1=polyfit(allMovePRE_responsive-allImmobPRE_responsive,allMoveUN_responsive-allImmobUN_responsive,1);

hold on
plot([-6:9],polyval(pp1,[-6:9]),'r')

axis([-6 9 -6 9])

%% locomotion influence on pre-stim window against locomotion influence on stim window
dPre=allMovePRE_responsive-allImmobPRE_responsive;
dDuring=allMoveUN_responsive-allImmobUN_responsive;

figure('Position',[50,50,1100,300]);

subplot(1,3,1)

dRespImmob=allImmobUN_responsive-allImmobPRE_responsive;
dRespMove=allMoveUN_responsive-allMovePRE_responsive;

bar([mean(dRespImmob) mean(dRespMove)]);
hold on
errorbar([mean(dRespImmob) mean(dRespMove)],[std(dRespImmob)/sqrt(length(dRespImmob)) std(dRespMove)/sqrt(length(dRespMove))] ,'.k','linewidth',2);
title('Effect of sound on neural activity')
set(gca,'xticklabel',{'Immobility','Locomotion'})


subplot(1,3,2)
plot(dPre,dDuring,'ko','markerfacecolor','m')
[rr pp]=corrcoef(dPre,dDuring);
pp1=polyfit(dPre,dDuring,1);
hold on
xlabel('Locomotion influence on ongoing (prestim) activity');
ylabel('Locomotion influence on evoked (stim) activity');
plot([-6:9],polyval(pp1,[-6:9]),'r')
axis([-6 9 -6 9])


plot([-9 9],[-9 9],'k')
p5=signrank(dPre,dDuring)
title(['R= ' num2str(rr(1,2)) ' Pcorr= ' num2str(pp(1,2)) ' Pdiff= ' num2str(p5)])

subplot(1,3,3)
difMovEffectOnSpontMinusEvoked=dPre-dDuring;
boxplot([difMovEffectOnSpontMinusEvoked],'Symbol','');
hold on
plot([-1 3],[0 0],'r--')
ylabel('Pre minus during')
ylim([-3 3])

title('Effect of locomotion on neural activity')
set(gca,'xticklabel',{'pre-stim window','stim window'})


%% locomotion influence on pre-stim window against locomotion influence on S/N
dPre=allMovePRE_responsive-allImmobPRE_responsive;
dSN=allMove_responsive-allImmob_responsive;

figure('Position',[50,50,1100,300]);
plot(dPre,dSN,'ko','markerfacecolor','m')



%% Sound-evoked response in immobility against change in spont
% There is a trend but no significant correlation between the magnitude of sound-evoked
% response (in immobility) and locomotion-modulation of spontaneous activity
% i.e., it is not that the strongly responsive neurons all decrease spont
% during locomotion
figure;
plot(allMovePRE_responsive-allImmobPRE_responsive,allImmob_responsive,'ko','markerfacecolor','m')
xlabel('Locomotion-modulation  of spont')
ylabel('NORMALIZED sound-evoked response in immobility')
[rr pp]=corrcoef(allMovePRE_responsive-allImmobPRE_responsive,allImmob_responsive);
%% SPONT (PRE-Sound) in immobility vs locomotion (all excitatory cells)
figure('Position',[50,50,1100,300]);
subplot(1,3,1)

plot(allImmobPRE,allMovePRE,'ko');hold on;
plot(allImmobPRE(allsigsPRE==0),allMovePRE(allsigsPRE==0),'ko','markerfacecolor','b');hold on;
plot(allImmobPRE(allsigsPRE==1&allImmobPRE<allMovePRE),allMovePRE(allsigsPRE==1&allImmobPRE<allMovePRE),'ko','markerfacecolor','g');hold on;
plot(allImmobPRE(allsigsPRE==1&allImmobPRE>allMovePRE),allMovePRE(allsigsPRE==1&allImmobPRE>allMovePRE),'ko','markerfacecolor','r');hold on;

plot([-1 9],[-1 9],'k--')
axis([-1 9 -1 9])
xlabel('PRE response during immobility')
ylabel('PRE response during locomotion')
title('All excitatory cells')

subplot(1,3,2)
b=bar([mean(allImmobPRE) mean(allMovePRE)]);
b.FaceColor = 'flat';
b.CData(1,:)=[1 0 0];
b.CData(2,:)=[0 1 0];

hold on
errorbar([mean(allImmobPRE) mean(allMovePRE)],[std(allImmobPRE)/sqrt(length(allImmobPRE)) std(allMovePRE)/sqrt(length(allMovePRE))],'k.','linewidth',2);
ylabel('PREnormalized soPREd-evoked response (spikes)')
set(gca,'xticklabel',{'Immobility','Locomotion'})
p2=signrank(allImmobPRE,allMovePRE);
title(['P= ' num2str(p2)])

subplot(1,3,3)
boxplot([allImmobPRE' allMovePRE'],'symbol','o','outliersize',4);xlim([-1 4])



%% Spiking-speed corr

spikingspeedcorr=[];
celltypes=[];
alltrackaudresp=[];
allinds1=[];
for xx=1:length(combinedSBTrackData)
    if any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKALLBINNEDspiking'))& any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKmat'))
        if combinedSBTrackData{xx}.celltype==1
            curspiking=combinedSBTrackData{xx}.TRACKALLBINNEDspiking;
            curspeeds=combinedSBTrackData{xx}.TRACKALLBINNEDspeeds;
            [rd pd]=corrcoef(curspiking,curspeeds);
            spikingspeedcorr=[spikingspeedcorr rd(1,2)];
            allinds1=[allinds1;combinedSBTrackData{xx}.cellind];
        end
    end
end

figure;
histogram(spikingspeedcorr,[-1:0.1:1],'Normalization','probability')
xlabel('dff-speed correlation');ylabel('cell count')

xlim([-1 1])
skewness(spikingspeedcorr)
%% ENSEMBLE stim-detection during immobility and locomotion
torun=1;
if torun
    rng(77)
    minnumtrials=20;
    curindlist=[];
    allcelldata={};
    allEnsembleMoveLoss=[];
    allEnsembleImmobLoss=[];
    allEnsembleStateLoss=[];
    allEnsembleStateSig=[];
    allEnsembleStateLossShuf=[];
    allEnsembleMoveSig=[];
    allEnsembleImmobSig=[];
    
    allEnsembleStateLossDuringStim=[];
    allEnsembleStateStimLossDuringStim=[];
    allEnsembleStateLossShufDuringStim=[];
    allEnsembleStateStimLossShufDuringStim=[];
    allEnsembleStateSigDuringStim=[];
    allEnsembleStateStimSigDuringStim=[];
    allEnsembleMoveLossShuf=[];
    allEnsembleImmobLossShuf=[];
    allensembleinds={};
    x=1;
   
    
    for xx=1:length(combinedSBTrackData)
        if any(strcmp(fieldnames(combinedSBTrackData{xx}),'TRACKmat')) & combinedSBTrackData{xx}.celltype==1
            curindlist=[curindlist;combinedSBTrackData{xx}.cellind];
            currespmat=combinedSBTrackData{xx}.TRACKmat;
            trackaudresp=sum(currespmat(:,soundstart:soundend),2);
            trackaudpreresp=sum(currespmat(:,prestart:preend),2);
            
            sp=combinedSBTrackData{xx}.TRACKspeeds;
            allcelldata{x}.ind=combinedSBTrackData{xx}.cellind;
            allcelldata{x}.immobResps=trackaudresp(sp<speedthresh);
            allcelldata{x}.moveResps=trackaudresp(sp>=speedthresh);
            allcelldata{x}.immobPreResps=trackaudpreresp(sp<speedthresh);
            allcelldata{x}.movePreResps=trackaudpreresp(sp>=speedthresh);
            allcelldata{x}.TRACKmat=combinedSBTrackData{xx}.TRACKmat;
            x=x+1;
        end
    end
    allanimdays=unique(curindlist(:,1:2),'rows');
    for j=1:size(allanimdays,1)
        curanimday=allanimdays(j,:);
        curensembleimmobResps=[];
        curensembleimmobPreResps=[];
        curensemblemoveResps=[];
        curensemblemovePreResps=[];
        curensembleinds=find(ismember(curindlist(:,1:2),curanimday,'rows'));
     
        mxlength=0;
        for k=1:length(curensembleinds)
            mxlength=max(mxlength,length(allcelldata{curensembleinds(k)}.immobResps));
        end
        ensembleinds=[];
        for k=1:length(curensembleinds)
            % adding only cells that have the max num of responses
            if length(allcelldata{curensembleinds(k)}.immobResps)==mxlength
                curensembleimmobResps=[curensembleimmobResps allcelldata{curensembleinds(k)}.immobResps];
                curensembleimmobPreResps=[curensembleimmobPreResps allcelldata{curensembleinds(k)}.immobPreResps];
                curensemblemoveResps=[curensemblemoveResps allcelldata{curensembleinds(k)}.moveResps];
                curensemblemovePreResps=[curensemblemovePreResps allcelldata{curensembleinds(k)}.movePreResps];
                ensembleinds=[ensembleinds;allcelldata{curensembleinds(k)}.ind];
            end
        end
        if size(curensembleimmobResps,1)>minnumtrials & size(curensemblemoveResps,1)>minnumtrials
            allensembleinds{end+1}=ensembleinds;
            minlength1=min(size(curensembleimmobResps,1),size(curensemblemoveResps,1));
            
            curensembleimmobResps=curensembleimmobResps(1:minlength1,:);
            curensembleimmobPreResps=curensembleimmobPreResps(1:minlength1,:);
            curensemblemoveResps=curensemblemoveResps(1:minlength1,:);
            curensemblemovePreResps=curensemblemovePreResps(1:minlength1,:);
            
            
            
            immobiletags=[zeros(size(curensembleimmobResps,1),1);ones(size(curensembleimmobPreResps,1),1)];
            MdlImmob = fitcsvm([curensembleimmobResps; curensembleimmobPreResps],immobiletags);
            CVSVMModelImmob = crossval(MdlImmob);
            classLossImmob = kfoldLoss(CVSVMModelImmob);
            
            movetags=[zeros(size(curensemblemoveResps,1),1);ones(size(curensemblemovePreResps,1),1)];
            MdlMove = fitcsvm([curensemblemoveResps; curensemblemovePreResps],movetags);
            CVSVMModelMove = crossval(MdlMove);
            classLossMove = kfoldLoss(CVSVMModelMove);
            
            % decoding state by pre-stim window
            statetags=[zeros(size(curensembleimmobPreResps,1),1);ones(size(curensemblemovePreResps,1),1)];
            MdlState = fitcsvm([curensembleimmobPreResps; curensemblemovePreResps],statetags);
            CVSVMModelState = crossval(MdlState);
            classLossStateEnsemble = kfoldLoss(CVSVMModelState);
            
            % decoding state by stim window
            statetagsDuringStim=[zeros(size(curensembleimmobResps,1),1);ones(size(curensemblemoveResps,1),1)];
            MdlStateDuringStim = fitcsvm([curensembleimmobResps; curensemblemoveResps],statetagsDuringStim);
            CVSVMModelStateDuringStim = crossval(MdlStateDuringStim);
            classLossStateDuringStimEnsemble = kfoldLoss(CVSVMModelStateDuringStim);
            
            % decoding state and stim 
            statestimtagsDuringStim=[zeros(size(curensembleimmobPreResps,1),1);ones(size(curensemblemovePreResps,1),1);2*ones(size(curensembleimmobResps,1),1);3*ones(size(curensemblemoveResps,1),1)];
            MdlStateStimDuringStim = fitcdiscr([curensembleimmobPreResps;curensemblemovePreResps;curensembleimmobResps; curensemblemoveResps],statestimtagsDuringStim,'DiscrimType','pseudolinear');
            CVSVMModelStateStimDuringStim = crossval(MdlStateStimDuringStim);
            classLossStateStimDuringStimEnsemble = kfoldLoss(CVSVMModelStateStimDuringStim);
            
            
            numshufs1=200;
            allclassLossImmobShuf=[];
            allclassLossMoveShuf=[];
            allStateLossShufCurrEnsemble=[];
            allStateLossShufCurrDuringStimEnsemble=[];
            allStateStimLossShufCurrDuringStimEnsemble=[];
            
            for mm=1:numshufs1
                ['Day ' num2str(j) ' out of ' num2str(size(allanimdays,1)) ' shuf ' num2str(mm)]
                MdlImmobShuf = fitcsvm([curensembleimmobResps; curensembleimmobPreResps],immobiletags(randperm(length(immobiletags))));
                CVSVMModelImmobshuf = crossval(MdlImmobShuf);
                classLossImmobShuf = kfoldLoss(CVSVMModelImmobshuf);
                allclassLossImmobShuf=[allclassLossImmobShuf classLossImmobShuf];
                
                
                MdlMoveShuf = fitcsvm([curensemblemoveResps; curensemblemovePreResps],movetags(randperm(length(movetags))));
                CVSVMModelMoveshuf = crossval(MdlMoveShuf);
                classLossMoveShuf = kfoldLoss(CVSVMModelMoveshuf);
                allclassLossMoveShuf=[allclassLossMoveShuf classLossMoveShuf];
                
                
                MdlStateShuf = fitcsvm([curensembleimmobPreResps; curensemblemovePreResps],statetags(randperm(length(statetags))));
                CVSVMModelStateshuf = crossval(MdlStateShuf);
                classLossStateShuf = kfoldLoss(CVSVMModelStateshuf);
                allStateLossShufCurrEnsemble=[allStateLossShufCurrEnsemble classLossStateShuf];
                
                MdlStateShufDuringStim = fitcsvm([curensembleimmobResps; curensemblemoveResps],statetagsDuringStim(randperm(length(statetagsDuringStim))));
                CVSVMModelStateshufDuringStim = crossval(MdlStateShufDuringStim);
                classLossStateShufDuringStim = kfoldLoss(CVSVMModelStateshufDuringStim);
                allStateLossShufCurrDuringStimEnsemble=[allStateLossShufCurrDuringStimEnsemble classLossStateShufDuringStim];
                
                MdlStateStimDuringStimShuf = fitcdiscr([curensembleimmobPreResps;curensemblemovePreResps;curensembleimmobResps; curensemblemoveResps],statestimtagsDuringStim(randperm(length(statestimtagsDuringStim))),'DiscrimType','pseudolinear');
                CVSVMModelStateStimDuringStimShuf = crossval(MdlStateStimDuringStimShuf);
                classLossStateStimDuringStimShuf = kfoldLoss(CVSVMModelStateStimDuringStimShuf);
                allStateStimLossShufCurrDuringStimEnsemble=[allStateStimLossShufCurrDuringStimEnsemble classLossStateStimDuringStimShuf];
                
                
            end
            
            pImmob=mean(classLossImmob>allclassLossImmobShuf);
            pMove=mean(classLossMove>allclassLossMoveShuf);
            pState=mean(classLossStateEnsemble>allStateLossShufCurrEnsemble);
            pStateDuringStim=mean(classLossStateDuringStimEnsemble>allStateLossShufCurrDuringStimEnsemble);
            pStateStimDuringStim=mean(classLossStateStimDuringStimEnsemble>allStateStimLossShufCurrDuringStimEnsemble);
            
            
            allEnsembleMoveLoss=[allEnsembleMoveLoss classLossMove];
            allEnsembleImmobLoss=[allEnsembleImmobLoss classLossImmob];
            allEnsembleStateLoss=[allEnsembleStateLoss classLossStateEnsemble];
            allEnsembleStateLossDuringStim=[allEnsembleStateLossDuringStim classLossStateDuringStimEnsemble];
            allEnsembleStateStimLossDuringStim=[allEnsembleStateStimLossDuringStim classLossStateStimDuringStimEnsemble];
            
            allEnsembleMoveLossShuf=[allEnsembleMoveLossShuf mean(allclassLossMoveShuf)];
            allEnsembleImmobLossShuf=[allEnsembleImmobLossShuf mean(allclassLossImmobShuf)];
            allEnsembleStateLossShuf=[allEnsembleStateLossShuf mean(allStateLossShufCurrEnsemble)];
            allEnsembleStateLossShufDuringStim=[allEnsembleStateLossShufDuringStim mean(allStateLossShufCurrDuringStimEnsemble)];
            allEnsembleStateStimLossShufDuringStim=[allEnsembleStateStimLossShufDuringStim mean(allStateStimLossShufCurrDuringStimEnsemble)];
            
            allEnsembleMoveSig=[allEnsembleMoveSig pMove];
            allEnsembleImmobSig=[allEnsembleImmobSig pImmob];
            allEnsembleStateSig=[allEnsembleStateSig pState];
            allEnsembleStateSigDuringStim=[allEnsembleStateSigDuringStim pStateDuringStim];
            allEnsembleStateStimSigDuringStim=[allEnsembleStateStimSigDuringStim pStateStimDuringStim];
            
        end
    end
    
    figure;
    
    plot(allEnsembleImmobLoss,allEnsembleMoveLoss,'ko','markerfacecolor','m');hold on;
    plot([0 1],[0 1],'r--')
    xlabel('Tone-evoked detection error during immobility')
    ylabel('Tone-evoked detection error during locomotion')
else
    load dataAfterNewElectroRun041422
end

%
%% numbers of cells per ensemble in decoding
allnumcellsindecoding=[];
for i=1:length(allensembleinds)
    allnumcellsindecoding=[allnumcellsindecoding size(allensembleinds{i},1)];
end
mean(allnumcellsindecoding)
std(allnumcellsindecoding)/sqrt(length(allnumcellsindecoding))

%% stim coding in locomotion against state coding from prestim window
% === using this one
allEnsembleStateGain=1-allEnsembleStateLoss;
allEnsembleMoveGain=1-allEnsembleMoveLoss;


meanState=mean(allEnsembleStateGain);
semState=std(allEnsembleStateGain);%/sqrt(length(allEnsembleStateGain));

meanStim=mean(allEnsembleMoveGain);
semStim=std(allEnsembleMoveGain);%/sqrt(length(allEnsembleMoveGain));


figure('Position',[50,50,1100,300]);

subplot(1,3,1)
plot([0 1],[0.5 0.5],'r--')
hold on
plot([0.5 0.5],[0 1],'r--')
plot(allEnsembleStateGain,allEnsembleMoveGain,'ko')
hold on
plot(allEnsembleStateGain(allEnsembleStateSig<=0.05),allEnsembleMoveGain(allEnsembleStateSig<=0.05),'ko','markerfacecolor',[0 0.8 1])
plot(allEnsembleStateGain(allEnsembleMoveSig<=0.05),allEnsembleMoveGain(allEnsembleMoveSig<=0.05),'o','color',[1 0.2 0],'linewidth',2)
plot([meanState-semState meanState+semState],[meanStim meanStim],'k','linewidth',3)
plot([meanState meanState],[meanStim-semStim meanStim+semStim],'k','linewidth',3)

axis([0 1 0 1])
xlabel('State decoding from pre stim window')
ylabel('Stim detection in locomotion')


%% summary stats for above
popPsoundInMov=signrank(allEnsembleMoveGain-0.5)
popPstate=signrank(allEnsembleStateGain-0.5)

%% STATE AND STIM coding by ensembles based on STIM window
% USING THIS
figure('Position',[50,50,1500,400]);

subplot(1,3,1)

rng(77)
allEnsembleStateStimGainDuringStim=1-allEnsembleStateStimLossDuringStim;
bar(mean(allEnsembleStateStimGainDuringStim));
hold on

rr=0.3*rand(1,length(allEnsembleStateStimGainDuringStim));
plot(.85+rr,allEnsembleStateStimGainDuringStim,'ko');hold on;
plot([0 2],[0.25 0.25],'r--')
ylabel('State and stim decoding GAIN')
hold on
plot(.85+rr(allEnsembleStateStimSigDuringStim<0.05),allEnsembleStateStimGainDuringStim(allEnsembleStateStimSigDuringStim<0.05),'ko','markerfacecolor','m');hold on;
errorbar(mean(allEnsembleStateStimGainDuringStim),std(allEnsembleStateStimGainDuringStim)/sqrt(length(allEnsembleStateStimGainDuringStim)),'k','linewidth',2);

title('State AND stim decoding from ensembles, STIM win ')
ylim([0 1])
%% Population significance
signrank(allEnsembleStateStimGainDuringStim,0.25)
