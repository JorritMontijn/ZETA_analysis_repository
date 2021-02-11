function sMIMI = fitMIMI(intCoeffsL1,intCoeffsG1,vecSpikeTimes,vecStimOnTime,dblTrialDur,boolVerbose)
	%fitMIMI Fits multiplicative inhomogeneous markov interval model
	%    sMIMI = fitMIMI(intCoeffsL1,intCoeffsG1,vecSpikeTimes,vecStimOnTime,dblTrialDur)
	%Inputs:
	%intCoeffsL1 = number of coefficients for linear poisson part
	%intCoeffsG1 = number of coefficients for non-linear part
	%vecSpikeTimes = spike time vector
	%vecStimOnTime = stimulus onset times
	%dblTrialDur = trial duration to fit 
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%% get spike times in trials
	vecSpikeTimes = sort(vecSpikeTimes);
	vecStimOnTime = sort(vecStimOnTime);
	dblStartT = max([vecSpikeTimes(1) vecStimOnTime(1)-dblTrialDur]);
	dblStopT = vecStimOnTime(end)+dblTrialDur*2;
	dblStepSize = 1/1000;
	vecAllT = dblStartT:dblStepSize:dblStopT;
	vecAll_ISI = nan(size(vecAllT));
	for intSpikeIdx=2:numel(vecSpikeTimes)
		dblPrevSpikeT = vecSpikeTimes(intSpikeIdx-1);
		dblSpikeT = vecSpikeTimes(intSpikeIdx);
		if dblSpikeT < dblStartT,continue;end
		if dblPrevSpikeT > dblStopT,break;end
		vecIdx = vecAllT < dblSpikeT & vecAllT >= dblPrevSpikeT;
		vecAssign = 0:dblStepSize:(dblSpikeT-dblPrevSpikeT);
		if dblSpikeT > dblStopT
			vecAll_ISI(vecIdx) = vecAssign(1:sum(vecIdx));
		else
			vecAll_ISI(vecIdx) = vecAssign((end+1-sum(vecIdx)):end);
		end
	end
	
	%assign to trials
	intTrials=numel(vecStimOnTime);
	dblBinSize = 1/1000;
	intBins = ceil(dblTrialDur/dblBinSize);
	matISI = zeros(intTrials,intBins);
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,vecStimOnTime,dblTrialDur);
	cellSpikePerTrial = cell(1,intTrials);
	matSpikeBins1ms = false(intTrials,intBins);
	for intTrial=1:intTrials
		cellSpikePerTrial{intTrial} = vecTimePerSpike(vecTrialPerSpike==intTrial);
		vecSpikeBins = ceil(vecTimePerSpike(vecTrialPerSpike==intTrial)*1000);
		vecSpikeBins(vecSpikeBins>intBins)=[];
		matSpikeBins1ms(intTrial,vecSpikeBins) = true;
		intStartISI = find(vecAllT>vecStimOnTime(intTrial),1);
		if isempty(intStartISI),continue;end
		matISI(intTrial,:) = vecAll_ISI(intStartISI:(intStartISI+intBins-1));
	end
	
	%% assign to globals
	global gMatISI;
	gMatISI = matISI;
	global gHyperParams_mIMI;
	gHyperParams_mIMI = [intCoeffsL1 intCoeffsG1];
	
	%% build vars
	vecX = (1:intBins)*dblBinSize;
	vecY = mean(matSpikeBins1ms/dblBinSize,1);
	vecY_sem = std(matSpikeBins1ms/dblBinSize,[],1)/sqrt(intTrials);
	vecCoeffsL1 = mean(vecY)*ones(1,intCoeffsL1);
	vecCoeffsG1 = zeros(1,intCoeffsG1);
	vecCoeffs0 = cat(2,vecCoeffsL1,vecCoeffsG1);
	
	sOptions = struct;
	if boolVerbose
		sOptions.Display = 'iter-detailed';
	else
		sOptions.Display = 'off';
	end
	
	%% fit
	[vecFitCoeffs,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@mIMI,vecCoeffs0,vecX,vecY,-10000*ones(size(vecCoeffs0)),10000*ones(size(vecCoeffs0)),sOptions);
	vecFitY=mIMI(vecFitCoeffs,vecX);
	%ci = nlparci(vecFitCoeffs,residual,'jacobian',jacobian);
	
	%% assign output
	sMIMI.vecX = vecX;
	sMIMI.vecY = vecY;
	sMIMI.vecFitY = vecFitY;
	sMIMI.FitCoeffs = vecFitCoeffs;
	sMIMI.resnorm = resnorm;
	sMIMI.residual = residual;
	sMIMI.exitflag = exitflag;
	sMIMI.output = output;
	sMIMI.lambda = lambda;
	sMIMI.jacobian = jacobian;
	
end

