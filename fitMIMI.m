function sMIMI = fitMIMI(intCoeffsL1,intCoeffsG1,vecSpikeTimes,vecStimOnTime,dblUseMaxDur,boolVerbose,vecCoeffs0)
	%fitMIMI Fits multiplicative inhomogeneous markov interval model
	%    sMIMI = fitMIMI(intCoeffsL1,intCoeffsG1,vecSpikeTimes,vecStimOnTime,dblTrialDur,boolVerbose,vecCoeffs0)
	%Inputs:
	%intCoeffsL1 = number of coefficients for linear poisson part
	%intCoeffsG1 = number of coefficients for non-linear part
	%vecSpikeTimes = spike time vector
	%vecStimOnTime = stimulus onset times
	%dblUseMaxDur = trial duration to fit 
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%% get spike times in trials
	vecSpikeTimes = sort(vecSpikeTimes);
	vecStimOnTime = sort(vecStimOnTime);
	dblStartT = max([vecSpikeTimes(1) vecStimOnTime(1)-dblUseMaxDur]);
	dblStopT = vecStimOnTime(end)+dblUseMaxDur*2;
	vecSpikeTimes(vecSpikeTimes < dblStartT | vecSpikeTimes > dblStopT) = [];
	
	%% check spike #
	sMIMI = [];
	if numel(vecSpikeTimes) < 3 || numel(vecStimOnTime) < 3
		return
	end
	
	%% build matrix
	dblStepSize = 1/1000;
	vecAllT = dblStartT:dblStepSize:dblStopT;
	vecTimeSinceSpike = nan(size(vecAllT));
	for intSpikeIdx=2:numel(vecSpikeTimes)
		dblPrevSpikeT = vecSpikeTimes(intSpikeIdx-1);
		dblSpikeT = vecSpikeTimes(intSpikeIdx);
		if dblSpikeT < dblStartT,continue;end
		if dblPrevSpikeT > dblStopT,break;end
		vecIdx = vecAllT < dblSpikeT & vecAllT >= dblPrevSpikeT;
		vecAssign = 0:dblStepSize:(dblSpikeT-dblPrevSpikeT);
		if dblSpikeT > dblStopT
			vecTimeSinceSpike(vecIdx) = vecAssign(1:sum(vecIdx));
		else
			vecTimeSinceSpike(vecIdx) = vecAssign((end+1-sum(vecIdx)):end);
		end
	end
	
	%assign to trials
	intTrials=numel(vecStimOnTime);
	intBins = ceil(dblUseMaxDur/dblStepSize);
	if intBins < (intCoeffsL1 + intCoeffsG1) %more coefficients than bins
		return
	end
	matTimeSinceSpike = zeros(intTrials,intBins);
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,vecStimOnTime,dblUseMaxDur);
	cellSpikePerTrial = cell(1,intTrials);
	matSpikeBins1ms = false(intTrials,intBins);
	for intTrial=1:intTrials
		cellSpikePerTrial{intTrial} = vecTimePerSpike(vecTrialPerSpike==intTrial);
		vecSpikeBins = ceil(vecTimePerSpike(vecTrialPerSpike==intTrial)*1000);
		vecSpikeBins(vecSpikeBins>intBins)=[];
		matSpikeBins1ms(intTrial,vecSpikeBins) = true;
		intStartISI = find(vecAllT>vecStimOnTime(intTrial),1);
		if isempty(intStartISI),continue;end
		matTimeSinceSpike(intTrial,:) = vecTimeSinceSpike(intStartISI:(intStartISI+intBins-1));
	end
	
	%% assign to globals
	global gMatISI;
	gMatISI = matTimeSinceSpike;
	global gHyperParams_mIMI;
	gHyperParams_mIMI = [intCoeffsL1 intCoeffsG1];
	
	%% build vars
	vecX = (1:intBins)*dblStepSize;
	vecY = nanmean(matSpikeBins1ms/dblStepSize,1);
	vecY_sem = nanstd(matSpikeBins1ms/dblStepSize,[],1)/sqrt(intTrials);
	vecCoeffsL1 = nanmean(vecY)*ones(1,intCoeffsL1);
	vecCoeffsG1 = zeros(1,intCoeffsG1);
	if ~exist('vecCoeffs0','var') || numel(vecCoeffs0) ~= (intCoeffsL1 + intCoeffsG1)
		vecCoeffs0 = cat(2,vecCoeffsL1,vecCoeffsG1);
	end
	
	sOptions = struct;
	if boolVerbose
		sOptions.Display = 'iter-detailed';
	else
		sOptions.Display = 'off';
	end
	
	%% fit
	[vecFitCoeffs,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@mIMI,vecCoeffs0,vecX,vecY,-10000*ones(size(vecCoeffs0)),10000*ones(size(vecCoeffs0)),sOptions);
	vecFitY=mIMI(vecFitCoeffs,vecX);
	
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

