function vecL = mIMI(vecCoeffs,vecX)
	%mIMI multiplicative inhomogeneous markov interval model
	%   vecY = mIMI(vecCoeffs,vecX)
	%
	%Specify the hyperparameters by creating a global gHyperParams_mIMI:
	%gHyperParams_mIMI(1) = intCoeffsL1
	%gHyperParams_mIMI(2) = intCoeffsG1
	%gHyperParams_mIMI(3) = intKnotsL1 (optional)
	%gHyperParams_mIMI(4) = intKnotsG1 (optional)
	%
	%Specify the training set
	%gMatSpikeBins1ms = [trials x bins] logical matrix
	%gMatISI = [trials x bins] ISI matrix
	%
	%Optionally specify a specific knot-sequence (not to be optimized) by
	%creating a global gHyperParams_mIMI_vecKnotsL1 containing the knots
	%Note that gHyperParams_mIMI(3) must be empty for
	%gHyperParams_mIMI_vecKnotsL1 to be used. Same for
	%%gHyperParams_mIMI_vecKnotsG1
	
	
	%% inputs
	global gMatISI
	global gHyperParams_mIMI;
	if isempty(gHyperParams_mIMI)
		error([mfilename ':HyperParamsMissing'],sprintf('Hyperparameters in gHyperParams_mIMI not initialized; please create a global first'));
	end
	intCoeffsL1 = gHyperParams_mIMI(1);
	intCoeffsG1 = gHyperParams_mIMI(2);
	if numel(gHyperParams_mIMI) > 2
		intKnotsL1 = gHyperParams_mIMI(3);
	else
		intKnotsL1 = [];
	end
	if numel(gHyperParams_mIMI) > 3
		intKnotsG1 = gHyperParams_mIMI(4);
	else
		intKnotsG1 = [];
	end
	
	%split inputs
	vecCoeffsL1 = vecCoeffs(1:intCoeffsL1);
	vecCoeffsG1 = vecCoeffs((intCoeffsL1+1):(intCoeffsL1+intCoeffsG1));
	
	%% retrieve knots or build uniform knot sequence if none supplied
	if isempty(intKnotsL1)
		global gHyperParams_mIMI_vecKnotsL1;
		if isempty(gHyperParams_mIMI_vecKnotsL1)
			intImposeOrder = 4;
			intKnotsL1 = numel(vecCoeffsL1)+intImposeOrder;
			vecKnotsL1 = min(vecX):(range(vecX)/(intKnotsL1-1-(intImposeOrder-1)*2)):max(vecX);
		else
			vecKnotsL1 = gHyperParams_mIMI_vecKnotsL1;
		end
	else
		vecKnotsL1 = vecCoeffsL1(sum(gHyperParams_mIMI(1:2)+1):end);
	end
	if isempty(intKnotsG1)
		global gHyperParams_mIMI_vecKnotsG1;
		if isempty(gHyperParams_mIMI_vecKnotsG1)
			intImposeOrder = 4;
			intKnotsG1 = numel(vecCoeffsG1)+intImposeOrder;
			dblRangeG1 = 0.5;
			vecKnotsG1 = 0:(dblRangeG1/(intKnotsG1-1-(intImposeOrder-1)*2)):dblRangeG1;
		else
			vecKnotsG1 = gHyperParams_mIMI_vecKnotsG1;
		end
	else
		vecKnotsG1 = vecCoeffsG1(sum(gHyperParams_mIMI(1:2)+1):end);
	end
	
	%% calc L1 component
	intCoeffsL1 = numel(vecCoeffsL1);
	vecKnotsAugL1 = augknt(vecKnotsL1, 4);
	intOrder = numel(vecKnotsAugL1) - intCoeffsL1; %cubic
	if intOrder ~= 4
		error([mfilename ':WrongOrder'],sprintf('%d coeffs with %d knots would result in order %d after knot augmentation',...
			intCoeffsL1,numel(vecKnotsL1),intOrder));
	end
	
	%build splines and evaluate Poisson part
	vecY_L1 = fnval(spmak(vecKnotsAugL1,vecCoeffsL1),vecX);
	
	%% calc G1 component; build splines and evaluate non-linear part
	intCoeffsG1 = numel(vecCoeffsG1);
	vecKnotsAugG1 = augknt(vecKnotsG1, 4);
	intOrder = numel(vecKnotsAugG1) - intCoeffsG1; %cubic
	if intOrder ~= 4
		error([mfilename ':WrongOrder'],sprintf('%d coeffs with %d knots would result in order %d after knot augmentation',...
			intCoeffsG1,numel(vecKnotsG1),intOrder));
	end
	intTrials = size(gMatISI,1);
	sSplines = spmak(vecKnotsAugG1,vecCoeffsG1);
	matY_G1 = fnval(sSplines,gMatISI);
	
	%% multiply by addition in log domain
	vecL = exp(nanmean(matY_G1,1) + vecY_L1);
end

