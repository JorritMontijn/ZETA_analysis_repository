function UpdateWaitbar(~)
	global intWaitbarCounter;
	global intWaitbarTotal;
	global ptrWaitbarHandle;
	if ~exist('ptrWaitbarHandle','var') || isempty(ptrWaitbarHandle)
		ptrWaitbarHandle = waitbar(0, 'Please wait ...');
	end
	if ~exist('intWaitbarCounter','var') || isempty(intWaitbarCounter)
		intWaitbarCounter = 0;
	end
	waitbar(intWaitbarCounter/intWaitbarTotal, ptrWaitbarHandle,sprintf('Please wait ... Finished %d/%d',intWaitbarCounter,intWaitbarTotal));
	intWaitbarCounter = intWaitbarCounter + 1;
	if intWaitbarCounter == intWaitbarTotal
		delete(ptrWaitbarHandle);
	end
end