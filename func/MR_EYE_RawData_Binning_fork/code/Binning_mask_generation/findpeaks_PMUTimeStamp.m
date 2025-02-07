function [RRpeaks,RRperiods] = findpeaks_PMUTimeStamp(PMUTimeStamp,TimeStamp)
% FINDPEAKS_PMUTIMESTAMP
% Corrected version of the code to find the peaks in the PMUTimeStamp:
% it takes into account the "blind-spots" we have in the non-freerunning acquisitions
%
% INPUTS
% PMUTimeSTamp: temporal period from the last ECG trigger
% TimeStamp:    temporal period from the beginning of the acquisition
%
% OUTPUTS
% RRpeaks:      temporal position of the peaks
% RRperiods:    temporal interval between two consecutive peaks


% find the "lower peaks"
% findpeaks: Find local maxima of input -> find the minima of PMUTime
[ ~, idxPeaks ] = findpeaks( -PMUTimeStamp );
% disp('idxPeaks')
% disp(idxPeaks)
% subtract from each peak the time passed from the last TRUE cardiac trigger
% disp('TimeStamp( idxPeaks )')
% disp(TimeStamp( idxPeaks ))
% disp('PMUTimeStamp( idxPeaks )')
% disp(PMUTimeStamp( idxPeaks ))
% plot(PMUTimeStamp)
figure(1)
plot(idxPeaks, TimeStamp( idxPeaks ))
xlabel('idxPeaks');
ylabel('TimeStamp(idxPeaks)');

figure(2)
plot(idxPeaks, PMUTimeStamp( idxPeaks ))
xlabel('idxPeaks');
ylabel('PMUTimeStamp(idxPeaks)');

disp('Thinking: why the PMUTimeStamps are not alway the minima: 0?')
RRpeaks = TimeStamp( idxPeaks ) - PMUTimeStamp( idxPeaks );

figure(3)
plot(idxPeaks, RRpeaks)
xlabel('idxPeaks');
ylabel('RRpeaks');
% compute the interval between two consecutive corrected peaks
RRperiods = diff( RRpeaks );
figure(4)
plot((1:length(idxPeaks)-1), RRperiods)
xlabel('idx');
ylabel('RRperiods');
disp('The salient element is diff(RRpeaks(89)-RRpeaks(88))')
% plot(RRperiods, RRpeaks)

end
