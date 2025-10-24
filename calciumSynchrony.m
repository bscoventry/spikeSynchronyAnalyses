%--------------------------------------------------------------------------
% Author: Brandon S Coventry, PhD    Wisconsin Institute for Translational
% Neuroengineering
% Date: 10/24/25
% Purpose: Implements Vector strength calculations to measure phase
% synchrony, ie the timing of an event related to a single period of the
% signal. Also implements circular statistical measures to assess
% significantly modulated responses
% Revision History: Code modified from some of my earlier scripts.
%--------------------------------------------------------------------------
function [vs,phasemean,rayleigh] = vectorstrength(ts,period)
    %{
    Calculates the vector strength of a time series. 
    Inputs:
        ts - array or matrix of time points corresponding to spikes. ts
        can be a single time series array, where columns correspond to
        spike times. Alternatively, ts can also be a N-size cell, where
        each row cell value is an array containing an individual cells
        spike times. For example
        spikes = cell{Nx1}
        spikes{1} = [0.2, 0.34, 0.5, 0.8]
        spikes{2} = [0.22, 0.34, 0.83]
        period - Dominant period of the signal. 
    Outputs:
        vs - If array, single value output of the vector strength. If cell,
        an output array containing vector strengths for N cells contained
        in put cell.
        phasemean - Output mean phase of the time series. If cell, outputs
        array of length N that contains the mean phase for each cell.
        raleigh - Output rayleigh statistic. Significance is a rayleigh
        statistic above 13.8 (p<0.05).
    %}
    if iscell(ts)
        nr = length(ts);
        vs = empty(nr,1);
        phasemean = empty(nr,1);
        rayleigh = empty(nr,1);
        for ck = 1:nr
            curts = ts{ck};
            if isempty(curts);
                vs(ck) = NaN;
                phasemean(ck) = NaN;
                rayleigh(ck) = NaN;
            else
                %Calculate phase response
                phases = mod(curts, period)./(period*2*pi);
                %Convert to complex exponential. Can do Euler's ID too
                comEx = sum(exp(-j*phases));
                vs(ck) = abs(comEx)/length(curts);
                curpm = angle(comEx);
                if curpm < 0       %Catch this so it's from 0-360, not 180-180
                    curpm = curpm + 2 * pi;
                end
                phasemean(ck) = curpm;
                %Now calculate Rayleigh statistic
                rayleigh(ck) = length(curts)*vs(ck).^2;
            end
        end
    else
        %This is the case where you feed in just a single cell array
        if isempty(ts);
            vs = NaN;
            phasemean = NaN;
            rayleigh = NaN;
        else
            nc = length(ts);
            %Calculate phase response
            phases = mod(ts, period)./(period*2*pi);
            %Convert to complex exponential. Can do Euler's ID too
            comEx = sum(exp(-j*phases));
            vs = abs(comEx)/length(nc);
            curpm = angle(comEx);
            if curpm < 0       %Catch this so it's from 0-360, not 180-180
                curpm = curpm + 2 * pi;
            end
            phasemean = curpm;
            %Now calculate Rayleigh statistic
            rayleigh = nc*vs.^2;
        end
    end
end
%{
Okay, this is great, but what if we don't know the dominant period? Very
space bin cross correlation can get us there! This function will create
distributions for all time series cells. 
%}
function [meanPeriod,totalPeriods] = predictPeriod(ts);
%{
This is a function to estimate periods within and between spike trains.
Uses derivative measures to observe a mean shift in response as well as
changes.
Inputs
 ts - input time series of spike times. Input ts as a cell from the vs
 function above if you want to analyze multiple cells.
Outputs
 meanPeriod - Outputs mean period between all events. If input cell, output
 array of mean periods for each cell
 totalPeriod - Ouputs an array of the total period calculations. Useful for
 if the cell synchrony drifts. If input is a cell array, output is a cell
 array with each index corresponding to an array for each cell.
%}
    if iscell(ts)
        nr = length(ts);
        meanPeriod = empty(nr,1);
        totalPeriods = cell(nr,1);
        for ck = 1:nr
            curResp = ts{ck};
            if length(curResp)>1
                diffPeriod = diff(curResp);
                meanPeriod(ck) = mean(diffPeriod);
                totalPeriods{ck} = diffPeriod;
            else
                meanPeriod(ck) = NaN;
                totalPeriods{ck} = NaN;
            end
        end
    else
        if length(ts)>1
            diffPeriod = diff(ts);
            meanPeriod = mean(diffPeriod);
            totalPeriods = diffPeriod;
        else
            meanPeriod = NaN;
            totalPeriods = NaN;
        end
    end
end
%Use of these functions.
%This is a short tutorial on how to perform this analysis.
%{
Start by setting up your calcium events as calcium event times. If you want
to analyze all cells, organize this into a matlab cell (confusing on using
cell twice here, I agree). I'll call biocells neurons here. Organize your
neurons into the following matlab cell structure
neuronResponses = cell(Nneuronsx1)
neuronResponses{1} = spike times neuron 1
neuronResponses{2} = spike times neuron 2 etc.

Estimate the period of synchronous activity.
[meanPeriod,totalPeriods] = predictPeriod(neuronResponses);
This will output the mean periods of each neuron, and then the total
difference time series for all neurons. Evaluate if there is a dominate
period and if that stays consistent. Now, calculate vector strength

meanmeanPeriod = mean(meanPeriod) to "quess" at the average synchony.

[vs, phasemean, rayleigh] = vectorstrength(NeuronResponses,meanmeanPeriod)
You'll get the output vector strengths and each neurons mean phase
difference from the overall synchony period. Rayleigh is a statistic, from
a circular distribution. Effectively a significance test for when things
are mapped around a circle, such as phase response. Tests the null
hypothesis that the time series is not synchronized to any particular
phase. Rayleigh statistic above 13.8 is significant (p<0.05). This can be
used to determine how cells loose synchony and how they may regain it
again.
%}