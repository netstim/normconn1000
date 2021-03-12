function bptc=ea_bpfilter(tc,TR,sampleLength)
    lp_HighCutoff=0.08;
    hp_LowCutoff=0.009;

%     disp('Bandpass-filtering...');
    % sampleFreq   = 1/TR;
    paddedLength = rest_nextpow2_one35(sampleLength); %2^nextpow2(sampleLength);

    mask=ones(size(tc,1),1);
    mask(:)=1;
    maskLowPass =	repmat(mask, [1, paddedLength]);
    maskHighPass=	maskLowPass;

    % GENERATE LOW PASS WINDOW	20070514, reference: fourior_filter.c in AFNI
    % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band.
    % Low pass, such as freq < 0.08 Hz
    idxCutoff	=round(lp_HighCutoff *paddedLength *TR + 1); % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band. %idxCutoff	=round(ALowPass_HighCutoff *paddedLength *TR);
    idxCutoff2	=paddedLength+2 -idxCutoff;				%Center Index =(paddedLength/2 +1)
    maskLowPass(:,idxCutoff+1:idxCutoff2-1)=0; %High eliminate

    % GENERATE HIGH PASS WINDOW
    % high pass, such as freq > 0.01 Hz
    idxCutoff	=round(hp_LowCutoff *paddedLength *TR + 1); % Revised by YAN Chao-Gan, 100420. Fixed the bug in calculating the frequency band. %idxCutoff	=round(AHighPass_LowCutoff *paddedLength *TR);
    idxCutoff2	=paddedLength+2 -idxCutoff;				%Center Index =(paddedLength/2 +1)
    maskHighPass(:,1:idxCutoff-1)=0;	%Low eliminate
    maskHighPass(:,idxCutoff2+1:paddedLength)=0;	%Low eliminate

    % 20070513	remove trend --> FFT --> filter --> inverse FFT --> retrend
    % YAN Chao-Gan, 100401. remove the mean --> FFT --> filter --> inverse FFT --> add mean back
    fftw('dwisdom');

    theMean=mean(tc,2);
    tc=tc-repmat(theMean,[1, sampleLength]);
    tc=cat(2,tc,zeros(size(tc,1),paddedLength-sampleLength));

    %FFT
    tc =fft(tc, [], 2);

    %Apply the filter Low Pass
    tc(~maskLowPass)=0;


    %Apply the filter High Pass
    tc(~maskHighPass)=0;

    %inverse FFT
    tc =ifft(tc, [], 2);
    tc =tc(:, 1:sampleLength);%remove the padded parts

    % Add the mean back after filter.
    bptc=tc+repmat(theMean,[1, sampleLength]);

%     disp('Done.');

%     if vizz
%         subplot(3,2,pcnt)
%         plot(tc(round(1:size(tc,1)/1000:size(tc,1)),:)');
%         title('Bandpass filtered time series.');
%     end


function Result = rest_nextpow2_one35(n)
%Compute the min length for FFT according to AFNI's algorithm, By Xiao-Wei Song
%------------------------------------------------------------------------------------------------------------------------------
%	Copyright(c) 2007~2010
%	State Key Laboratory of Cognitive Neuroscience and Learning in Beijing Normal University
%	Written by Xiao-Wei Song
%	http://resting-fmri.sourceforge.net
% 	<a href="Dawnwei.Song@gmail.com">Mail to Author</a>: Xiaowei Song
%	Version=1.0;
%	Release=20070903;

if length(n)>1
    n = cast(length(n),class(n));
end
if n<16
    Result =2^nextpow2(n);
    return;
end

limit =nextpow2(n);             %n=134, limit=8
tbl=[2^(limit-1):2^limit];      %tbl =128, 129, ... , 256
tbl =tbl(find(tbl>=n));          %tbl =134, 135, ... , 256
for x=1:length(tbl)
    Result =tbl(x);
    [f,p]=log2(Result);
    if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
        return;
    end
    if mod(Result,3*5)==0
        y= Result /(3*5);
        [f,p]=log2(y);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
    end
    if mod(Result,3)==0
        y= Result /3;
        [f,p]=log2(y);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
    end
    if mod(Result,5)==0
        y= Result /5;
        [f,p]=log2(y);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
    end
end
Result =NaN;    % Should not reach, except when n=1
