
% The MAD database example codes
% Citation: 
% Dong Huang, Yi Wang, Shitong Yao and F. De la Torre. Sequential Max-Margin Event Detectors, ECCV 2014
% The script is modified by yan.zhang@uni-ulm.de, to evaluate the segmentation results
function Result= funEvalDetection(gtlabE, tslab, thr, is_show)
% INPUTS:
%
% gtlab: frame-level ground truth label (obtain by loading a true label file)
% tslab: frame-level label obtained by your algorithm
% thr: threshold of overlap ratio between 
% 
% OUTPUTS:
% Result.tru_N: total number of events 
%       .dct_N: total number of detected events
%       .dct_NT:dct_NT: number of correctly detection events
%       .Prec: correctly detected events over all detected events (dct_NT/dct_N)
%       .Rec: correctly detected events over all ground truth events (dct_NT/tru_N)

%tslab(tslab==0)=36;
gtlabE(gtlabE(:,1)==36, 1)=0;
% 
class_N=length(unique(gtlabE(:,1)));
e = cumsum(gtlabE(:,2));
s = [1; e(1:end-1)+1];
seglab = [s e];

sel = (gtlabE(:,1)~=0);
gtlab = gtlabE(sel,:);
seglab = seglab(sel,:);

used_label = [];
dct_NT = 0;
tru_N = sum(sel);
for i = 1:size(seglab,1)
    tsseg = tslab(seglab(i,1):seglab(i,2));
    mj_label = mode(tsseg);
    if mj_label~=0
    if isempty(find(used_label==mj_label))
        tsseg(tsseg==mj_label)=gtlab(i,1);
        ratio = sum(tsseg==gtlab(i,1)) / gtlab(i,2);
        if ratio > thr
            dct_NT = dct_NT + 1;
        end
        used_label = [used_label,mj_label];
    end
    end
        
end


dct_N=1; %% it was 0 in the original work, but it gave errors when the result has no segment boundaries
changeframe=1; 
for i=1:length(tslab)-1
    if (abs(tslab(i+1)-tslab(i))>0) % a change point
       if  ((i-changeframe)>2) % min length
           changeframe=i; 
           if tslab(changeframe)~=0 % not null class
              dct_N=dct_N+1; 
           end
       end
    end
    
end

gt = [];
for i = 1:size(gtlabE, 1)
    gt = [gt; repmat(gtlabE(i,1), gtlabE(i,2), 1)];
end

%%% when annotation is not correct, we unify the frame length
deltan = length(gt)-length(tslab);
if deltan > 0
    gt(end-deltan+1:end) = [];
elseif deltan < 0
    tslab(end+deltan+1:end) = [];
end


[acc] = compacc(tslab,gt);
% Output------
Result.tru_N= tru_N; %total number of events 
Result.dct_N= dct_N; %total number of detected events
Result.dct_NT=dct_NT; %number of correctly detection events
Result.Prec= dct_NT/dct_N; %correctly detected events over all detected events (dct_NT/dct_N)
Result.Rec=dct_NT/tru_N;%correctly detected events over all ground truth events (dct_NT/tru_N)
Result.Acc = acc;
% gt = [];
% for i = 1:size(gtlabE, 1)
%     gt = [gt; repmat(gtlabE(i,1), gtlabE(i,2), 1)];
% end
% gt = gt(1:end-1);
% 
% [pp,rr,ff] = evaluate(tslab, gt, delta);
% randIdx = calRandIdx(tslab, gt);
% Result.precision = tru_N;
% Result.recall = rr;
% Result.f_score = ff;
% Result.RandIdx = randIdx;


% Show Bar-------
if is_show
    f = figure('Units', 'normalized', 'Position', [0,0.5,.8,0.2]);

    param.height = 1;
    param.class_N = class_N;   
    map = [[255,106,88],
            [1,64,226],
            [137,213,0],
            [130,0,190],
            [0,220,142],
            [192,0,150],
            [0,121,35],
            [0,13,119],
            [169,156,0],
            [0,138,253],
            [245,87,0],
            [3,236,253],
            [192,0,37],
            [1,187,146],
            [116,0,129],
            [181,255,178],
            [197,0,78],
            [190,255,226],
            [255,69,84],
            [0,108,85],
            [255,149,251],
            [42,87,0],
            [2,76,165],
            [175,89,0],
            [182,176,255],
            [68,64,0],
            [229,211,255],
            [30,22,0],
            [177,243,255],
            [146,0,72],
            [0,125,133],
            [255,147,156],
            [0,13,56],
            [255,236,232],
            [70,0,44],
            [0,53,46]]/255;
    cmap = colormap(map);
    cmap(1,:) = .9*[1 1 1];
    colormap(cmap);

    im_true = labelConv(gtlabE, 'slab2flab');
    im_test = tslab;

    gt = subplot(2,1,1);
    imagesc(im_true);
    % ft1 = title('');
    % set(ft1, 'FontSize', 10);
    set(gt, 'XTick', []);
    set(get(gca,'XLabel'),'String','Frame')
    set(gt, 'XTickLabel', []);
    set(gt, 'YTick', []);
    set(get(gca,'YLabel'),'String','True')
    set(gt, 'Layer', 'bottom');
    axis on
    title(['Event-based Detection Results (',num2str(thr),' overlap): Total Events=',num2str(class_N),...
           ': Precision=', num2str(Result.Prec),...
           '; Recall=', num2str(Result.Rec)])

    ts = subplot(2,1,2);
    imagesc(im_test);
    % ft2 = title('');
    % set(ft2, 'FontSize', 10);
    set(gcf, 'Color','white');
    set(ts, 'XTick', []);
    set(get(gca,'XLabel'),'String','Frame')
    set(ts, 'XTickLabel', []);
    set(ts, 'YTick', []);
    set(get(gca,'YLabel'),'String','Detected')
    set(ts, 'Layer', 'bottom');
    axis on
    end
end

% function a = findMostFrequentNumber(vec)
% ele = unique(vec);
% num = zeros(length(ele),1);
% a = 0;
% for i = 1: length(ele)
%     num = sum(vec == ele(i));
%     if num > a
%         a = num;
%     end
% end
% 
% end


function randIdx = calRandIdx(prd, gt)
a = 0; 
b = 0; 
c = 0; 
d = 0;

for ii = 1:length(gt)
    for jj = 1:length(gt)
        if prd(ii)==prd(jj)
            if gt(ii) == gt(jj)
                a = a+1;
            else
                c = c+1;
            end
        else
            if gt(ii) == gt(jj)
                b = b+1;
            else
                d = d+1;
            end
        end
    end
end
randIdx=  (a+b)/(a+b+c+d);
end





function [pp,rr,ff] = evaluate(prd, gt, delta)

%%% first, we find the index of the segment boundaries
prdf = abs(diff(prd));
gtf = abs(diff(gt));

idx_prdf = find(abs(diff(prd))>0.1);
idx_gtf = find(abs(diff(gt))>0.1);

% fprintf('--- n_segment_gt=%f  n_segment_prd = %f  \n', length(idx_gtf), length(idx_prdf));


pp = 0;
tp = 0;

%%% compute recall
for i = 1:length(idx_gtf)
    idx = idx_gtf(i);
    lb = max(1, idx-delta);
    ub = min(length(prdf), idx+delta);
    score_seg = prdf(lb:ub);
    if ~isempty(find(score_seg > 0.1))
        tp=tp+1;
    end
end
rr = tp/length(idx_gtf);

%%% compute precision
% for i = 1:length(idx_prdf)
%     idx = idx_prdf(i);
%     lb = max(1, idx-delta);
%     ub = min(length(gtf), idx+delta);
%     score_seg = gtf(lb:ub);
%     if ~isempty(find(score_seg > 0.1))
%         pp=pp+1;
%     end
% end
pp = tp/length(idx_prdf);
if isempty(idx_prdf)
    pp = 0;
    rr = 0;
end

%%% compute f measure
ff = 2*pp*rr/(pp+rr);

end

function label = labelConv(lab, mode)
%
% Convert from frame-level label to segment-level label, or vice versa.
%
% Description 
% label = labelConv(lab, mode) convert between frame-level label and
% segment-level label according to the mode.
%
% Inputs ------------------------------------------------------------------
%   o lab  : Frame-level label or segment-level label. Segment-level label
%            must be N*2, the first column is the label, the second column
%            should be segment length.
%   o mode : 2 mode. 'flab2slab' or 'slab2flab'. 
% Outputs -----------------------------------------------------------------
%   o label: label after conversion
% 
% By: Shitong Yao  // yshtng(at)gmail.com    
% Last modified: 18 July 2012
% 
if nargin < 2
    error('Two input arguments required!'); 
elseif nargin > 2
    error('Too many input arguments!');
end

if strcmpi(mode, 'flab2slab')
    % Frame-level label to segment-level label
    lab = [lab NaN];
    slab = zeros(length(lab),2);
    frame_count = 0;
    seg_count = 0;
    for i = 1:length(lab)-1
        frame_count = frame_count + 1;        
        if lab(i) ~= lab(i+1)   
            seg_count = seg_count + 1;
            slab(seg_count,:) = horzcat(lab(i), frame_count);
            frame_count = 0;   
            if i+1 == length(lab)
                break; 
            end
        end
    end
    label = slab(1:seg_count,:);  
elseif strcmpi(mode, 'slab2flab')
    % Segment-level label to frame-level label
    flab = zeros(1, sum(lab(:,2)));
    m = 0;
    for i = 1:size(lab,1)
        flab(1,m+1:m+lab(i,2)) = repmat(lab(i,1), 1, lab(i,2));
        m = m + lab(i,2);
    end
    label = flab;
else
    error('No such mode!');
end

end
