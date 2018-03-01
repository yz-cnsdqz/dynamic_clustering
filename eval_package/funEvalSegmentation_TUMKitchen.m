
% Evaluation for action segmentation based on occurrence and boundary.
% The evaluation is implemented according to the following work: 
% Dong Huang, Yi Wang, Shitong Yao and F. De la Torre. Sequential Max-Margin Event Detectors, ECCV 2014
% The script is modified by yan.zhang@uni-ulm.de

function Result= funEvalSegmentation_TUMKitchen(gt, prd, thr, is_show, varargin)
% INPUTS:
%
% gt: frame-level ground truth label (obtain by loading a true label file)
% prd: segmentation results, labels are given by clustering
% thr: threshold of overlap ratio between. thr >=1, we match boundaries not
%      regions.
% 
% OUTPUTS:
% Result.tru_N: total number of gt segments
%       .dct_N: total number of detected segments
%       .dct_NT:number of correctly partitioned actions
%       .Prec: correctly detected actions over all detected actions (dct_NT/dct_N)
%       .Rec: correctly detected actions over all ground truth actions (dct_NT/tru_N)

%%%-------------------evaluation method---------------------------------%%%
%%% (1) Evaluation is based on recall and precision.
%%% (2) In contrast to CMUMAD that actions only appear once, in this dataset
%%%     an action can repeat several times. Therefore, we cannot define the
%%%     true positives as in CMUMAD.
%%% (3) In this case, we only match the regions of overlaping for computing
%%%     recall and precision if thr < 1 and match the boundary if thr > 1. Nothing more.

diff_gt = diff(gt);
diff_prd = diff(prd);

if thr < 1

    tru_N = length(find(diff_gt))+1;
    dct_N = length(find(diff_prd))+1;
    dct_NT = 0;
    n_frames = length(gt);
    changeframe = 1;
    idx_change_left = -1;
    idx_change_right = -1;
    label_left = -1;
    label_right = -1;
    for ii = 1:tru_N
        for jj = changeframe:n_frames-1
           if abs(diff_gt(jj)) > 1e-6
               seg_gt = gt(changeframe:jj);
               seg_prd = prd(changeframe:jj);

               mj_label = mode(seg_prd);
               label_right = mj_label;
               idx_change_right = find(seg_prd==mj_label);
               seg_prd(idx_change_right)=seg_gt(1);
               idx_change_right = idx_change_right+changeframe-1;
               ratio = sum(seg_prd==seg_gt(1))/length(seg_gt);
               if ratio > thr && (idx_change_right(1)>idx_change_left(end)+1 || label_right~=label_left)
                   dct_NT = dct_NT+1;
               end
               label_left = label_right;
               idx_change_left = idx_change_right;
               changeframe = jj+1;
               break;
           end
        end
    end
elseif thr >=1
    tru_idx = find(diff_gt);
    dec_idx = find(diff_prd);
    tru_N = length(find(diff_gt));
    dct_N = length(find(diff_prd));
    dct_NT = 0;
    n_frames = length(gt);
    for ii = 1:tru_N
        idx = tru_idx(ii);
        lb = max(1,idx-thr);
        ub = min(idx+thr, n_frames-1);
        seg_prd = diff_prd(lb:ub);
        if sum(find(seg_prd)) > 0
            dct_NT = dct_NT+1;
        end
    end
end
    
    
    
% Output------
Result.tru_N= tru_N; %total number of events 
Result.dct_N= dct_N; %total number of detected events
Result.dct_NT=dct_NT; %number of correctly detection events
Result.Prec= dct_NT/dct_N; %correctly detected events over all detected events (dct_NT/dct_N)
Result.Rec=dct_NT/tru_N;%correctly detected events over all ground truth events (dct_NT/tru_N)


% Show Bar-------
if is_show
    f = figure('Units', 'normalized', 'Position', [0,0.5,.8,0.2]);  
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
    
    colormap(cmap);

%     im_true = labelConv(gtlabE, 'slab2flab');
%     im_test = tslab;
    im_true = gt';
    im_test = prd';
    gtt = subplot(2,1,1);
    imagesc(im_true);
    % ft1 = title('');
    % set(ft1, 'FontSize', 10);
    set(gtt, 'XTick', []);
    set(get(gca,'XLabel'),'String','Frame')
    set(gtt, 'XTickLabel', []);
    set(gtt, 'YTick', []);
    set(get(gca,'YLabel'),'String','True')
    set(gtt, 'Layer', 'bottom');
    axis on
    title([' Results (',num2str(thr),' overlap): ',...
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
    if nargin == 8
        method = varargin{1};
        bodypart = varargin{2};
        feature = varargin{3};
        video_idx = varargin{4};
        figname = sprintf('%s_%s_%s_%d',method,bodypart, feature, video_idx);
        savefig(f, figname);
    end
    
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

% function label = labelConv(lab, mode)
% %
% % Convert from frame-level label to segment-level label, or vice versa.
% %
% % Description 
% % label = labelConv(lab, mode) convert between frame-level label and
% % segment-level label according to the mode.
% %
% % Inputs ------------------------------------------------------------------
% %   o lab  : Frame-level label or segment-level label. Segment-level label
% %            must be N*2, the first column is the label, the second column
% %            should be segment length.
% %   o mode : 2 mode. 'flab2slab' or 'slab2flab'. 
% % Outputs -----------------------------------------------------------------
% %   o label: label after conversion
% % 
% % By: Shitong Yao  // yshtng(at)gmail.com    
% % Last modified: 18 July 2012
% % 
% if nargin < 2
%     error('Two input arguments required!'); 
% elseif nargin > 2
%     error('Too many input arguments!');
% end
% 
% if strcmpi(mode, 'flab2slab')
%     % Frame-level label to segment-level label
%     lab = [lab NaN];
%     slab = zeros(length(lab),2);
%     frame_count = 0;
%     seg_count = 0;
%     for i = 1:length(lab)-1
%         frame_count = frame_count + 1;        
%         if lab(i) ~= lab(i+1)   
%             seg_count = seg_count + 1;
%             slab(seg_count,:) = horzcat(lab(i), frame_count);
%             frame_count = 0;   
%             if i+1 == length(lab)
%                 break; 
%             end
%         end
%     end
%     label = slab(1:seg_count,:);  
% elseif strcmpi(mode, 'slab2flab')
%     % Segment-level label to frame-level label
%     flab = zeros(1, sum(lab(:,2)));
%     m = 0;
%     for i = 1:size(lab,1)
%         flab(1,m+1:m+lab(i,2)) = repmat(lab(i,1), 1, lab(i,2));
%         m = m + lab(i,2);
%     end
%     label = flab;
% else
%     error('No such mode!');
% end
% 
% end
