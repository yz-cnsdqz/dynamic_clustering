
% used in our evaluation. Uncomment for reproduce the result in the paper
function dst = prdInterpolation(src, S, n_frames)
dst = zeros(n_frames,size(src,2));
n_snippets = size(src,1);
dst(1:n_snippets,:) = src;
delta = n_frames-n_snippets;
dst(n_snippets+1:n_frames,:) = dst(n_snippets:-1:n_snippets-delta+1, :);
dst(isnan(dst)) = 0;
end


% %% used in demo for better visualization
% function dst = prdInterpolation(src, S, n_frames)
% dst = zeros(n_frames,size(src,2));
% n_snippets = size(src,1);
% dst(end-n_snippets+1:end,:) = src;
% delta = n_frames-n_snippets;
% dst(1:delta,:) = dst(2*delta:-1:delta+1, :);
% dst(isnan(dst)) = 0;
% 
% end