function dst = extractGTFormat(src)

src(src(:,1)>=1000,1) = 36;

act_label = src(:,1);
n_frames = [src(1,3); diff(src(:,3))];
meaning1 = src(:,4);
meaning2 = src(:,5);

dst = [act_label n_frames meaning1 meaning2];
end