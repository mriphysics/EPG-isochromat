figure(121)
tmp=colormap(jet(11));
%%% Add some ones
tmp = cat(1,ones([1 3]),tmp);
% interpolate
jetfade = zeros([256 3]);
for ii=1:3
    jetfade(:,ii) = interp1(linspace(0,1,size(tmp,1)),tmp(:,ii),linspace(0,1,256));
end
close(121)