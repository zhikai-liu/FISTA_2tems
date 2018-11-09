function fista_autocorrelogram(X1_max,clust_index)
si=20;%this is the time per data point
show_ms=100;%prox time window for x axis, 100ms if 100
pad=show_ms*1e3/si;%pad length is data point length
time_res=1;% time resolution for the histogram count
bin=-show_ms:time_res:show_ms;
clust_num=max(clust_index);% Total cluster numbers
map=colormap(jet(clust_num));% Colormap for different clusters
clust=struct();
%% Separate different cluster events out
for i=1:clust_num
    clust(i).LM=X1_max(clust_index==i);
end
%% Correlogram
% The idea to auto-correlate or cross-correlate is to look at a single event, and count probability distribution of events nearby from
% either within the cluster or other clusters
figure;
subplot_num=1;
for j=1:clust_num
    spikes_j=zeros(1,max(X1_max)+2*pad);
    spikes_j(clust(j).LM+pad)=1;
    for i=1:clust_num 
        spikes_i=zeros(1,max(X1_max)+2*pad);
        spikes_i(clust(i).LM+pad)=1;
        dist_prox_index=[];
        for k=1:length(spikes_j)
            %% Look at every spike first
            if spikes_j(k)==1
                % For each spike, use it as a reference, and examine a small spike train within (pad)
                % ms for that spike
                bin_train=spikes_i(k-pad:k+pad);
                % Make the reference spike to 0 so it will not be counted
                bin_train(pad+1)=0;
                % Calculate the distance of other spikes to the reference
                % spike
                dist_prox_index=[dist_prox_index,find(bin_train==1)-pad-1];
            end
        end
        dist_prox_index=dist_prox_index.*si/1e3;
        subplot(clust_num,clust_num,subplot_num)
        if i==j
            map_color=map(j,:);% if auto-correlogram, color-coded by cluster colors
        else
            map_color=[0.3 0.3 0.3];% if cross-correlogram, colored by gray
        end
        histogram(dist_prox_index,bin,'Normalization','pdf','FaceColor',map_color,'EdgeColor','none')
        xlim([bin(1),bin(end)])
        subplot_num=subplot_num+1;
    end
end
% samexaxis('ytac','box','off');
xlabel('ms')
ylabel('Probability')

end 
