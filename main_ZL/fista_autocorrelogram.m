function fista_autocorrelogram(X1_max,recon_integral,chemical)
si=20;
pad=20*1e3/si;%pad length is 100ms
LM=X1_max(~chemical);
X1_X2_ele=recon_integral(~chemical,:);
clust_index=isosplit5(X1_X2_ele');

clust_num=max(clust_index);
map=colormap(jet(clust_num));
figure;
hold on;
for i=1:clust_num
    scatter(X1_X2_ele(clust_index==i,1),X1_X2_ele(clust_index==i,2),[],map(i,:));
end
hold off;

clust=struct();
for i=1:clust_num
    clust(i).LM=LM(clust_index==i);
end
figure;
subplot_num=1;
for j=1:clust_num
    spikes_j=zeros(1,max(LM)+2*pad);
    spikes_j(clust(j).LM+pad)=1;
    for i=1:clust_num 
        spikes_i=zeros(1,max(LM)+2*pad);
        spikes_i(clust(i).LM+pad)=1;
        dist_prox_index=[];
        for k=1:length(spikes_j)
            if spikes_j(k)==1
                bin_train=spikes_i(k-pad:k+pad);
                bin_train(pad+1)=0;
                dist_prox_index=[dist_prox_index,find(bin_train==1)-pad-1];
            end
        end
        dist_prox_index=dist_prox_index.*si/1e3;
        bin=-20:1:20;
        subplot(clust_num,clust_num,subplot_num)
        if i==j
            map_color=map(j,:);
        else
            map_color=[0.3 0.3 0.3];
        end
        histogram(dist_prox_index,bin,'Normalization','pdf','FaceColor',map_color,'EdgeColor','none')
        xlim([-20,20])
        subplot_num=subplot_num+1;
    end
end
% samexaxis('ytac','box','off');
xlabel('ms')
ylabel('Probability')

figure;
subplot_num=1;
spikes_i=zeros(1,max(X1_max(chemical))+2*pad);
spikes_i(X1_max(chemical)+pad)=1;
for j=1:clust_num 
        spikes_j=zeros(1,max(LM)+2*pad);
        spikes_j(clust(j).LM+pad)=1;
        dist_prox_index=[];
        for k=1:length(spikes_j)
            if spikes_j(k)==1&&k+pad<length(spikes_i)
                bin_train=spikes_i(k-pad:k+pad);
                bin_train(pad+1)=0;
                dist_prox_index=[dist_prox_index,find(bin_train==1)-pad-1];
            end
        end
        dist_prox_index=dist_prox_index.*si/1e3;
        bin=-20:1:20;
        subplot(clust_num,1,subplot_num)
        map_color=map(j,:);
        histogram(dist_prox_index,bin,'Normalization','pdf','FaceColor',map_color,'EdgeColor','none')
        xlim([-20,20])
        subplot_num=subplot_num+1;
end

end 

function AxisFormat()
    A=gca;
    set(A.XAxis,'FontSize',20,'FontWeight','bold','LineWidth',1.2,'Color','k');
    set(A.YAxis,'FontSize',20,'FontWeight','bold','LineWidth',1.2,'Color','k');
    set(A,'box','off')
end
