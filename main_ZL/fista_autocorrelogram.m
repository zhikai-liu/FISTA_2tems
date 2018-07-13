si=20;
pad=20*1e3/si;%pad length is 100ms
LM=X1_max(~chemical);
X1_X2_ele=X1_X2_ratio(~chemical,:);
clust_index=isosplit5(X1_X2_ele');

clust_num=max(clust_index);
map=colormap(jet(clust_num));
figure;
hold on;
for i=1:clust_num
    scatter(X1_X2_ele(clust_index==i,1),X1_X2_ele(clust_index==i,2),[],map(i,:));
end
hold off;
figure;
for j=1:clust_num
    cross_corr=zeros(sum(clust_index==j),2*pad+1);
    count=1;
    subplot(clust_num,1,j)
    c_xdata=LM(clust_index==j);
    spikes=zeros(1,max(c_xdata)+2*pad);
    spikes(c_xdata+pad)=1;
    dist_prox_index=[];
    for k=1:length(spikes)
        if spikes(k)==1
            bin_train=spikes(k-pad:k+pad);
            bin_train(pad+1)=0;
            dist_prox_index=[dist_prox_index,find(bin_train==1)-pad-1];
            cross_corr(count,:)=bin_train;
            count=count+1;
        end
    end
    dist_prox_index=dist_prox_index.*si/1e3;
    bin=-20:1:20;
    histogram(dist_prox_index,bin,'Normalization','pdf','FaceColor',map(j,:),'EdgeColor','none')
    xlim([-20,20])
    %ylim([0 0.3])
    AxisFormat;
end
samexaxis('ytac','box','off');
xlabel('ms')
ylabel('Probability')
AxisFormat;
function AxisFormat()
    A=gca;
    set(A.XAxis,'FontSize',20,'FontWeight','bold','LineWidth',1.2,'Color','k');
    set(A.YAxis,'FontSize',20,'FontWeight','bold','LineWidth',1.2,'Color','k');
    set(A,'box','off')
end
