X1_recon=conv(X1,template1);
X2_recon=conv(X2,template2);


X1_max=get_local_maxima_above_threshold(X1,0.005,1);
X1_X2_ratio=zeros(length(X1_max),2);
figure;
plot(signal);
hold on;
scatter(X1_max,signal(X1_max),'r*')
figure;
subplot(3,1,1)
hold on;
for i=1:length(X1_max)
    if X1_max(i)+5<length(X1)&&X1_max(i)-5>0
        plot(X1(X1_max(i)-5:X1_max(i)+5),'k')
        X1_X2_ratio(i,1)=sum(X1_recon(X1_max(i):X1_max(i)+50))-51*max(X1_recon(X1_max(i):X1_max(i)+50));
    end
end
hold off;
subplot(3,1,2)
hold on;
for i=1:length(X1_max)
    if X1_max(i)+5<length(X2)&&X1_max(i)-5>0
        plot(X2(X1_max(i)-5:X1_max(i)+5),'k')
         X1_X2_ratio(i,2)=sum(X2_recon(X1_max(i):X1_max(i)+50))-51*max(X2_recon(X1_max(i):X1_max(i)+50));
    end
end
hold off;
subplot(3,1,3)
hold on;
for i=1:length(X1_max)
    if X1_max(i)+5<length(X2)&&X2(X1_max(i))<-0.005
        plot(X1(X1_max(i)-5:X1_max(i)+5)+X2(X1_max(i)-5:X1_max(i)+5),'k')
    end
end
figure;
scatter(-X1_X2_ratio(:,1),-X1_X2_ratio(:,2),'k');
chemical=-X1_X2_ratio(:,2)>-X1_X2_ratio(:,1);


figure;
subplot(4,1,1)
plot(X1);
hold on;
scatter(X1_max(chemical),X1(X1_max(chemical)),'b')
scatter(X1_max(~chemical),X1(X1_max(~chemical)),'g')
hold off;
subplot(4,1,2)
plot(X2);
hold on;
scatter(X1_max(chemical),X2(X1_max(chemical)),'b')
scatter(X1_max(~chemical),X2(X1_max(~chemical)),'g')
hold off;
subplot(4,1,3)
plot(signal,'k')
hold on;
plot(conv(X1,template1)+conv(X2,template2),'r')
scatter(X1_max(chemical),signal(X1_max(chemical)),'b')
scatter(X1_max(~chemical),signal(X1_max(~chemical)),'g')
hold off;
subplot(4,1,4)
plot(X1_recon,'b')
hold on;
plot(X2_recon,'c')
hold off;
samexaxis('abc','xmt','on','ytac','join','yld',1);