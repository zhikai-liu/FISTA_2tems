load('all_traces_padded.mat');
load('EPSC_templates.mat');
Y=smooth(data_pad'-mean(data_pad));
template1=fast_EPSC';
%template2=fast_EPSC(1:length(slow_EPSC))'-slow_EPSC';
template2=slow_EPSC';
opts.backtracking=true;
opts.verbose=false;
% opts.lambda1=4000;
% opts.lambda2=9000;
opts.pos=true;
Xinit=[];
%l=6.1e6:6.7e6;
l=6e4:12e4;
%[~,Xinit]=signal_deconv(Y(1:length(Y)/100,1), template,5e4,50,2000);
%X = fista_lasso_backtracking_template(Y(l,1), template, Xinit, opts);
% [X1,X2,cost_matrix] = fista_lasso_backtracking_2tems(Y(l,1), template1,template2, Xinit,Xinit, opts);
%% The lasso problem is: argmin( 1/2(Y-A*x1-B*x2)^2+lambda1*|x1|+lambda2*|x2|)

cost_matrix=zeros(10,10);
for i=1:10
    for j=1:10
        opts.lambda1=i*2000;
        opts.lambda2=j*2000;
        [X1,X2,cost_matrix(i,j)] = fista_lasso_backtracking_2tems(Y(l,1), template1,template2, Xinit,Xinit, opts);
    end
end
% figure;
% subplot(4,1,1)
% plot(X1);
% subplot(4,1,2)
% plot(X2);
% subplot(4,1,3)
% plot(Y(l,1),'k')
% hold on;
% plot(conv(X1,template1)+conv(X2,template2),'r')
% hold off;
% subplot(4,1,4)
% plot(conv(X1,template1),'b')
% hold on;
% plot(conv(X2,template2),'c')
% hold off;
% samexaxis('abc','xmt','on','ytac','join','yld',1);
heatmap(cost_matrix)