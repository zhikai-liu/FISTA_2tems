load('all_traces_padded.mat');
load('EPSC_templates.mat');
Y=smooth(data_pad'-median(data_pad));
% l=6e4:12e4;
l=1:length(Y);
%l=6.48e6:6.6e6;
signal=Y(l,1);
EPSC_w1=fast_EPSC(1:441)';
EPSC_w2=slow_EPSC';
template1=EPSC_w1;
%template2=fast_EPSC(1:length(slow_EPSC))'-slow_EPSC';
alpha=EPSC_w1'*EPSC_w1/(EPSC_w1'*EPSC_w2);
template2=EPSC_w1-alpha.*EPSC_w2;
opts.backtracking=true;
opts.verbose=true;
opts.lambda1=rms(signal).*max(abs(template1)).*norminv(0.99);
opts.lambda2=rms(signal).*max(abs(template2)).*norminv(0.99);
opts.pos=false;
Xinit=[];

%[~,Xinit]=signal_deconv(Y(1:length(Y)/100,1), template,5e4,50,2000);
%X = fista_lasso_backtracking_template(Y(l,1), template, Xinit, opts);
[X1,X2,cost_matrix] = fista_lasso_backtracking_2tems(Y(l,1), template1,template2, Xinit,Xinit, opts);
%% The lasso problem is: argmin( 1/2(Y-A*x1-B*x2)^2+lambda1*|x1|+lambda2*|x2|)

% cost_matrix=zeros(10,10);
% for i=1:10
%     for j=1:10
%         opts.lambda1=i*1000;
%         opts.lambda2=j*1000;
%         [X1,X2,cost_matrix(i,j)] = fista_lasso_backtracking_2tems(Y(l,1), template1,template2, Xinit,Xinit, opts);
%     end
% end
