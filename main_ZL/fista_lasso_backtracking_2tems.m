function [X1,X2,cost_total] = fista_lasso_backtracking_2tems(Y, template1,template2, Xinit1,Xinit2, opts)
%% Modified from the orignal FISTA code, which uses only one template.
%% Y is n*1 matrix,template is m*1 matrix, therefore X is (n-m+1)*1
    if ~isfield(opts, 'backtracking')
        opts.backtracking = false;
    end 
    opts.regul = 'l1';
    opts = initOpts(opts);
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
   
    if numel(Xinit1) == 0
        Xinit1 = zeros(size(Y,1)-size(template1,1)+1, size(Y,2));
    end
    if numel(Xinit2) == 0
        Xinit2 = zeros(size(Y,1)-size(template2,1)+1, size(Y,2));
    end
    %% cost f
    function cost = calc_f(X1,X2)
         cost = 1/2 *normF2(Y - conv(X1,template1)- conv(X2,template2));
    end 
    %% cost function 
    function cost = calc_F(X1,X2)
            cost = calc_f(X1,X2) + lambda1*norm1(X1) + lambda2*norm1(X2);  
    end 
    %% gradient
    %sigma_f(X1)=At(A*X1+B*X2-Y)
    %sigma_f(X2)=Bt(A*X1+B*X2-Y)
    function [res1,res2] = grad(X1,X2) 
        E=conv(X1,template1) + conv(X2,template2) -Y;
        Tem1E_corr=xcorr(template1,E);
        Tem2E_corr=xcorr(template2,E);
        res1=Tem1E_corr(length(E):-1:length(template1));
        res2=Tem2E_corr(length(E):-1:length(template2));
    end 
    % Checking gradient 
    if opts.check_grad
        check_grad(@calc_f, @grad, Xinit);
    end 

    opts.max_iter = 500;
    % for backtracking, we need to optimize one by one 
    X1 = zeros(size(Xinit1));
    X2 = zeros(size(Xinit2));
    for i = 1:size(X1, 2) 
       [X1(:, i),X2(:, i)] = fista_backtracking_2tems(@calc_f, @grad, Xinit1,Xinit2, opts, ...
                                        @calc_F);
    end
    cost_total=calc_F(X1,X2);
end 