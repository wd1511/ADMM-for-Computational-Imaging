function X = fista_row_sparsity1(Y, A, AT,PHI,PSI,Xinit, opts,Nvx,Nvy,Nvz)
    opts = initOpts(opts);
    lambda = opts.lambda;

    %% cost f
    function cost = calc_f(X)
        cost = 1/2 *normF2(Y - A(X));
    end 
    %% cost function 
    function cost = calc_F(X)
        cost = calc_f(X) + 0.0025*PHI(X);
    end 
    %% gradient
    function res = grad(X) 
        res = -AT(Y)+AT(A(X));
    end 
    %% Checking gradient 
    if opts.check_grad
        check_grad(@calc_f, @grad, Xinit);
    end 
    %% Lipschitz constant 
    L = 1;
    %% Use fista 
    opts.max_iter = 200;
    %calc_F(Xinit)
    %grad(Xinit)
    [X, ~, ~] = fista_general(@grad, @proj_l12, Xinit, L, opts, @calc_F,PSI);
    
    % we can also replace @proj_l12 by a function provided by SPAMS:
    % mexProximalFlat(U, opts)
    %%
    % opts.lambda = opts.lambda
    % opts.regul = 'l1l2'
    %[X, ~, ~] = fista_general(@grad, @proj_l12, Xinit, L, opts, @calc_F);    
    
end 