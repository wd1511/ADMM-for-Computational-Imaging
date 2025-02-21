function [X, iter, min_cost] = fista_general(grad, proj, Xinit, L, opts, calc_F,Psi)   
% function [X, iter, min_cost] = fista_general(grad,proj, Xinit, L, opts, calc_F)   
% * A Fast Iterative Shrinkage-Thresholding Algorithm for 
% Linear Inverse Problems.
% * Solve the problem: X = arg min_X F(X) = f(X) + lambda*g(X) where:
%   - X: variable, can be a matrix.
%   - f(X): a smooth convex function with continuously differentiable 
%       with Lipschitz continuous gradient `L(f)` (Lipschitz constant of 
%       the gradient of `f`).
%  INPUT:
%       grad   : a function calculating gradient of f(X) given X.
%       proj   : a function calculating pL(x) -- projection
%       Xinit  : a matrix -- initial guess.
%       L      : a scalar the Lipschitz constant of the gradient of f(X).
%       opts   : a struct
%           opts.lambda  : a regularization parameter, can be either a scalar or
%                           a weighted matrix.
%           opts.max_iter: maximum iterations of the algorithm. 
%                           Default 300.
%           opts.tol     : a tolerance, the algorithm will stop if difference 
%                           between two successive X is smaller than this value. 
%                           Default 1e-8.
%           opts.verbose : showing F(X) after each iteration or not. 
%                           Default false. 
%       calc_F: optional, a function calculating value of F at X 
%               via feval(calc_F, X). 
%  OUTPUT:
%      X        : solution
%      iter     : number of run iterations
%      min_cost : the achieved cost
% Modifications:
% 06/17/2016: set default value for opts.pos = false
% -------------------------------------
% Author: Tiep Vu, thv102, 4/6/2016
% (http://www.personal.psu.edu/thv102/)
% -------------------------------------
%     opts = initOpts(opts);
    if ~isfield(opts, 'max_iter')
        opts.max_iter = 500;
    end
    if ~isfield(opts, 'regul')
        opts.regul = 'l1';
    end     
    if ~isfield(opts, 'pos')
        opts.pos = false;
    end
    
    if ~isfield(opts, 'tol')
        opts.tol = 1e-8;
    end
    
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end
    Linv = 1/L;    
    lambdaLiv = opts.lambda*Linv;
    % opts_shrinkage = opts;
    % opts_shrinkage.lambda = lambdaLiv;
    x_old = Xinit;
    y_old = Xinit;
    x_min = Xinit;
    t_old = 1;
    iter = 0;
    max_svd = 1;
    beta  = 0.01/(opts.lambda+1);
    t0 = cputime;
    times(1) = cputime - t0;
    cost_old = feval(calc_F, x_old);
    cost_min = cost_old;
    %% MAIN LOOP
    
    opts_proj = opts;
    opts_proj.lambda = lambdaLiv;
    while  iter < opts.max_iter
        iter = iter + 1;
        %fprintf('%d\n',iter);
        %x_new = feval(proj, y_old - Linv*feval(grad, y_old), opts_proj);
        x_new = Psi(y_old-Linv*feval(grad, y_old)/max_svd,0.005/max_svd);
        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        y_new = x_new + (t_old - 1)/t_new * (x_new - x_old);
        %% check stop criteria
        e = norm1(x_new - x_old)/numel(x_new);
        if e < opts.tol
            break;
        end
        %% update
        x_old = x_new;
        t_old = t_new;
        y_old = y_new;
        times(iter) = cputime-t0;
        cost1 = feval(calc_F,x_new);
        %while cost1 > cost_min
            %max_svd = 2*max_svd;
            %x_old = (1-beta)*x_old + beta*Psi(x_old-Linv*feval(grad, y_old)/max_svd,0.005/max_svd);
            %cost1 = feval(calc_F,x_old);
            %fprintf('%d\n',cost1);
        %end
        %cost_min = cost1;
        %x_min = x_new;
        %max_svd = 0.9*max_svd;
        if cost1 <= cost_min
            cost_min = cost1;
            x_min = x_new;
            max_svd = 0.9*max_svd;
        else
            max_svd = 2*max_svd;
        end
        %% show progress
        if opts.verbose
            if nargin ~= 0
                cost_new = feval(calc_F, x_new);
                if cost_new <= cost_old 
                    stt = 'YES.';
                else 
                    stt = 'NO .';
                end
                fprintf('iter = %3d, cost = %f, cost decreases? %s', ...
                    iter, cost_new, stt);
                fprintf(' ,times=%d\n',times(iter));
                cost_old = cost_new;
            else 
                if mod(iter, 5) == 0
                    fprintf('.');
                end
                if mod(iter, 10) == 0 
                   fprintf('%d', iter);
                end     
            end        
        end 
    end
    X = x_min;
    if nargout == 3 
        min_cost = feval(calc_F, X);
    end 
end 