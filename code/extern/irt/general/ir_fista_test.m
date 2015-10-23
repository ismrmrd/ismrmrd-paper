% ir_fista_test
% test FISTA using the quadratic programming example from O'Donoghue 2015
% min_x 1/2 x' Q x + q' x sub to a <= x <= b

if ~isvar('step')
	rng(0)
	np = 500;
	Q = randn(np+0,np);
	Q = Q' * Q;
	pr cond(Q) % 4e5
	a = -ones(np,1); % lower bound
	b = ones(np,1); % upper bound
	q = randn(np,1);
	step = 1 / eigs(Q,1)
	f.cost = @(x) 0.5 * sum(x .* (Q * x), 1) + q' * x;
	f.grad = @(x) Q * x + q;
	f.prox = @(x) max(min(x,b), a);

	x0 = zeros(np,1); % initial guess
	cost0 = f.cost(x0);
	f.niter = 2000;
end

if ~isvar('xgp'), printm 'gradient projection'
	xgp = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 0, ...
		'niter', f.niter, 'isave', 'all', 'momentum' , 0);
	cost.gp = f.cost(xgp);
end

if ~isvar('xfg'), printm 'fista'
	xfg = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 0, ...
		'niter', f.niter, 'isave', 'all');
	cost.fg = f.cost(xfg);
end

% todo: why is this not helping?
if ~isvar('xar'), printm 'adaptive restart'
	xar = ir_fista(x0, 'gradfun', f.grad, 'step', step, 'restart', 1, ...
		'niter', f.niter, 'isave', 'all', 'chat', 1);
	cost.ar = f.cost(xar);
end

ii = 0:f.niter;
plot(ii, cost.gp, '-o', ii, cost.fg, '-+', ii, cost.ar, '-s')
legend('gradient projection', 'fista', 'adaptive restart')
