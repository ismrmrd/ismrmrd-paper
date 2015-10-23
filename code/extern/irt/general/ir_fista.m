 function [xs, info] = ir_fista(x, varargin)
%function [xs, info] = ir_fista(x, varargin)
%|
%| Generic FISTA method
%|
%| in
%|	x	[np 1]		initial estimate
%|
%| required
%|	'gradfun' @()		function returning gradient: gradfun(x)
%|	'step'			user-selected step size (1/Lipschitz)
%|
%| option
%|	'proxfun' @()		function performing proximal step
%|					(default: @(x)=x)
%|					e.g. @(x) = max(x,0)
%|					for nonnegativity constraint
%|	'niter'			# total iterations (default 1)
%|	'chat'			verbosity (default 0)
%|	'restart'	0|1	1 to use restart (default: 1)
%|	'momentum'	0|1	set to 0 to eliminate momentum (default: 1)
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 3]	step, time, cost
%|
%| Copyright 2015-07-01, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('ir_fista_test'), return
end
if nargin < 3, help(mfilename), error(mfilename), end

arg.gradfun = [];
arg.proxfun = @(x) x;
arg.niter = 1;
arg.step = [];
arg.chat = 0;
arg.isave = [];
arg.restart = true;
arg.alpha_restart = -cos(80*pi/180); % 80 degrees
arg.momentum = 1;

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

if isempty(arg.gradfun), fail('gradfun required'), end
if isempty(arg.step), fail('step required'), end

if ~isa(arg.gradfun, 'function_handle'), error 'gradfun not function handle?', end
if ~isa(arg.proxfun, 'function_handle'), error 'proxfun not function handle?', end

xs = zeros(numel(x), length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

info = zeros(arg.niter, 3);

% fast iterative shrinkage/thresholding algorithm (FISTA)
v = x;
xold = x; told = 1;
cpu etic
for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	grad = arg.gradfun(v);

	x = v - arg.step * grad;
	x = arg.proxfun(x);

	% adaptive restart: todo - why not working?
	tmp = ((v - x)' * (x - xold)) / norm(v - x) / norm(x - xold);
%	pr tmp
	if arg.restart && (tmp > arg.alpha_restart)
		t = 1;
		v = x;
		if arg.chat
			printm('restart at %d', iter)
		end
	else
		t = (1 + sqrt(1 + 4 * told^2)) / 2;
		frac = arg.momentum * (told-1) / t;
		v = x + frac * (x - xold);
	end

	xold = x; told = t;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end

	info(iter,:) = cpu('etoc');
end
