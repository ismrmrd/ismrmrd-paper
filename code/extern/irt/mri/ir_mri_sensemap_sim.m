 function [smap x y z] = ir_mri_sensemap_sim(varargin)
%function [smap x y z] = ir_mri_sensemap_sim(varargin)
%|
%| Simulate 2D or 3D sensitivity maps for sensitivity-encoded MRI
%| based grivich:00:tmf doi:10.1119/1.19461
%|
%| option
%|	nx, ny, nz		image size (default: [64 64 1])
%|	dx, dy, dz		pixel/voxel dimensions (default: [3 3 3])
%|	ncoil			# of coils total (default: 4)
%|	nring			# of rings of coils (default: 1)
%|	rcoil			coil radius
%|	coil_distance		distance of coil center from isocenter for
%|				central ring of coils as a multiple of FOVx,
%|				where FOVx=nx*dx (default: 1.2)
%|	orbit			default: 360
%|
%| out
%|	smap	[nx ny nz ncoil]	simulated sensitivity maps (complex!)
%|
%| all length parameters must have same units (e.g., mm or cm)
%|
%| Copyright 2005-6-20, Jeff Fessler and Amanda Funai, University of Michigan
%| 2014-08-19 JF more testing, verifying phase is correct
%| modified for 3D by Mai 2014-09-09

if nargin == 1 % tests
	switch(varargin{1})
	case 'test1'
		ir_mri_sensemap_sim_test1 % basic test
	case 'test2'
		ir_mri_sensemap_sim_test2 % 2D test
	case 'test3'
		ir_mri_sensemap_sim_test3 % 3D test
	case 'test'
		ir_mri_sensemap_sim_test1
		ir_mri_sensemap_sim_test2
		ir_mri_sensemap_sim_test3
	otherwise
		fail('bad argument "%s"', varargin{1})
	end
return
end

arg.nx = 64;
arg.ny = [];
arg.nz = 1; % 2D
arg.dx = 3; % pixel size in mm
arg.dy = [];
arg.dz = [];
arg.ncoil = 4; % # of coils
arg.coils_per_ring = 4;
arg.nring = 1;
arg.rcoil = 100; % coil radius
arg.orbit = 360;
arg.coil_distance = 1.2; % multiplies fov/2
arg.flag_old = false; % old way based on wang:00:dop
arg.chat = nargout == 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.dy), arg.dy = arg.dx; end
if isempty(arg.dz), arg.dz = arg.dx; end
if isempty(arg.ny), arg.ny = arg.nx; end
%if isempty(arg.nz), arg.nz = arg.nx; end
if isempty(arg.rcoil), arg.rcoil = arg.dx * arg.nx / 2 * 0.50; end

coils_per_ring = arg.ncoil / nring;
assert(mod(coils_per_ringring,1) == 0, 'nring must be divisor of ncoil');
%nring = arg.ncoil/arg.coils_per_ring;
%assert(mod(nring,1) == 0, 'coils_per_ring needs to be integer factor of ncoil');

[ring_smap x y z] = ir_mri_sensemap_sim_do(arg.nx, arg.ny, arg.nz, ...
	arg.dx, arg.dy, arg.dz, ...
	coils_per_ring, coils_per_ring, arg.rcoil, ...
	arg.orbit, arg.coil_distance, ...
	arg.flag_old, arg.chat);

full_ring_fov = cat(3, ring_smap, flipdim(ring_smap(:,:,1:end-1,:),3));

block_size = arg.nz/(nring - 1);

if mod(block_size,1) ~= 0 % if any non integer
	display('warning: coil z-coordinates rounded to nearest integer');
end

smap = zeros(arg.nx, arg.ny, arg.nz, arg.ncoil);
smap(:,:,:,1:arg.coils_per_ring) = ring_smap;
for ii = 1:nring - 1
	coil_ndx = (1:arg.coils_per_ring) + ii*arg.coils_per_ring;
	full_fov_ndx = round(ii*block_size : arg.nz + ii*block_size - 1);
	%keyboard;
	smap(:,:,:, coil_ndx) = full_ring_fov(:,:,full_fov_ndx,:);
end

if ~nargout, clear, end


% ir_mri_sensemap_sim_do()
% note: can only handle coil center locations along unit sphere
function [smap x y z] = ir_mri_sensemap_sim_do(nx, ny, nz, ...
		dx, dy, dz, ncoil, ncoilpr, rcoil, ...
		orbit, coil_distance, flag_old, chat)

nring = ncoil/ncoilpr;
rlist = rcoil * ones(ncoil, 1); % coil radii

plist = zeros(ncoil,3); % position of coil center [x y z]
nlist = zeros(ncoil,3); % normal vector (inward) from coil center

% cylindrical coil configuration, like abdominal coils
alist = deg2rad(orbit)/ncoilpr * [0:(ncoilpr-1)]; % list of coil angles in radians
zcoords = linspace(0,1,nring);
%zcoords = linspace(-1,1,nring); % alternate implementation for 2 ring case
for jj = 1:nring
	for ii = 1:ncoilpr
		phi = alist(ii);
		Rad = max(nx/2 * dx, ny/2 * dy) * coil_distance;
		plist((jj-1)*ncoilpr + ii,:) = Rad * [cos(phi) sin(phi) zcoords(jj)];
		nlist((jj-1)*ncoilpr + ii,:) = -[cos(phi) sin(phi) zcoords(jj)];
		olist((jj-1)*ncoilpr + ii,:) = [-sin(phi) cos(phi) 0]; % unit vector orthogonal to nlist
		% okay to leave last coord zero
	end
end

% object coordinates
x = ([1:nx] - (nx+1)/2) * dx;
y = ([1:ny] - (ny+1)/2) * dy;
z = ([1:nz] - (nz+1)/2) * dz;
[xx yy zz] = ndgrid(x,y,z);

smap = zeros(nx, ny, nz, ncoil, 'single');
for jj = 1:nring
	for ii=1:ncoilpr
		% rotate coordinates to correspond to coil orientation
		zr =	(xx - plist((jj-1)*ncoilpr + ii,1)) .* nlist((jj-1)*ncoilpr + ii,1) + ...
			(yy - plist((jj-1)*ncoilpr + ii,2)) .* nlist((jj-1)*ncoilpr + ii,2) + ...
			(zz - plist((jj-1)*ncoilpr + ii,3)) .* nlist((jj-1)*ncoilpr + ii,3);
		xr =	xx .* nlist((jj-1)*ncoilpr + ii,2) - yy .* nlist((jj-1)*ncoilpr + ii,1);

		yr = 0;

		if 0 % see coordinates
			im plc 1 2
			im(1, x, y, xr), xlabel x, ylabel y
			im(2, x, y, zr)
			keyboard
		end

		[sx sy sz] = ir_mri_smap1(xr, yr, zr, rlist((jj-1)*ncoilpr + ii)); % in coil coordinates

		if 0 % see field components
			im plc 2 2
			im(1, x, y, sx), cbar
			im(2, x, y, sy), cbar
			im(3, x, y, sz), cbar
			im subplot 2
			tmp = sqrt(sx.^2 + sz.^2);
			quiver(x, y, (sx./tmp)', (sz./tmp)', 0), axis square
		end

		if flag_old
			smap(:,:,:,ii) = sz; % old way (wrong!) based only on smap_z

		else
			% assume z component of plist and nlist are 0 % WRONG do not assume that!
			bx = sz * nlist((jj-1)*ncoilpr + ii,1) + sx * olist((jj-1)*ncoilpr + ii,1);
			by = sz * nlist((jj-1)*ncoilpr + ii,2) + sx * olist((jj-1)*ncoilpr + ii,2);
			bz = sz * nlist((jj-1)*ncoilpr + ii,3) + sx * olist((jj-1)*ncoilpr + ii,3);
			smap(:,:,:,(jj-1)*ncoilpr + ii) = bx + 1i * by;

			if 0 % see final field components vs phase
				im subplot 4
				bb = sqrt(bx.^2 + by.^2);
				quiver(x, y, (bx./bb)', (by./bb)', 0), axis square
				im(2, x, y, angle(smap(:,:,ii))), cbar
				keyboard
			end
		end
		%if ii == 1 keyboard; end % DEBUG HERE
	end
end
smap = smap * rlist(1) / (2*pi); % trick: scale so maximum is near unity

if chat && im % show smap and array geometry in z=0 plane
	mri_sensemap_sim_show(smap, x, y, z, dx, dy, dz, nlist, plist, rlist)
end


% ir_mri_sensemap_sim_show()
function ir_mri_sensemap_sim_show(smap, x, y, z, dx, dy, dz, nlist, plist, rlist)
[nx ny ncoil] = size(smap);
nshow = min(max(ncoil,2),4);
im('plc', 3, nshow)
for ii=1:min(ncoil,nshow)
	clim = [0 max(abs(smap(:)))];
	tmp = smap(:,:,ii);
	im(ii, x, y, abs(tmp), clim, 'Magnitude'), cbar
%	xmax = max(max(abs(x)), max(plist(:,1)));
%	ymax = max(max(abs(y)), max(plist(:,2)));
	xmax = max([max(abs(x)) max(abs(y)) max(col(plist(:,[1 2])))]);
%	axis([-xmax xmax -ymax ymax]*1.05)
	axis(xmax * [-1 1 -1 1] * 1.1)
	xtick([-1 0 1] * nx/2 * dx)
	ytick([-1 0 1] * ny/2 * dy)

	hold on
	plot(0,0,'.', plist(:,1), plist(:,2), 'o')
	xdir = nlist(ii,2);
	ydir = nlist(ii,1);
	r = rlist(ii);
	plot(plist(ii,1)+r*xdir*[-1 1], plist(ii,2)+r*ydir*[1 -1], '-')
	hold off

%	ph = ir_unwrap(angle(tmp)); % trick: unwrap phase for pretty display
	ph = angle(tmp); % show raw phase (understandable with hsv colormap)
	im(ii+nshow, x, y, ph, [-pi pi], 'Phase'), cbar
	axis(xmax*[-1 1 -1 1]*1.1)
	xtick([-1 0 1] * nx/2 * dx)
	ytick([-1 0 1] * ny/2 * dy)
end

ssos = sqrt(sum(abs(smap), 3));
%if arg.chat
	minmax(ssos)
%end
ssos = ssos / ssos(end/2,end/2);
im(2*nshow+1, x, y, ssos, 'SSoS (normalized)'), cbar
xtick([-1 0 1] * nx/2 * dx)
ytick([-1 0 1] * ny/2 * dy)

if ncoil == 1
	subplot(122)
	bx = real(smap);
	by = imag(smap);
	quiver(x, y, bx', by'), title 'Field pattern in x-y plane'
	axis equal, axis tight
end

% clf, im(angle(smap(:,:,2))), cbar


% mri_smap_r(r, z)
% function for testing near 0
function out = mri_smap_r(r, z)
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);
out = 2 * z ./ r .* ((1 + r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1 - r).^2 + z.^2) .* E - K);


% ir_mri_smap1()
% based on grivich:00:tmf
% for a circular coil in "x-y plane" of radius a
% note that coil x-y plane is not same as object x-y plane!
function [smap_x smap_y smap_z] = ir_mri_smap1(x, y, z, a)
x = x ./ a; % normalized units
y = y ./ a;
z = z ./ a;
r = sqrt(x.^2 + y.^2);
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);
if ir_is_octave
	K = reshape(K, size(M)); % ellipke shape bug in octave
	E = reshape(E, size(M));
end

% the following is B_z in eqn (18) in grivich:00:tmf
% and same as eqn [10] in wang:00:dop to within constant scale factor
smap_z = 2 * ((1 + r).^2 + z.^2).^(-0.5) .* ...
	(K + (1 - r.^2 - z.^2) ./ ((1 - r).^2 + z.^2) .* E);
smap_z = smap_z / a;

if 0 && any(r(:) == 0) % test code to explore when r is near 0
	r0 = linspace(0,5e-7,101);
	z0 = 0.4;
	t0 = mri_smap_r(r0, z0);
	slope = 3*pi * z0 / ((1+z0^2)^2.5);
	clf, plot(r0, t0, '-', r0, slope * r0, '--'); grid, prompt
end

% the following is B_r in eqn (17) in grivich:00:tmf
smap_r = 2 * z ./ r .* ((1+r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1-r).^2 + z.^2) .* E - K);
bad = abs(r) < 1e-6;
smap_r(bad) = 3 * pi * z(bad) ./ ((1 + z(bad).^2).^2.5) .* r(bad);
smap_r = smap_r / a;

if any(isnan(smap_r(:))) || any(isnan(smap_z(:)))
	keyboard
end

smap_x = smap_r .* div0(x, r);
smap_y = smap_r .* div0(y, r);
%smap_z = smap_r .* div0(z, r);

%phi = atan2(y, x);
%smap_x = smap_r .* cos(phi);
%smap_y = smap_r .* sin(phi);


% ir_mri_sensemap_sim_test1
% test ir_mri_smap1 routine, cf Fig. 4 of grivich:00:tmf
function ir_mri_sensemap_sim_test1
a = 1;
x = linspace(-2,2,99);
y = linspace(-2,2,97);
zlist = [0.1 0.2 0.5 1.0];
[xx yy zz] = ndgrid(x, y, zlist);
[smap_x smap_y smap_z] = ir_mri_smap1(xx, yy, zz, a);
im('plc', 4, numel(zlist))
mri_sensemap_sim_test1_show(smap_x, x, y, 0, zlist, 'x')
mri_sensemap_sim_test1_show(smap_y, x, y, 4, zlist, 'y')
mri_sensemap_sim_test1_show(smap_z, x, y, 8, zlist, 'z')
smap_b = sqrt(smap_x.^2 + smap_y.^2);
mri_sensemap_sim_test1_show(smap_b, x, y, 12, zlist, 'b')


% ir_mri_sensemap_sim_test1_show()
function ir_mri_sensemap_sim_test1_show(map, x, y, offset, zlist, leg)
clim = [-20 20];
for iz = 1:numel(zlist)
	p = (iz-1) * 4;
	im(offset+iz, x, y, map(:,:,iz), leg, clim), cbar
	if iz == 1
		xtick([-2 2]), ytick([-2 2])
	else
		xtick off, ytick off
	end
	if zlist(iz) < 0.5
		blim = [7 12 19];
	else
		blim = [1 3 5];
	end
	hold on
	contour(x, y, abs(map(:,:,iz))', blim, 'b-')
	contour(x, y, abs(map(:,:,iz))', [0 0]+0.001, 'g-')
	hold off
end
drawnow


% ir_mri_sensemap_sim_test2
function ir_mri_sensemap_sim_test2

[smap x y] = ir_mri_sensemap_sim('chat', 1, 'nx', 32, ...
	'rcoil', [], 'ncoil', 4, 'coil_distance', 1.2);
colormap(hsv) % make it easier to see phase

if 1 % check rotational symmetry in 4-coil case
	for ic=2:4
		tmp = rot90(smap(:,:,1), (ic-1));
		equivs(abs(tmp), abs(smap(:,:,ic)))
	%	im(8+ic, x, y, abs(tmp)), cbar
		p1 = angle(tmp) + (ic-1) * pi/2; % add pi/2 to rotated
		p2 = angle(smap(:,:,ic));
		tmp = exp(1i * (p2 - p1));
		equivs(tmp, ones(size(tmp))) % trick: equivs mod 2*pi
	end
end

if 0 && im % see ellipke
	m = linspace(0,1,101);
	[k e] = ellipke(m);
	clf, plot(m, k, '-', m, e, '--'), legend('k', 'e')
	yaxis_pi('0 p/2 p 3*p/2')
end
