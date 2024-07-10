function [somm_disc] = get_sommerfeld_disc(zks, xylim, tol, opts)
%  This subroutine returns the sommerfeld discretization for 
%  evaluating corrections to half-space greens' function
%  using alpert corrections
%
%  Input arguments:
%    zks: complex(2)
%      wave numbers in the upper and lower half planes
%      respectively
%    xylim: double (2,2)
%      [xmin, xmax; ymin, ymax]
%      x extent computed using xmax - xmin
%      ymin determines closest source, and ymax
%      determines the furthest out sommerfeld integral
%      is evaluated
%    tol: double 
%      tolerance
%    opts: options struct (optional)
%      opts.norder (16) - order of alpert rule correction
%      opts.nover (20) - oversampling paarameter 
%

   if nargin < 4
     opts = [];
   end

   norder = 16;
   if isfield(opts, 'norder')
      norder = opts.norder;
   end

   nover = 20;
   if isfield(opts, 'nover')
     nover = opts.nover;
   end

   if nargin < 3
     tol = 1e-8;
   end

   somm_disc = [];
   xmin = xylim(1,1);
   xmax = xylim(1,2);

   ymin = xylim(2,1);
   ymax = xylim(2,2);

   xc = 0.5*(xmin + xmax);
   xext = 0.5*(xmax - xmin);

   nmax = 100000;
   xnodes = zeros(nmax,1);
   weights1 = zeros(nmax,1);
   xfac1 = zeros(nmax,1);
   yfac1 = zeros(nmax,1);

   weights2 = zeros(nmax,1);
   xfac2 = zeros(nmax,1);
   yfac2 = zeros(nmax,1);

   zk1 = zks(1);
   zk2 = zks(2);

   tmax = 1.3d0*real(zk2) + abs(log(tol))/ymin + 100;
   ntot = 0;

   mex_id_ = 'getsommerquad_2d2mi(i dcomplex, i dcomplex, i dcomplex, i double, i int, i double, i double, i double, i int, io int[x], i int, io double[x], io dcomplex[x], io dcomplex[x], io dcomplex[x])';
[ntot, xnodes, weights1, xfac1, yfac1] = somm2drouts(mex_id_, zk1, zk2, zk1, tmax, norder, tol, xext, ymax, nover, ntot, nmax, xnodes, weights1, xfac1, yfac1, 1, nmax, nmax, nmax, nmax);
   mex_id_ = 'getsommerquad_2d2mi(i dcomplex, i dcomplex, i dcomplex, i double, i int, i double, i double, i double, i int, io int[x], i int, io double[x], io dcomplex[x], io dcomplex[x], io dcomplex[x])';
[ntot, xnodes, weights2, xfac2, yfac2] = somm2drouts(mex_id_, zk1, zk2, zk2, tmax, norder, tol, xext, ymax, nover, ntot, nmax, xnodes, weights2, xfac2, yfac2, 1, nmax, nmax, nmax, nmax);

   somm_disc.n = ntot;
   somm_disc.x = xnodes(1:ntot) + xc;
   somm_disc.w1 = weights1(1:ntot);
   somm_disc.xfac1 = xfac1(1:ntot);
   somm_disc.yfac1 = yfac1(1:ntot);

   somm_disc.w2 = weights2(1:ntot);
   somm_disc.xfac2 = xfac2(1:ntot);
   somm_disc.yfac2 = yfac2(1:ntot);


end
