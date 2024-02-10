c
c
      subroutine getsommerquad_2d2(zk1,zk2,zk,tmax,norder,tol,rx,ry,
     1       nover,ntot,nmax,xnodes,weights,xfac,yfac)
c
c     Discretize Sommerfeld integral to respect two different Helmholtz
c     parameters: zk1, zk2 with real(zk2) > real(zk1), assuming
c     zk is wavenumber of interest (could be zk1 or zk2). 
c     The scheme respects the sqrt singularities at +/- Re(zk1,zk2)
c     and uses 5 intervals with Alpert gauss-trapezoidal weights
c     on 
c     [-tmax,-Re(zk2)], [-Re(zk2),-Re(zk1)]
c     [-Re(zk1),Re(zk1)], [Re(zk1),Re(zk2)],[Re(zk2),tmax].
c
c     ASSUMES REAL(ZK2) > REAL(ZK1),
c
c     Since the scheme respects the sqrt singularities at 
c     +/- Re(zk1,zk2), it will provide poor precision for other 
c     values of zk.
c
c        
c g_k = 1/(4 pi) \int exp(-\sqrt{lambda^2-zk^2}y exp(i lambda x) dlambda 
c                     -----------------------------------------
c                              \sqrt{lambda^2-zk^2}
c        
c
c     INPUT:
c
c     zk1     Helmholtz parameter for upper half space
c     zk2     Helmholtz parameter for lower half space
c             with Re(zk2) > Re(zk1).
c     zk      Helmholtz parameter for specific integral being computed
c     tmax    limits on integration in lambda (in evanescent region)
c     norder  order of Alpert rule corrections
c     tol     desired precision
c     rx      max excursion in x
c     ry      max excursion in y
c     nover   oversampling parameter (currently set by user)
c     nmax    dimensions of input arrays xnodes,weights,xfac,yfac
c
c     OUTPUT:
c
c     ntot    total number of quadrature nodes in Sommerfeld integral
c     xnodes  node locations
c     weights quad. weights, includes 1/(4 pi) factor and weight for
c             dlambda but NOT 1/\sqrt{lambda^2-k^2}
c             factor in denominator
c             THUS - ACTUALLY REAL BUT DIMENSIONED AS COMPLEX
c             IN CASE OF FUTURE ADJUSTMENT.
c     xfac    in this case, simply discrete \lambda values,
c             (a copy of xnodes). Keeping the array in case of 
c             some future change of variables
c     yfac    the expression \sqrt{k^2 - lambda^2} at the lambda nodes
c             arranged so that branch cuts are the standard ones
c             (see notes).
c             
c             More precisely, the numerator of the integrand above is
c             cdexp(eye*yfac(i)*y)*cdexp(eye*xfac(i)*x)
c
c     The routines manychgdp_somm2d and somm2d_eval are compatible 
c     with this formalism.
c        
c

      implicit none
      integer ntot,nquad,norder,i,ii,nmax,ntot_loc,iflag,ntemp,nover
      integer nintervals,nminquad,nquad11,nquad12,nquad_decay
      real *8 xnodes(nmax)
      complex *16 weights(nmax)
      real *8 rkm1,rkp1,tmax,pi,rk,aa,bb
      real *8 rkm2,rkp2,tol,rx,ry
      complex *16 xfac(nmax)
      complex *16 yfac(nmax)
      complex *16 zk1,zk2,eye,zarg1,zarg2,zk
      real *8, allocatable :: xnodes_loc(:)
      real *8, allocatable :: weights_loc(:)
c
      allocate(xnodes_loc(2*nmax))
      allocate(weights_loc(2*nmax))
      nminquad = 24
      eye = dcmplx(0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      rkm1 = -real(zk1)
      rkp1 = real(zk1)
      rkm2 = -real(zk2)
      rkp2 = real(zk2)
c
c     ray to - \infty
c
      aa = -tmax
      bb =  rkm2
      nquad_decay = max(24,nint((bb-aa)*nover*(rx+ry)/(2*pi)))
      call prinf(' nquad_decay is *',nquad_decay,1)
      call alpert_sqrtb(aa,bb,nquad_decay,norder,
     1     xnodes_loc,weights_loc)
      do i = 1,nquad_decay
         xnodes(i) = xnodes_loc(i)
         rk = xnodes(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(i) = zarg1*zarg2
         xfac(i) = rk
ccc         weights(i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = nquad_decay
c
c     the interval [-real(zk2),-real(zk1)]
c
      aa = rkm2
      bb = rkm1
c
      nquad12 = max(24,nint((bb-aa)*nover*(rx+ry)/(2*pi)))
      call prinf(' nquad12 is *',nquad12,1)
      call alpert_sqrtab(aa,bb,nquad12,norder,
     1     xnodes_loc,weights_loc)

      do i = 1,nquad12
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
ccc         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntemp + nquad12
c
c     the interval [-real(zk1),real(zk1)]
c
      aa = rkm1
      bb = rkp1
      nquad11 = max(24,nint((bb-aa)*nover*(rx+ry)/(2*pi)))
      call prinf(' nquad11 is *',nquad11,1)
      call alpert_sqrtab(aa,bb,nquad11,norder,
     1     xnodes_loc,weights_loc)
      do i = 1,nquad11
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
ccc         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntemp + nquad11
c
c     the interval [real(zk1),real(zk2)]
c
      aa = rkp1
      bb = rkp2
      call alpert_sqrtab(aa,bb,nquad12,norder,
     1     xnodes_loc,weights_loc)
      do i = 1,nquad12
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
ccc         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntemp + nquad12
c
c     the ray to + \infty
c 
      aa = rkp2
      bb = tmax
      call alpert_sqrta(aa,bb,nquad_decay,norder,
     1     xnodes_loc,weights_loc)
      do i = 1,nquad_decay
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
ccc         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntot = ntemp + nquad_decay
c
      return
      end
c
c
c
        subroutine onechgdp_somm2d(zk,iup,source,chg,dipstr,dipvec,
     1                center,xfac,yfac,weights,sommer_exp,nquad)
c
c        
c g_k = 1/(4 pi) \int exp(-iup*\sqrt{lambda^2-k^2}y exp(i lambda x) dlambda 
c                     -----------------------------------------
c                              \sqrt{lambda^2-k^2}
c 
c       This routine converts a charge and dipole into a Sommerfeld
c       expansion valid for targets with y-coordinate greater
c       than the y-coordinate of the source (if iup = 1) and
c       for targets with y-coordinate 
c       less than the y-coordinate of the source (if iup = -1).
c
c       INPUT:
c
c       zk       Helmholtz parameter
C       iup      1 = up, -1 = down  (defines decay direction)
c       source   source location
c       chg      charge strength
c       dipstr   dipole strength
c       dipvec   dipole vector
c       center   Sommerfeld expansion origin
c       xfac, yfac  discretization arrays from quadrature generator
c       weights  weights from quadrature generator
c       nquad    number of quadrature nodes
c
c       OUTPUT:
c
c       sommer_exp  moment function for discretized Sommerfeld integral 
c       THIS IS A FORMULA AS A FUNCTION OF LAMBDA AND DOES NOT
c       INVOLVE QUAD WEIGHTS
c        
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,iup
      real *8 center(2),source(2),dipvec(2) 
      complex *16 chg,dipstr
      complex *16 xfac(nquad),yfac(nquad),sommer_exp(nquad)
      complex *16 weights(nquad)
      complex *16 zk,zinc,eye,zinc2,zp,zm
      complex *16 zarg1,zarg2,zarg,zdot
c
      pi = 4.0d0*datan(1.0d0)
      eye = (0,1)
      do i = 1,nquad
         zdot = dipvec(1)*eye*xfac(i) + iup*eye*dipvec(2)*yfac(i)
         zinc = cdexp(-eye*xfac(i)*(source(1)-center(1)))*
     1        cdexp(-iup*eye*yfac(i)*(source(2)-center(2)))
ccc         zinc = zinc*weights(i)/(4.0d0*pi)
         sommer_exp(i) = chg*zinc
         sommer_exp(i) = sommer_exp(i) - dipstr*zinc*zdot
      enddo
      return
      end
c
c
c
      subroutine somm2d_eval(zk,iup,ztrg,xnodes,center,xfac,yfac,
     1                weights,sighat,ntot,zsum,grad)
c
c        
c zsum = 1/(4 pi) \int exp(-SGN*\sqrt{lambda^2-k^2}y 
c                             exp(i lambda x) * sighat dlambda 
c                     -----------------------------------------
c                              \sqrt{lambda^2-k^2}
c 
c     This routine evaluates the Sommerfeld integral
c     at ztrg = (x,y) with SGN= iup (input parameter)
c
c     INPUT:
c
c     zk       Helmholtz parameter
C     iup      1 = up, -1 = down  (defines decay direction)
c     ztrg     target location
c     xnodes   node locations
c     center   Sommerfeld expansion origin
c     xfac, yfac  discretization arrays from quadrature generator
c     sighat   moment function for discretized Sommerfeld integral 
c              which can be thought of as Fourier transform of SLP
c     ntot     number of quadrature nodes 
c
c     OUTPUT:
c
c     zsum     value of Sommerfeld integral
c        
      implicit none
      integer ntot,i,iup,ii
      real *8 xnodes(ntot)
      real *8 center(2)
      real *8 ztrg(2),x,y
      complex *16 zk,zsum,eye
      complex *16 grad(2)
      complex *16 weights(ntot)
      complex *16 sighat(ntot)
      complex *16 xfac(ntot)
      complex *16 yfac(ntot)
c
      write(6,*) 'sighat ',sighat(ntot/2)
      write(6,*) 'weights ',weights(ntot/2)
      eye = dcmplx(0.0d0,1.0d0)
      x = ztrg(1) -center(1)
      y = ztrg(2) -center(2)
      zsum = 0.0d0
      grad(1) = 0.0d0
      grad(2) = 0.0d0
      do i = 1,ntot
         zsum = zsum + sighat(i)*weights(i)*
     1      cdexp(eye*iup*yfac(i)*y)*
     1       cdexp(eye*xfac(i)*x)/(-eye*yfac(i))
         grad(1) = grad(1) + sighat(i)*weights(i)*
     1      cdexp(eye*iup*yfac(i)*y)*eye*xfac(i)*
     1       cdexp(eye*xfac(i)*x)/(-eye*yfac(i))
         grad(2) = grad(2) + sighat(i)*weights(i)*
     1      iup*eye*yfac(i)*cdexp(eye*iup*yfac(i)*y)*
     1       cdexp(eye*xfac(i)*x)/(-eye*yfac(i))
      enddo
c
      return
      end
c
c
c
c
c
      subroutine somm2d_eval_dp(zk,iup,ztrg,xnodes,center,xfac,yfac,
     1                weights,zmuhat,ntot,zsum,grad)
c
c        
c     zsum = 1/(4 pi) *SGN* \int exp(-SGN*\sqrt{lambda^2-k^2}y 
c                               exp(i lambda x) * zmuhat dlambda 
c 
c     This routine evaluates the Sommerfeld representation of 
c     double layer at ztrg = (x,y) with SGN= iup (input parameter)
c
c     INPUT:
c
c     zk       Helmholtz parameter
C     iup      1 = up, -1 = down  (defines decay direction)
c     ztrg     target location
c     xnodes   node locations
c     center   Sommerfeld expansion origin
c     xfac, yfac  discretization arrays from quadrature generator
c     zmuhat    moment function for discretized Sommerfeld integral 
c              which can be thought of as Fourier transform of DLP
c     ntot     number of quadrature nodes 
c
c     OUTPUT:
c
c     zsum     value of Sommerfeld integral
c     grad     of zsum
c        
      implicit none
      integer ntot,i,iup,ii
      real *8 xnodes(ntot)
      real *8 center(2)
      complex *16 weights(ntot)
      real *8 ztrg(2),x,y
      complex *16 zk,zsum,eye
      complex *16 grad(2)
      complex *16 zmuhat(ntot)
      complex *16 xfac(ntot)
      complex *16 yfac(ntot)
c
      eye = dcmplx(0.0d0,1.0d0)
      x = ztrg(1) -center(1)
      y = ztrg(2) -center(2)
      zsum = 0.0d0
      grad(1) = 0.0d0
      grad(2) = 0.0d0
      do i = 1,ntot
         zsum = zsum + iup*zmuhat(i)*weights(i)*
     1       cdexp(eye*iup*yfac(i)*y)*
     1       cdexp(eye*xfac(i)*x)
         grad(1) = grad(1) + iup*zmuhat(i)*weights(i)*
     1       cdexp(eye*iup*yfac(i)*y)*eye*xfac(i)*
     1       cdexp(eye*xfac(i)*x)
         grad(2) = grad(2) + iup*zmuhat(i)*weights(i)*
     1       cdexp(eye*iup*yfac(i)*y)*eye*iup*yfac(i)* 
     1       cdexp(eye*xfac(i)*x)
      enddo
c
      return
      end
c

c
        subroutine manychgdp_somm2d(zk,iup,ns,source,chg,dipstr,
     1             dipvec,center,xfac,yfac,weights,sommer_exp,nquad)
c
c        
c       Based on Sommerfeld representation:
c
c g_k = 1/(4 pi) \int exp(-iup*\sqrt{lambda^2-k^2}y exp(i lambda x) dlambda 
c                     -----------------------------------------
c                              \sqrt{lambda^2-k^2}
c 
c       This routine converts a collection of charges and dipoles into
c       a Sommerfeld expansion valid for targets with y-coordinate 
c       greater than the y-coordinate of the source (if iup = 1) and
c       for targets with y-coordinate 
c       less than the y-coordinate of the source (if iup = -1).
c
c       This is a slow, direct algorithm. If there are many sources,
c       this can be accelerated with anterpolation and nufft.
c
c       INPUT:
c
c       zk       Helmholtz parameter
C       iup      1 = up, -1 = down  (defines decay direction)
c       ns       number of sources
c       source   source locations
c       chg      charge strengths
c       dipstr   dipole strengths
c       dipvec   dipole vectors
c       center   Sommerfeld expansion origin
c       xfac, yfac  discretization arrays from quadrature generator
c       weights  weights from quadrature generator
c       nquad    number of quadrature nodes
c
c       OUTPUT:
c
c       sommer_exp  moment function for discretized Sommerfeld integral 
c
c       That is, the field due to the charges/dipoles is
c
c       u = 1/(4 pi) \int exp(-iup*\sqrt{lambda^2-k^2}y 
c                          exp(i lambda x) sommer_exp(lambda) dlambda 
c                         -----------------------------------
c                              \sqrt{lambda^2-k^2}
c 
c       sommer_exp IS A FUNCTION AT DISCRETE LAMBDA VALUES AND DOES 
c       NOT INVOLVE QUAD WEIGHTS. weights has been left in the 
c       calling sequence for possible use at a later date.
c
c

 
      implicit none
      integer nterms,nquad,iup, ns,i,j
      real *8 center(2),source(2,ns),dipvec(2,ns),pi
      complex *16 chg(ns),dipstr(ns)
      complex *16 xfac(nquad),yfac(nquad),sommer_exp(nquad)
      complex *16 weights(nquad)
      complex *16 zk,zinc,eye,zinc2,zp,zm
      complex *16 zarg1,zarg2,zarg,zdot
c
      pi = 4.0d0*datan(1.0d0)
      eye = (0,1)
      do i = 1,nquad
         sommer_exp(i) = 0.0d0
      enddo
c
      do i = 1,nquad
         do j = 1,ns
            zdot=dipvec(1,j)*eye*xfac(i)+iup*eye*dipvec(2,j)*yfac(i)
            zinc = cdexp(-eye*xfac(i)*(source(1,j)-center(1)))*
     1        cdexp(-iup*eye*yfac(i)*(source(2,j)-center(2)))
ccc            zinc = zinc*weights(i)/(4.0d0*pi)
            sommer_exp(i) = sommer_exp(i) + chg(j)*zinc
            sommer_exp(i) = sommer_exp(i) - dipstr(j)*zinc*zdot
         enddo
      enddo
      return
      end
c
c
c
      subroutine getsommerquad_2d2mi(zk1,zk2,zk,tmax,norder,tol,rx,ry,
     1       nover,ntot,nmax,xnodes,weights,xfac,yfac)
c
c     Discretize Sommerfeld integral to respect two different Helmholtz
c     parameters: zk1, zk2 with real(zk2) > real(zk1), assuming
c     zk is wavenumber of interest (could be zk1 or zk2). 

c     Unlike getsommerquad_2d2, this routine uses an adaptive
c     grid that is refined toward Real(zk1),Real(zk2) and also
c     uses Alpert quadrature rules for 1/sqrt singularity.
c     The *mi* suffix in the name is meant to indicate that 
c     MULTPLE INTERVALS are used.
c     The intervals [-tmax,-Re(zk2)], [-Re(zk2),-Re(zk1)]
c     [-Re(zk1),Re(zk1)], [Re(zk1),Re(zk2)],[Re(zk2),tmax] are
c     all refined toward the singular points.
c             
c     ASSUMES REAL(ZK2) > REAL(ZK1),
c
c     Since the scheme respects the sqrt singularities at 
c     +/- Re(zk1,zk2), it will provide poor precision for other 
c     values of zk.
c
c        
c g_k = 1/(4 pi) \int exp(-\sqrt{lambda^2-zk^2}y exp(i lambda x) dlambda 
c                     -----------------------------------------
c                              \sqrt{lambda^2-zk^2}
c        
c     INPUT:
c
c     zk1     Helmholtz parameter for upper half space
c     zk2     Helmholtz parameter for lower half space
c             with Re(zk2) > Re(zk1).
c     zk      Helmholtz parameter for specific integral being computed
c     tmax    limits on integration in lambda (in evanescent region)
c     norder  order of Alpert rule corrections
c     tol     desired precision
c     rx      max excursion in x
c     ry      max excursion in y
c     nover   oversampling parameter (currently set by user)
c     nmax    dimensions of input arrays xnodes,weights,xfac,yfac
c
c     OUTPUT:
c
c     ntot    total number of quadrature nodes in Sommerfeld integral
c     xnodes  node locations
c     weights quad. weights, includes 1/(4 pi) factor and weight for
c             dlambda but NOT 1/\sqrt{lambda^2-k^2}
c             factor in denominator
c             THUS - ACTUALLY REAL BUT DIMENSIONED AS COMPLEX
c             IN CASE OF FUTURE ADJUSTMENT.
c     xfac    in this case, simply discrete \lambda values,
c             (a copy of xnodes). Keeping the array in case of 
c             some future change of variables
c     yfac    the expression \sqrt{k^2 - lambda^2} at the lambda nodes
c             arranged so that branch cuts are the standard ones
c             (see notes).
c             
c             More precisely, the numerator of the integrand above is
c             cdexp(eye*yfac(i)*y)*cdexp(eye*xfac(i)*x)
c
c     The routines manychgdp_somm2d and somm2d_eval are compatible 
c     with this formalism.
c
      implicit none
      integer ntot,nquad,norder,i,ii,nmax,ntot_loc,iflag,ntemp,nover
      integer nintervals,nminquad
      real *8 xnodes(nmax)
      complex *16 weights(nmax)
      real *8 rkm1,rkp1,tmax,pi,rk,rkx,rky,rkreal,aa,bb
      real *8 rkm2,rkp2,tol,rx,ry
      complex *16 xfac(nmax)
      complex *16 yfac(nmax)
      complex *16 zk1,zk2,eye,zarg1,zarg2,zk
      integer, allocatable :: nquads(:)
      real *8, allocatable :: xnodes_loc(:)
      real *8, allocatable :: weights_loc(:)
      real *8, allocatable :: breaks(:)
c
      allocate(nquads(nmax))
      allocate(breaks(nmax))
      allocate(xnodes_loc(2*nmax))
      allocate(weights_loc(2*nmax))
      nminquad = 24
      eye = dcmplx(0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      rkm1 = -real(zk1)
      rkp1 = real(zk1)
      rkm2 = -real(zk2)
      rkp2 = real(zk2)
c
      aa = -tmax
      bb =  rkm2
      call mkintervalsb(aa,bb,rx,ry,tol,nover,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot_loc)
      write(6,*) ' nintervals is',nintervals
      iflag=2
      call mk_aquad_comp(iflag,breaks,nquads,nintervals,norder,
     1     ntot_loc,xnodes_loc,weights_loc)
      do i = 1,ntot_loc
         xnodes(i) = xnodes_loc(i)
         rk = xnodes(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(i) = zarg1*zarg2
         xfac(i) = rk
         weights(i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntot_loc
c
      aa = rkm2
      bb = rkm1
      call mkintervalsab(aa,bb,rx,ry,tol,nover,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot_loc)
      write(6,*) ' nintervals is',nintervals
      iflag=3
      call mk_aquad_comp(iflag,breaks,nquads,nintervals,norder,
     1     ntot_loc,xnodes_loc,weights_loc)
      do i = 1,ntot_loc
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntemp + ntot_loc
c
      aa = rkm1
      bb = rkp1
      iflag=3
      call mkintervalsab(aa,bb,rx,ry,tol,nover,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot_loc)
      write(6,*) ' nintervals is',nintervals
      call mk_aquad_comp(iflag,breaks,nquads,nintervals,norder,
     1     ntot_loc,xnodes_loc,weights_loc)
      do i = 1,ntot_loc
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntemp + ntot_loc
c
      aa = rkp1
      bb = rkp2
      iflag=3
      call mkintervalsab(aa,bb,rx,ry,tol,nover,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot_loc)
      write(6,*) ' nintervals is',nintervals
      call mk_aquad_comp(iflag,breaks,nquads,nintervals,norder,
     1     ntot_loc,xnodes_loc,weights_loc)
      do i = 1,ntot_loc
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntemp = ntemp + ntot_loc
c 
      aa = rkp2
      bb = tmax
      iflag=1
      call mkintervalsa(aa,bb,rx,ry,tol,nover,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot_loc)
      write(6,*) ' nintervals is',nintervals
      call mk_aquad_comp(iflag,breaks,nquads,nintervals,norder,
     1     ntot_loc,xnodes_loc,weights_loc)
      do i = 1,ntot_loc
         xnodes(ntemp+i) = xnodes_loc(i)
         rk = xnodes_loc(i)
         zarg1 = sqrt(zk-rk)
         zarg2 = sqrt(zk+rk)
         yfac(ntemp+i) = zarg1*zarg2
         xfac(ntemp+i) = rk
         weights(ntemp+i) = eye*weights_loc(i)/(zarg1*zarg2)
         weights(ntemp+i) = weights_loc(i)/(4*pi)
      enddo
      ntot = ntemp + ntot_loc
c
      return
      end
c
c

