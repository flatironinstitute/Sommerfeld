c
c       using Alpert quadrature to evaluate Sommerfeld integral
c
c       We assume zk1,zk2 are Helmholtz parameters for upper and
c       lower domains, respectively with real parts rk1 < rk2.
c
c        
c g_k = 1/(4 pi) \int exp(-\sqrt{lambda^2-k^2}y exp(i lambda x) dlambda 
c                     -----------------------------------------
c                              \sqrt{lambda^2-k^2}
c        
c       ------x--------x----.----x--------x----------
c             |        |         |        |        lambda axis
c            -rk2k   -rk1       +rk1     +rk2
c
c     Use one-sided inverse sqrt quadrature on the two semi-infinite 
c     intervals and two sided inverse sqrt quadrature for the three 
c     finite segments. The idea is to have a single set of discrete 
c     nodes in lambdathat can be uses to compute the Sommerfeld 
c     representation for both sides.
c
      implicit none
      integer nquad,ntot,ifexpon,iup,norder,nover
      integer nmax,i,ns,ifgrad,ifhess,nd,nt
      parameter(nmax= 100000)
      real *8 xnodes(nmax)
      complex *16 weights(nmax)
      complex *16 sommer_exp(nmax)
      complex *16 xfac(nmax)
      complex *16 yfac(nmax)
      complex *16 chg,dipstr
      real *8 x,y,z,proj,rr
      real *8 y1,y2,y3,rl1,rl2,rl3
      real *8 source(2),ztrg(2),dipvec(2),center(2)
      real *8 pi,d,u,v,rlam,tmax,a,b,t,h
      real *8 rx,ry,tol,rxmax,rymin,rymax,thresh
c
      complex *16 h0,h1,pot,grad(2),hess(3)
      complex *16 pot1,pot2,pot3
      complex *16 eye,zsum,zz,zk1,zk2,zk
c
      data eye/(0.0d0,1.0d0)/
c
      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
      zk1 = 100.0d0*pi
      zk2 = 200.0d0*2*pi*sqrt(1.0d0+0.001d0*eye)
c
      center(1) = 0.0d0
      center(2) = 0.0d0
      source(1) = -1.0d0
      source(2) = 0.0d0
ccc      ztrg(2) = 10*2.0d0*pi/dreal(zk2)
ccc      x = ztrg(1) - source(1)
ccc      y = ztrg(2) - source(2)
      zk = zk2
c
c     rxmax is max excursion in x
c     rymin is min separation of source and target in y
c     rymax is max separation of source and target in y
c
c     rxmax,rymax determine how oscillatory the integral is
c     rymin controls rate of decay of integrand in Sommerfeld integral
c
c     currently set in terms of zk2. May need more careful thinking...
c
      rxmax = 2.0d0
      rymin = 2.0d0*pi/dreal(zk2)/2.0d0
      rymax = 10.0d0*2.0d0*pi/dreal(zk2)

      call prin2(' rxmax is *',rxmax,1)
      call prin2(' rymin is *',rymin,1)
      call prin2(' rymax is *',rymax,1)
c
c     exp(- \sqrt{\lambda^2- k^2} y = tol ->
c     once \lambda > k -> exponential decay of order exp(- \lambda y)
c     \lambda approx (1+C)*real(zk2) + log(tol)/rymin
c
      tol = 1.0d-4
      tmax = 1.3d0*real(zk2) + abs(log(tol))/rymin
      write(6,*) ' tmax = ',tmax
      norder = 16
      nover = 6
c
      call getsommerquad_2d2(zk1,zk2,zk,tmax,norder,tol,rxmax,rymax,
     1       nover,ntot,nmax,xnodes,weights,xfac,yfac)
      call prinf(' ntot is *',ntot,1)
ccc      call prin2(' xnodes is *',xnodes,ntot)
ccc      call prin2(' weights is *',weights,ntot)
ccc      open(unit=11,file='xnodes.m',status='unknown')
ccc         write(11,*) ' xn = ['
ccc      do i = 1,ntot
ccc         write(11,*) xnodes(i)
ccc      enddo
ccc         write(11,*) '];'

      chg = 1.0d0
      dipstr = 1.0d0
      dipvec(1) = 0.0d0
      dipvec(2) = 1.0d0
      iup = -1
      call onechgdp_somm2d(zk,iup,source,chg,dipstr,dipvec,
     1                center,xfac,yfac,weights,sommer_exp,ntot)

c
c
c     test at some target points
c
      ifexpon=1
      ifgrad = 0
      ifhess = 0
      ns = 1
      nd = 1
      nt = 1
      thresh = 1.0d-12
c
      ztrg(1) = 1.0d0
      ztrg(2) = 2.0d0*pi/dreal(zk2)/2.0d0
      ztrg(2) = -rymax
      call prin2(' ztrg is *',ztrg,2)
      call somm2d_eval(zk,iup,ztrg,xnodes,center,xfac,yfac,
     1                weights,sommer_exp,ntot,zsum,grad)
      call prin2(' rymax = 10 wavelength zsum is *',zsum,2)
ccc      call hpotgrad2dall(ifgrad,ifhess,source,chg,ns,
ccc     1 ztrg,zk,pot,grad,hess)
      y1 = ztrg(2)
      pot1 = zsum
      call h2d_directcp(nd,zk,source,ns,chg,
     $           ztrg,nt,pot,thresh)
      call h2d_directdp(nd,zk,source,ns,dipstr,dipvec,
     $           ztrg,nt,pot,thresh)

      call prin2(' pot is *',pot,2)
      call prin2(' ratio is *',pot/zsum,2)
      call prin2(' err is *',abs(pot-zsum),1)
c
c
      ztrg(2) = -2.0d0*pi/dreal(zk2)/2.0d0
      call prin2(' ztrg is *',ztrg,2)
      call somm2d_eval(zk,iup,ztrg,xnodes,center,xfac,yfac,
     1                weights,sommer_exp,ntot,zsum,grad)
      call prin2(' 1/2 wave zsum is *',zsum,2)
      y2 = ztrg(2)
      pot2 = zsum
ccc      call hpotgrad2dall(ifgrad,ifhess,source,chg,ns,
ccc     1 ztrg,zk,pot,grad,hess)
      pot = 0.0d0
      call h2d_directcp(nd,zk,source,ns,chg,
     $           ztrg,nt,pot,thresh)
      call h2d_directdp(nd,zk,source,ns,dipstr,dipvec,
     $           ztrg,nt,pot,thresh)
      call prin2(' pot is *',pot,2)
      call prin2(' ratio is *',pot/zsum,2)
      call prin2(' err is *',abs(pot-zsum),1)
c
      ztrg(2) = -2.0d0*pi/dreal(zk2)/4.0d0
      call prin2(' ztrg is *',ztrg,2)
      call somm2d_eval(zk,iup,ztrg,xnodes,center,xfac,yfac,
     1                weights,sommer_exp,ntot,zsum,grad)
      call prin2(' 1/4 wave zsum is *',zsum,2)
ccc      call hpotgrad2dall(ifgrad,ifhess,source,chg,ns,
ccc     1 ztrg,zk,pot,grad,hess)
      pot = 0.0d0
      call h2d_directcp(nd,zk,source,ns,chg,
     $           ztrg,nt,pot,thresh)
      call h2d_directdp(nd,zk,source,ns,dipstr,dipvec,
     $           ztrg,nt,pot,thresh)
      call prin2(' pot is *',pot,2)
      call prin2(' ratio is *',pot/zsum,2)
      call prin2(' err is *',abs(pot-zsum),1)
      stop
      y3 = ztrg(2)
      pot3 = zsum
      pot = 0.0d0
c
      ztrg(2) = 2.0d0*pi/dreal(zk2)/1000.0d0
      rl1 = (ztrg(2) - y2)*(ztrg(2)-y3)/((y1-y2)*(y1-y3))
      rl2 = (ztrg(2) - y1)*(ztrg(2)-y3)/((y2-y1)*(y2-y3))
      rl3 = (ztrg(2) - y1)*(ztrg(2)-y2)/((y3-y1)*(y3-y2))
      pot3 = rl1*pot1 + rl2*pot2 + rl3*pot3
      call prin2(' pot3 is *',pot3,2)
      call h2d_directcp(nd,zk,source,ns,chg,
     $           ztrg,nt,pot,thresh)
      call h2d_directdp(nd,zk,source,ns,dipstr,dipvec,
     $           ztrg,nt,pot,thresh)
      call prin2(' pot is *',pot,2)
      call prin2(' ratio is *',pot/pot3,2)
      call prin2(' err is *',abs(pot-pot3),1)
c
      stop
      end

