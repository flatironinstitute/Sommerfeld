c     Testing code for dielectric interface routines.
c
c     Suppose
c     u_1 = uin_1 + \int exp(-sq1*y)/sq1 exp(i lambda x) sig(lambda)
c                 + \int exp(-sq1*y) exp(i lambda x) mu(lambda)
c     u_2 = uin_2 + \int exp(+sq2*y)/sq2  exp(i lambda x) sig(lambda)
c                 - \int exp(+sq2*y) exp(i lambda x) mu(lambda)
c
c     Continuity conditions are
c
c     [ alpha u] = [ beta dudn] = 0  
c     alpha = (a1,a2) in up/down domain
c     beta  = (b1,b2) in up/down domain
c
c     0 = alpha1 u_1 - alpha2 u_2 = a1 uin_1 - a2 uin_2 + 
c        int exp(i lambda x) [a1 (sig/sq1 + mu) - a2(sig/sq2 - mu)]
c
c     0 = beta1 u_1_y - beta2 u_2_y = b1 uin_1_y - b2 uin_2_y + 
c        int exp(i lambda x) [b1(-sig - sq1*mu) - b2(sig - sq2*mu)]
c
c     f = - (a1 uin_1- a2 uin_2), g = - (b1 uin_1_y- b2 uin_2_y),
c     
c     [f] =  [(a1+a2)  (a1/sq1-a2/sq2) ][mu ]
c     [g]    [(b2sq2-b1sq1)  -(b1+b2)  ][sig] 
c                                            
c     [mu] =   [-(b1+b2)     (a2/sq2-a1/sq1) ] [f]  /det
c     [sig]    [(b1sq1-b2sq2)   (a1+a2)      ] [g]  /det.
c
c      det = -(a1+a2)*(b1+b2) -(a1b2sq2/sq1-a2b2-a1b1+a2b1sq1/sq2) 
c      det = -(a1b2+a2b1) -(a1b2sq2/sq1+a2b1sq1/sq2) 
c
c     When matching modes, be careful to use same scalings for
c     uin_1,uin_2 and layer potentials.
c
c     The sigma integral above is the Fourier transform of the
c     single layer potential with density ( F^{-1} sigma ) or
c     perhaps ( F^{-1} 2 sigma ) and the mu integral is the Fourier 
c     transform of the double layer potential with density 
c     ( F^{-1} mu ) or perhaps ( F^{-1} 2 mu ). Since we don't think
c     about the physical space densities, we can view the above as
c     definitions of the SLP and DLP. 
c
c     This code uses a manufactured solution consisting of the field
c     doe to some sources and dipoles in the lower half plane to define
c     u_1 andnsome sources and dipoles in the upper half plane to 
c     define u_2.
c
      implicit none
      integer nquad,ntot,ifexpon,iup,norder,nover,nsmax
      integer nmax,i,ns,ifgrad,ifhess,nd,nt
      parameter(nsmax= 10)
      parameter(nmax= 100000)
      real *8 xnodes(nmax)
      complex *16 weights1(nmax),weights2(nmax)
      complex *16 sommer_exp1(nmax)
      complex *16 sommer_exp2(nmax)
      complex *16 ff(nmax)
      complex *16 gg(nmax)
      complex *16 sigma(nmax)
      complex *16 rmu(nmax)
      complex *16 xfac1(nmax),yfac1(nmax)
      complex *16 xfac2(nmax),yfac2(nmax)
      complex *16 chg1(nsmax),dipstr1(nsmax)
      complex *16 chg2(nsmax),dipstr2(nsmax),sq1,sq2
      complex *16 a1,a2,b1,b2
      real *8 x,y,z
      real *8 y1,y2,y3,rl1,rl2,rl3
      real *8 source1(2,nsmax),ztrg1(2),dipvec1(2,nsmax),center(2)
      real *8 dipvec2(2,nsmax)
      real *8 source2(2,nsmax),ztrg2(2)
      real *8 pi,d,u,v,rlam,tmax,a,b,t,h
      real *8 rx,ry,tol,rxmax,rymin,rymax,thresh
c
      complex *16 h0,h1,pot,potd,grad(2),hess(3)
      complex *16 gradloc(2),graddir(2),graddir_loc(2)
      complex *16 pot1,pot2,pot3
      complex *16 eye,zsum,zz,zk1,zk2,zsum1,zsum2,det
c
      data eye/(0.0d0,1.0d0)/
c
      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
c
c     Set Helmholtz freq.
c
      zk1 = 50.0d0*2*pi
      zk2 = 100.0d0*2*pi*sqrt(1.0d0+0.001d0*eye)
c
      rlam = 2*pi/dreal(zk2)
      a1 = 2.0d0
      a2 = 1.2d0
      b1 = 10.0d0
      b2 = 3.0d0
      center(1) = 0.0d0
      center(2) = 0.0d0
      source1(1,1) = -1.1d0
      source1(2,1) = -1.2d0
      source1(1,2) = -0.9d0
      source1(2,2) = -1.12d0
      chg1(1) = dcmplx(0.56d0,0.2d0)
      chg1(2) = dcmplx(0.16d0,0.13d0)
      dipstr1(1) = dcmplx(1.16d0,0.1d0)
      dipstr1(2) = dcmplx(2.06d0,0.33d0)
      dipvec1(1,1) = 0.3d0
      dipvec1(2,1) = 0.2d0
      dipvec1(1,2) = 0.3d0
      dipvec1(2,2) = -0.2d0
      ztrg1(1) = 1.0d0
      ztrg1(2) = 1.0d0*rlam
      source2(1,1) = -0.97d0
      source2(2,1) = 1.2d0
      source2(1,2) = 0.7d0
      source2(2,2) = 0.94d0
      chg2(1) = dcmplx(0.256d0,0.42d0)
      chg2(2) = dcmplx(0.46d0,0.313d0)
      dipstr2(1) = dcmplx(1.46d0,0.1d0)
      dipstr2(2) = dcmplx(1.3d0,0.43d0)
      dipvec2(1,1) = -0.4d0
      dipvec2(2,1) = 0.1d0
      dipvec2(1,2) = 0.3d0
      dipvec2(2,2) = 0.2d0
      ztrg2(1) = 1.1d0
      ztrg2(2) = -1.2d0*rlam
      ns = 2
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
      rxmax = 2.2d0
      rymin = 1.9d0
      rymin = rlam
      rymax = 10*rymin

      tol = 1.0d-9
      tmax = 1.3d0*real(zk2) + abs(log(1.0d-6))/rymin + 100
      write(6,*) ' tmax = ',tmax
      norder = 16
      nover = 20
      write(6,*) ' zk1 = ',zk1
      write(6,*) ' zk2 = ',zk2
c
      call getsommerquad_2d2mi(zk1,zk2,zk1,tmax,norder,tol,rxmax,rymax,
     1       nover,ntot,nmax,xnodes,weights1,xfac1,yfac1)
      call getsommerquad_2d2mi(zk1,zk2,zk2,tmax,norder,tol,rxmax,rymax,
     1       nover,ntot,nmax,xnodes,weights2,xfac2,yfac2)
      call prinf(' ntot is *',ntot,1)
c
      write(6,*) ' zk1 = ',zk1
      write(6,*) ' zk2 = ',zk2
      iup = 1
      call manychgdp_somm2d(zk1,iup,ns,source1,chg1,dipstr1,dipvec1,
     1                center,xfac1,yfac1,weights1,sommer_exp1,ntot)
c
      iup = -1
      call manychgdp_somm2d(zk2,iup,ns,source2,chg2,dipstr2,dipvec2,
     1                center,xfac2,yfac2,weights2,sommer_exp2,ntot)
c
      write(6,*) ' zk1 = ',zk1
      write(6,*) ' zk2 = ',zk2
      do i = 1,ntot
         sq1 = -eye*yfac1(i)
         sq2 = -eye*yfac2(i)
         ff(i) = a1*sommer_exp1(i)/sq1 - a2*sommer_exp2(i)/sq2
         gg(i) = -b1*sommer_exp1(i)-b2*sommer_exp2(i)
      enddo
c
c     
c     f = - (a1 uin_1- a2 uin_2), g = - (b1 uin_1_y- b2 uin_2_y),
c
c     [mu] =   [-(b1+b2)     (a2/sq2-a1/sq1) ] [f]  /det
c     [sig]    [(b1sq1-b2sq2)   (a1+a2)      ] [g]  /det.
c      det = -(a1b2+a2b1) -(a1b2sq2/sq1+a2b1sq1/sq2) 
c

      do i = 1,ntot
         sq1 = -eye*yfac1(i)
         sq2 = -eye*yfac2(i)
ccc         write(6,*) ' sq1,sq2',sq1,sq2
         det = -(a1*b2+a2*b1)-(a1*b2*sq2/sq1+a2*b1*sq1/sq2) 
ccc         write(6,*) ' det',det
         rmu(i) = -(b1+b2)*ff(i) + (a2/sq2-a1/sq1)*gg(i)
         sigma(i) = (b1*sq1-b2*sq2)*ff(i) + (a1+a2)*gg(i)
         rmu(i) = rmu(i)/det
         sigma(i) = sigma(i)/det
ccc         rmu(i) = -ff(i)/2/sq1
ccc         sigma(i) = gg(i)/2/sq1
      enddo
c
ccc      write(6,*) ' rmu(1) = ',(rmu(i),i=1,ntot)
      write(6,*) ' zk1 = ',zk1
      write(6,*) ' zk2 = ',zk2
c
c     test at some target points
c
      ifgrad = 1
      ifhess = 0
      nd = 1
      nt = 1
      thresh = 1.0d-12
c
c     test target in upper half space
c
      do i = 1,ntot
         sq1 = -eye*yfac1(i)
         sq2 = -eye*yfac2(i)
ccc         rmu(i) = rmu(i)/sq1
ccc         sigma(i) = sigma(i)/sq1
      enddo
      write(6,*) ' rmu(1) = ',(rmu(i),i=1,6)

      call prin2(' ztrg is *',ztrg1,2)
      iup = 1
      call somm2d_eval(zk1,iup,ztrg1,xnodes,center,xfac1,yfac1,
     1                weights1,sigma,ntot,zsum,gradloc)
      call prin2(' zsum is *',zsum,2)
      zsum1= zsum
      grad(1)= gradloc(1)
      grad(2)= gradloc(2)
      call somm2d_eval_dp(zk1,iup,ztrg1,xnodes,center,xfac1,yfac1,
     1                weights1,rmu,ntot,zsum,gradloc)
      call prin2(' zsum is *',zsum,2)
      zsum1= zsum1+zsum
      grad(1)= grad(1)+gradloc(1)
      grad(2)= grad(2)+gradloc(2)
      call prin2(' zsum1 is *',zsum1,2)
      call prin2(' grad is *',grad,4)
      pot = 0.0d0
      potd = 0.0d0
      graddir_loc(1) = 0.0d0
      graddir_loc(2) = 0.0d0
      call h2d_directcg(nd,zk1,source1,ns,chg1,
     1           ztrg1,nt,pot,graddir_loc,thresh)
      call prin2(' chg pot is *',pot,2)
      call prin2(' ztrg is *',ztrg1,2)
      graddir(1) = graddir_loc(1)
      graddir(2) = graddir_loc(2)
      graddir_loc(1) = 0.0d0
      graddir_loc(2) = 0.0d0
      call h2d_directdg(nd,zk1,source1,ns,dipstr1,dipvec1,
     1           ztrg1,nt,potd,graddir_loc,thresh)
      pot = pot+potd
      graddir(1) = graddir(1)+graddir_loc(1)
      graddir(2) = graddir(2)+graddir_loc(2)
      call prin2(' potd is *',potd,2)
      call prin2(' pot is *',pot,2)
      call prin2(' graddir is *',graddir,4)
      call prin2(' ratio is *',pot/zsum1,2)
      call prin2(' err is *',abs(pot-zsum1),1)
      call prin2(' err1 is *',abs(graddir(1)-grad(1)),1)
      call prin2(' err2 is *',abs(graddir(2)-grad(2)),1)
c
c
c     test target in lower half space
c
      do i = 1,ntot
         sq1 = -eye*yfac1(i)
         sq2 = -eye*yfac2(i)
ccc         rmu(i) = sq1*rmu(i)/sq2
ccc         sigma(i) = sq1*sigma(i)/sq2
      enddo

      call prin2(' ztrg is *',ztrg1,2)
      iup = -1
      call somm2d_eval(zk2,iup,ztrg2,xnodes,center,xfac2,yfac2,
     1                weights2,sigma,ntot,zsum,gradloc)
      call prin2(' zsum is *',zsum,2)
      zsum1= zsum
      grad(1)= gradloc(1)
      grad(2)= gradloc(2)
      call somm2d_eval_dp(zk2,iup,ztrg2,xnodes,center,xfac2,yfac2,
     1                weights2,rmu,ntot,zsum,gradloc)
      call prin2(' zsum is *',zsum,2)
      zsum1= zsum1+zsum
      grad(1)= grad(1)+gradloc(1)
      grad(2)= grad(2)+gradloc(2)
      call prin2(' zsum1 is *',zsum1,2)
      call prin2(' grad is *',grad,4)
      pot = 0.0d0
      potd = 0.0d0
      graddir_loc(1) = 0.0d0
      graddir_loc(2) = 0.0d0
      call h2d_directcg(nd,zk2,source2,ns,chg2,
     1           ztrg2,nt,pot,graddir_loc,thresh)
      graddir(1) = graddir_loc(1)
      graddir(2) = graddir_loc(2)
      graddir_loc(1) = 0.0d0
      graddir_loc(2) = 0.0d0
      call h2d_directdg(nd,zk2,source2,ns,dipstr2,dipvec2,
     1           ztrg2,nt,potd,graddir_loc,thresh)
      pot = pot+potd
      graddir(1) = graddir(1)+graddir_loc(1)
      graddir(2) = graddir(2)+graddir_loc(2)
      call prin2(' potd is *',potd,2)
      call prin2(' pot is *',pot,2)
      call prin2(' graddir is *',graddir,4)
      call prin2(' ratio is *',pot/zsum1,2)
      call prin2(' err is *',abs(pot-zsum1),1)
      call prin2(' err1 is *',abs(graddir(1)-grad(1)),1)
      call prin2(' err2 is *',abs(graddir(2)-grad(2)),1)
      call prinf(' ntot is *',ntot,1)
c     
      stop
      end

