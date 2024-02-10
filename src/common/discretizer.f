c
c     Adaptive discretization of interval [a,b] with refinement to 
c     the left or right boundary or both until the smallest interval
c     is of length less than sqrt(tol). 
c     It then uses an Alpert (Gauss-trapezoidal) rule on the first 
c     and/or last interval for 1/sqrt singularity, if requested.
c
c------------------------------------------------------------------
      subroutine mkintervalsab(a,b,rx,ry,tol,nperwave,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot)
c------------------------------------------------------------------
c
c     Discretize interval [a,b] adaptively toward both endpoints
c     until the smallest interval is of length less than sqrt(tol). 
c     Designed for application to 
c
c     int   exp(-\sqrt(k^2-omega^2)y exp(i k x) /\sqrt(k^2-omega^2)
c
c     INPUT
c     
c     a,b:       interval endpoints
c     rx         maximum excursion in x (target - source)
c     ry         maximum excursion in y (target - source)
c     tol        requested tolerance for adaptive refinement
c     nperwave   number of pts per wavelength away from singularities
c     nminquad   minimum number of quad nodes per subinterval
c     nmax       length of breaks,nquads (and bound on # intervals)
c                CURRENTLY NOT CHECKED.....
c     OUPUT
c
c     nintervals number of subintervals
c     breaks     subinterval endpoints [breaks(i),breaks(i+1)]
c     nquads     number of quadrature nodes on each subinterval
c     ntot       total number fo quadrature nodes
c
      implicit real *8 (a-h,o-z)
      integer nintervals,nmax,nperwave,ntot,nminquad
      integer nquads(nmax)
      real *8 breaks(nmax)
      real *8, allocatable :: btemp(:)
c
c
c
      allocate(btemp(nmax))
      pi = 4.0d0*datan(1.0d0)
      rlen = (b-a)/2.0d0
      write(6,*) ' b-a',b-a
      rmid = (b+a)/2.0d0
      btemp(1) = rmid
      do i = 1,nmax
         rlen = rlen/2.0d0
         btemp(i+1) = btemp(i) - rlen
         ninters = i+1
         if (rlen .lt. sqrt(tol)) then
             btemp(ninters+1) = a
             goto 111
         endif
      enddo
 111  continue
c
      breaks(1) = btemp(ninters+1)
      do i = 1,ninters
         breaks(i+1) = btemp(ninters+1-i)
         rlen = breaks(i+1) - breaks(i)
         nquads(i) = nperwave*sqrt(rlen)*(sqrt(2.0d0)-1.0d0)*ry/(2*pi)
         nquads(i) = nquads(i) + nperwave*rlen*rx/(2*pi)
ccc         write(6,*) ' i, rlen, nq',i,rlen,nquads(i)
         nquads(i) = max(nquads(i),nminquad)
      enddo
      rlen = (b-a)/2.0d0
      do i = 1,ninters-1
         rlen = rlen/2
         breaks(ninters+i+1) = breaks(ninters+i) + rlen
         nquads(ninters+i) = 
     1      nperwave*sqrt(rlen)*(sqrt(2.0d0)-1.0d0)*ry/(2*pi)
         nquads(ninters+i)=nquads(ninters+i)+nperwave*rlen*rx/(2*pi)
ccc         write(6,*) ' ninters+i, rlen, nq',
ccc     1        ninters+i,rlen,nquads(ninters+i)
         nquads(ninters+i) = max(nquads(ninters+i),nminquad)
      enddo
c
c    last interval
c
      breaks(2*ninters+1) = b
      nquads(2*ninters) = nquads(1)
c
      nintervals = 2*ninters
      ntot =0
      do i = 1,nintervals
         ntot = ntot + nquads(i)  
      enddo
c
      return
      end
c
c
c
c
c------------------------------------------------------------------
      subroutine mkintervalsa(a,b,rx,ry,tol,nperwave,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot)
c------------------------------------------------------------------
c
c     Discretize interval [a,b] adaptively only toward left endpoint
c     until the smallest interval is of length less than sqrt(tol). 
c     Designed for application to 
c
c     int   exp(-\sqrt(k^2-omega^2)y exp(i k x) /\sqrt(k^2-omega^2)
c
c     INPUT
c     
c     a,b:       interval endpoints
c     rx         maximum excursion in x (target - source)
c     ry         maximum excursion in y (target - source)
c     tol        requested tolerance for adaptive refinement
c     nperwave   number of pts per wavelength away from singularities
c     nminquad   minimum number of quad nodes per subinterval
c     nmax       length of breaks,nquads (and bound on # intervals)
c                CURRENTLY NOT CHECKED.....
c     OUPUT
c
c     nintervals number of subintervals
c     breaks     subinterval endpoints [breaks(i),breaks(i+1)]
c     nquads     number of quadrature nodes on each subinterval
c     ntot       total number fo quadrature nodes
c
      implicit real *8 (a-h,o-z)
      integer nintervals,nmax,nperwave,ntot,nminquad
      integer nquads(nmax)
      real *8 breaks(nmax)
      real *8, allocatable :: btemp(:)
c
c
c
      allocate(btemp(nmax))
      pi = 4.0d0*datan(1.0d0)
      rlen = (b-a)/2.0d0
      rmid = (b+a)/2.0d0
      btemp(1) = rmid
      do i = 1,nmax
         rlen = rlen/2.0d0
         btemp(i+1) = btemp(i) - rlen
         ninters = i+1
         if (rlen .lt. sqrt(tol)) then
             btemp(ninters+1) = a
             goto 111
         endif
      enddo
 111  continue
c
      breaks(1) = btemp(ninters+1)
      do i = 1,ninters
         breaks(i+1) = btemp(ninters+1-i)
         rlen = breaks(i+1) - breaks(i)
         nquads(i) = nperwave*sqrt(rlen)*(sqrt(2.0d0)-1.0d0)*ry/(2*pi)
         nquads(i) = nquads(i) + nperwave*rlen*rx/(2*pi)
ccc         write(6,*) ' i, rlen, nq',i,rlen,nquads(i)
         nquads(i) = max(nquads(i),nminquad)
      enddo
      rlen = (b-a)/2.0d0
      breaks(ninters+2) = b
      nquads(ninters+1) = 
     1      nperwave*sqrt(rlen)*(sqrt(2.0d0)-1.0d0)*ry/(2*pi)
      nquads(ninters+1) = nquads(ninters+1) + nperwave*rlen*rx/(2*pi)
      nquads(ninters+1) = max(nquads(ninters+1),nminquad)
c
      nintervals = ninters+1
      ntot =0
      do i = 1,nintervals
         ntot = ntot + nquads(i)  
      enddo
c
      return
      end
c
c
c
c------------------------------------------------------------------
      subroutine mkintervalsb(a,b,rx,ry,tol,nperwave,nminquad,nmax,
     1           nintervals,breaks,nquads,ntot)
c------------------------------------------------------------------
c
c     Discretize interval [a,b] adaptively only toward right endpoint
c     until the smallest interval is of length less than sqrt(tol). 
c     Designed for application to 
c
c     int   exp(-\sqrt(k^2-omega^2)y exp(i k x) /\sqrt(k^2-omega^2)
c
c     INPUT
c     
c     a,b:       interval endpoints
c     rx         maximum excursion in x (target - source)
c     ry         maximum excursion in y (target - source)
c     tol        requested tolerance for adaptive refinement
c     nperwave   number of pts per wavelength away from singularities
c     nminquad   minimum number of quad nodes per subinterval
c     nmax       length of breaks,nquads (and bound on # intervals)
c                CURRENTLY NOT CHECKED.....
c     OUPUT
c
c     nintervals number of subintervals
c     breaks     subinterval endpoints [breaks(i),breaks(i+1)]
c     nquads     number of quadrature nodes on each subinterval
c     ntot       total number fo quadrature nodes
c

      implicit real *8 (a-h,o-z)
      integer nintervals,nmax,nperwave,ntot,nminquad
      integer nquads(nmax)
      real *8 breaks(nmax)
c
c
c
      pi = 4.0d0*datan(1.0d0)
      rlen = (b-a)/2.0d0
      rmid = (b+a)/2.0d0
      breaks(1) = a
      breaks(2) = rmid
      do i = 2,nmax
         rlen = rlen/2.0d0
         breaks(i+1) = breaks(i) + rlen
      write(6,*) 'rlen ',rlen
         nintervals = i
         if (rlen .lt. sqrt(tol)) then
             nintervals=nintervals+1
             breaks(nintervals+1) = b
             goto 111
         endif
      enddo
 111  continue
c
      do i = 1,nintervals
         rlen = breaks(i+1) - breaks(i)
         nquads(i) = 
     1      nperwave*sqrt(rlen)*(sqrt(2.0d0)-1.0d0)*ry/(2*pi)
         nquads(i)=nquads(i)+nperwave*rlen*rx/(2*pi)
         nquads(i) = max(nquads(i),nminquad)
      enddo
c
      ntot =0
      do i = 1,nintervals
         ntot = ntot + nquads(i)  
      enddo
c
      return
      end
c
c
c
c
c
c------------------------------------------------------------------
      subroutine mk_aquad_comp(iflag,breaks,nquads,nintervals,norder,
     1     ntot,xnodes,weights)
c------------------------------------------------------------------
c
c     Apply composite Gauss-trapezoidal quadrature, with endpoint
c     corrections if requested.
c
c     INPUT
c     
c     iflag:     0  apply smooth rule
c                1  apply 1/sqrt corrections at left
c                2  apply 1/sqrt corrections at right
c                3  apply 1/sqrt corrections at both
c     breaks     subintervla endpoints [breaks(i),breaks(i+1)]
c     nquads     number of quadrature nodes on each subinterval
c     nintervals number of subintervals
c     norder     order of Alpert correction. nquads(i) must be at
c                least 2*norder.
c     ntot       total number fo quadrature nodes
c
c     OUPUT
c
c     xnodes     quadrature nodes, ordered left to right.
c     weights    quadrature weights.
c
      implicit real *8 (a-h,o-z)
      integer nintervals,norder,nqstep
      integer nquads(nintervals)
      real *8 breaks(nintervals+1)
      real *8 xnodes(ntot)
      real *8 weights(ntot)
c
      nqstep = 0
c
      if ((iflag.eq.1) .or. (iflag.eq.3)) then
         call alpert_sqrta(breaks(1),breaks(2),nquads(1),norder,
     1     xnodes(1),weights(1))
      else
         call alpert_smooth(breaks(1),breaks(2),nquads(1),norder,
     1     xnodes(1),weights(1))
      endif
      nqstep = nqstep + nquads(1)
c
      do i = 2,nintervals-1
         call alpert_smooth(breaks(i),breaks(i+1),nquads(i),norder,
     1     xnodes(nqstep+1),weights(nqstep+1))
         nqstep = nqstep + nquads(i)
      enddo
      if ((iflag.eq.2) .or. (iflag.eq.3)) then
         call alpert_sqrtb(breaks(nintervals),breaks(nintervals+1),
     1   nquads(nintervals),norder,xnodes(nqstep+1),weights(nqstep+1))
      else
         call alpert_smooth(breaks(nintervals),breaks(nintervals+1),
     1   nquads(nintervals),norder,xnodes(nqstep+1),weights(nqstep+1))
      endif
      nqstep = nqstep + nquads(nintervals)

      write(6,*) 'nqstep ',nqstep
c
      return
      end
c
