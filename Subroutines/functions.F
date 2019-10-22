! NOTES:
! function ran1(idum)
! function ran2(idum)
! function gasdev(idum)
! function mod_Magnus
      function ran1(idum)
c
c ----------- stolen from numerical recipes,p. 271
c
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,
     +           ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      if(idum .le. 0 .or. iy .eq. 0) then
         idum = max(-idum,1)
         do j=ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum - k*iq) - ir*k
            if(idum .lt. 0) idum = idum + im
            if(j .le. ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k     = idum/iq
      idum  = ia*(idum - k*iq) - ir*k
      if(idum .lt. 0) idum = idum + im
      j     = 1 + iy/ndiv
      iy    = iv(j)
      iv(j) = idum
      ran1  = min(am*iy, rnmx)
c
      return
      end
      function ranf()
      data inc /1/
      save inc, ix, ia, m, fm
      if(inc.eq.1) then
        inc = 2
        m = 2**20
        fm = float(m)
        ix = 566387
        ia = 2**10 + 3
      endif
      ix = mod(ia*ix,m)
      fx = float(ix)
      ranf = fx/fm
      return
      end

      function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/,iv/NTAB*0/,iy/0/

      if (idum .le. 0) then
          idum=max(-idum,1)
          idum2 = idum
          do j = NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum .lt. 0) idum=idum+IM1
             if (j .le. NTAB) iv(j) = idum
          end do
          iy=iv(1)
      end if
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum .lt. 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 .lt. 0) idum2=idum2+IM2
      j = 1+iy/NDIV
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy .lt. 1) iy = iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2

      function gasdev(idum)
       integer :: idum
       real :: gasdev
       integer :: iset
       real :: fac,gset,rsq,v1,v2,ran2
       save iset,gset
       data iset/0/

       if (idum .lt. 0) iset=0
       if (iset .eq. 0) then
 1000    v1 = 2.*ran2(idum)-1.0
         v2 = 2.*ran2(idum)-1.0
         rsq = v1**2+v2**2
         if ( (rsq .ge. 1) .or. (rsq .eq. 0)) goto 1000
         fac = sqrt(-2.0*log(rsq)/rsq)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
       else
         gasdev = gset
         iset = 0
       end if
       return
      end
      function mod_Magnus(T)
      implicit none

      !Take in T in Kelvin and return saturation vapor pressure using
      !function of Alduchov and Eskridge, 1996
      real,intent(in) :: T  
      real :: mod_Magnus

      
      mod_Magnus = 610.94 *exp((17.6257*(T-273.15))/
     +  (243.04+(T-273.15)))
       
      
      end function mod_Magnus
