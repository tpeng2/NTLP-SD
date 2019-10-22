      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)
c
c -------- get 2d fft using fftpack routines and parallel mpi
c          use fftpack storage a0, (a1,b1), (a2,b2),...,
c
c         isgn = -1 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn = -2 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
c
c         isgn =  1 do inverse transform, move to physical space
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn =  2 do inverse transform, move to physical space
c                   incoming array is at(ny,jxs:jxe,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
      real ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2),
     +     trigx(2*nx+15), trigc(4*ny+15),
     +     a2d(2,ny), a_wrk(nx)
      integer jx_s(0:np-1), jx_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      nxp2 = nx + 2
      if(isgn .lt. 0) then
         fn   = 1.0/(float(nx)*float(ny))
c
c ------ 1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               do ix=1,nx
                  a_wrk(ix) = ax(ix,iy,iz)*fn
               enddo
               call rfftf(nx,a_wrk(1),trigx(1))
               ax(1,iy,iz) = a_wrk(1)
               ax(2,iy,iz) = 0.0
               do ix=2,nx
                  ax(ix+1,iy,iz) = a_wrk(ix)
               enddo
               ax(nx+2,iy,iz) = 0.0
            enddo
         enddo
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftf(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
c
c ---- decide whether to transpose back or leave as is
c
         if(isgn .eq. -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
      else
c
c ---- decide whether to first transpose or leave as is
c
         if(isgn .eq. 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftb(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------  1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               a_wrk(1) = ax(1,iy,iz)
               do ix=2,nx
                  a_wrk(ix) = ax(ix+1,iy,iz)
               enddo
               call rfftb(nx,a_wrk(1),trigx(1))
               do ix=1,nx
                  ax(ix,iy,iz) = a_wrk(ix)
               enddo
            enddo
         enddo
      endif
c
      return
      end
        subroutine gauss_newton_2d(vnext,h,vec1,vec2,flag)
        use particles
        implicit none

        real, intent(in) :: vnext(3), h, vec1(2)
        real, intent(out) :: vec2(2)
        integer, intent(out) :: flag
        real :: error = 1E-8, fv1(2), fv2(2), v1(2), v_output(3), rel
        real :: diff, temp1(2), temp2(2), relax, coeff, correct(2)
        real, dimension(1:2, 1:2) :: J, fancy, inv, finalJ
        integer :: iterations, neg, counts

        iterations = 0
        flag = 0

        v1 = vec1
        fv2 = (/1., 1./)
        coeff = 0.1
        do while ((sqrt(dot_product(fv2, fv2)) > error).AND.
     +              (iterations<1000))


                !print*, 'conditions: ', part%Tf, part%qinf, 'error', 
                !+          sqrt(dot_product(fv2, fv2)), 'vel', part%uf

                iterations = iterations + 1
                !print*, 'vnext', vnext, 'v1', v1, 'v_output', v_output,
                !+          'h', h

                call ie_vrt_nd(vnext, v1(1), v1(2), v_output, fv1, h)
                !print*,'fv1', fv1
                call jacob_approx_2d(vnext, v1(1), v1(2), h, J)

                !print*,'error ', sqrt(dot_product(fv1, fv1))

                fancy = matmul(transpose(J), J)

                call inverse_finder_2d(fancy, inv)

                finalJ = matmul(inv, transpose(J))

                !print*,'corrective term', matmul(finalJ, fv1)
                correct = matmul(finalJ,fv1)
                vec2 = v1 - correct

                call ie_vrt_nd(vnext, v1(1), v1(2),v_output,temp1,h)
                call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,temp2,h)

                diff = sqrt(dot_product(temp1,temp1))-
     +                 sqrt(dot_product(temp2,temp2))

                if (sqrt(dot_product(correct,correct))<1E-8) then
                        !flag = 1
                        EXIT
                end if

                relax = 1.0
                counts = 0

               do while ((diff<0).OR.(vec2(1)<0).OR.(vec2(2)<0)
     +                   .OR.isnan(vec2(1)))
                        counts = counts + 1
                        coeff = 0.5
                        relax = relax * coeff
                        vec2 = v1-matmul(finalJ,fv1)*relax
                call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,temp2,h)
                        diff = sqrt(dot_product(temp1,temp1))-
     +                  sqrt(dot_product(temp2,temp2))

                        if (counts>10) EXIT
                end do

                v1 = vec2

                call ie_vrt_nd(vnext, vec2(1), vec2(2), v_output, fv2,h)
        end do
      if (iterations == 100) flag = 1
      if (isnan(vec2(1)) .OR. vec2(1)<0 .OR. isnan(vec2(2))
     + .OR. vec2(2)<0) flag = 1

      end subroutine gauss_newton_2d
      subroutine ie_vrt_nd(vnext, tempr, tempt, v_output,rt_output, h)
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      real, intent(in) :: vnext(3), tempr, tempt, h
      real, intent(out) :: v_output(3), rT_output(2)

      real :: esa, dnext,  m_w, rhop, Rep, taup,
     +         vprime(3), rprime, Tprime, qstr, Shp, Nup, dp, VolP,
     +         diff(3), diffnorm, Tnext, rnext, T
      real :: taup0, g(3), mod_Magnus


        taup0 = (((part%m_s)/((2./3.)*pi2*radius_init**3) + rhow) *
     +   (radius_init*2)**2)/(18*rhoa*nuf)
        g(1:3) = part_grav(1:3)

        ! quantities come in already non-dimensionalized, so must be
        ! converted back;
        ! velocity is not non-dimensionalized so no need to change
        rnext = tempr * part%radius
        Tnext = tempt * part%Tp
        dnext = rnext * 2.

        esa = mod_Magnus(part%Tf)
        VolP = (2./3.)*pi2*rnext**3
        rhop = (part%m_s + VolP*rhow) / VolP

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Velocity !!!
        diff(1:3) = part%uf - vnext
        diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
        Rep = dnext * diffnorm/nuf
        taup = (rhop * dnext**2)/(18.0*rhoa*nuf)
        vprime(1:3) = (1. + 0.15 * (Rep**0.687)) * (1./taup) *
     +  diff(1:3) - g(1:3)
        vprime(1:3) = vprime(1:3) * taup0 ** 2
        !!!!!!!!!!!!!!!!

        !!! Humidity !!!
        qstr = (Mw/(Ru*Tnext*rhoa)) * esa *
     +  exp(((Lv*Mw/Ru)*((1./part%Tf) - (1./Tnext))) +
     +  ((2.*Mw*Gam)/(Ru*rhow*rnext*Tnext)) -
     +  ((Ion*Os*part%m_s*(Mw/Ms))/(Volp*rhop-part%m_s)))
        !!!!!!!!!!!!!!!!!!

        !!! Radius !!!
        Shp = 2. + 0.6 * Rep**(1./2.) * Sc**(1./3.)
        rprime = (1./9.) * (Shp/Sc) * (rhop/rhow) *
     +  (rnext/taup) * (part%qinf - qstr)
        rprime = rprime * (taup0/part%radius)
        !!!!!!!!!!!!!!!!!

        !!! Temperature !!!
        Nup = 2. + 0.6*Rep**(1./2.)*Pra**(1./3.);

        Tprime = -(1./3.)*(Nup/Pra)*CpaCpp*(rhop/rhow)*
     +           (1./taup)*(Tnext-part%Tf) +
     +           3.*Lv*(1./(rnext*Cpp))*rprime*(part%radius/taup0)
        Tprime = Tprime * (taup0/part%Tp)
        !!!!!!!!!!!!!!!!!

        ! velocity is not non-dimensionalized so it does not need to be
        ! changed back
        v_output(1:3) = vnext(1:3) - part%vp(1:3) - h * vprime(1:3)
        rT_output(1) = rnext/part%radius - 1.0  - h*rprime
        rT_output(2) = Tnext/part%Tp - 1.0  - h*Tprime

      end subroutine ie_vrt_nd
      subroutine inverse_finder_2d(C, invC)
        implicit none
        real :: det
        real, dimension(1:2, 1:2), intent(in) :: C
        real, dimension(1:2, 1:2), intent(out) :: invC

        det = C(1, 1) * C(2, 2) - C(1, 2) * C(2, 1)

        invC = reshape((/C(2, 2), -C(2,1), -C(1, 2), C(1, 1)/),
     +   shape(invC))
        invC = (1./det)*invC

      end subroutine inverse_finder_2d
        subroutine jacob_approx_2d(vnext, rnext, tnext, h, J)
        use particles
        implicit none
        integer :: n

        real, intent(in) :: vnext(3), rnext, tnext, h
        real, intent(out), dimension(1:2, 1:2) :: J
        real :: diff = 0, v_output(3), rt_output(2),xper(2),
     +          fxper(2), ynext(2),xper2(2),fxper2(2)

        diff = 1E-12

        ynext(1) = rnext
        ynext(2) = tnext

        call ie_vrt_nd(vnext, rnext, tnext, v_output, rt_output, h)

        xper = ynext
        xper2 = ynext

        do n=1, 2
                xper(n) = xper(n) + diff
                xper2(n) = xper2(n) - diff
                call ie_vrt_nd(vnext, xper(1), xper(2),v_output,
     +           fxper,h)
                call ie_vrt_nd(vnext, xper2(1), xper2(2),v_output,
     +           fxper2,h)
                !              if (.NOT. 
                !+
                !sqrt(dot_product((fxper-rt_output),(fxper-rt_output)))
                !> 0)
                J(:, n) = (fxper-rt_output)/diff
                xper(n) = ynext(n)
                xper2(n) = ynext(n)
        end do

      end subroutine jacob_approx_2d
      subroutine lterp(n,zary,zpt,i,ip1,ratio)
c
c ---- linear interpolation for zpt in zary, where zary is 
c      monotonic increasing or decreasing function
c
      dimension zary(*)
      nm1 = n-1
      if(n.le.1) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
      endif
      if(zary(1) .lt. zary(2)) go to 1
                               go to 101
    1 continue
c
c **** monotonic increasing array
c
        if(zpt .lt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .gt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .ge. zary(j) .and.
     $           zpt .le. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
                 go to 999
              endif
        enddo
c
c **** decreasing array
c
  101 continue 
        if(zpt .gt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .lt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .le. zary(j) .and.
     $           zpt .ge. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
                 go to 999
              endif
        enddo
  999 continue
      return
      end
      function rlim(d1,d2,d3)
c
c ------------- Cees's kappa=1/3 scheme
c
      r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
      rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c
c ------------- Cees's kappa=-1 scheme
c
c     r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
c     rlim = (d2-d3)*amin1(abs(r),0.5)
c
c ------------- first order upwind
c
c     rlim = 0.0
c
c ------------- QUICK scheme
c
c     rlim = -0.25*d2 - 0.125*d3 + 0.375*d1
c
      return
      end
