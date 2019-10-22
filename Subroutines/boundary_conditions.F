      subroutine upper
c
c ---- set boundary condition on upper boundary iz=nnz
c      option for special radiation boundary condition
c                 index f(.,.,1)  indicates upper. 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      iz   = nnz
      izm1 = iz - 1
      izm2 = iz - 2
      izp1 = iz + 1
      izp2 = iz + 2
c
      if(ibcu .eq. 0) then
c
c --------- boundary conditions are gradient conditions
c
c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = dtdzf
c            wbc    = 0.0
c            ebc    = 0.0
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- first get average scalar gradient
c
           dtdzf(iscl) = 0.0
           do iy=iys,iye
           do ix=1,nnx
              dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +                      t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
           enddo
           enddo
           dtdzf(iscl) = dtdzf(iscl)*fnxy
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
      else if(ibcu .eq. 1) then
c
c ------------- special if iradup boundary condition
c               get estimate of w from continuity and 
c               linearized relation for pressure
c
      xmeanp = 0.0
      grad_ug = ug(nnz) - ug((nnz-1))
      do iy=iys,iye
      do ix=1,nnx
         wbc(ix,iy,1) = w(ix,iy,izm1)-
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz))*dzw(iz)
         pbc(ix,iy,1) = .5*(w(ix,iy,izm1)+wbc(ix,iy,1))
         ebc(ix,iy,1) = 0.0
         ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
         vbc(ix,iy,1) = v(ix,iy,iz)
         pbc2(ix,iy,1)=0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +
     +              0.25*(w(ix,iy,izm1)**2 + wbc(ix,iy,1)**2)
         xmeanp = xmeanp + pbc2(ix,iy,1)
      enddo
      enddo
      call mpi_sum_xy(xmeanp,myid,iss,ise,1)
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                          dtdzf(iscl)*dzu(nnzp1)
         enddo
         enddo
      enddo
      xmeanp = xmeanp*fnxy
      do iy=iys,iye
      do ix=1,nnx
         pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
      enddo
      enddo
c
c ---------- end if block
c
      endif
c
      do iy=iys,iye
      do ix=1,nnx
         w(ix,iy,iz)   = wbc(ix,iy,1)
         e(ix,iy,iz)   = ebc(ix,iy,1)
         r3(ix,iy,iz)  = 0.0
         r5(ix,iy,iz)  = 0.0
         u(ix,iy,izp1) = ubc(ix,iy,1)
         v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
         w(ix,iy,izp1) = wbc(ix,iy,1)
         e(ix,iy,izp1) = ebc(ix,iy,1)
         r3(ix,iy,izp1)= 0.0
         r5(ix,iy,izp1)= 0.0
c
c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc
c
         wx(ix,iy,izp1) = 0.0
         wy(ix,iy,izp1) = 0.0
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
c        ux(ix,iy,izp1) = ubc(ix,iy,1)
c        uy(ix,iy,izp1) = ubc(ix,iy,1)
c        vx(ix,iy,izp1) = vbc(ix,iy,1)
c        vy(ix,iy,izp1) = vbc(ix,iy,1)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
            t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine lower_dns

      !Make lower BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izm1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(1)

c ------------ initialize u, v, w, t and derivatives at izm1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izm1)  = -2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izm1)  = -v(ix,iy,iz) 
         w(ix,iy,izm1)  = 0.0 
         r3(ix,iy,izm1) = 0.0
         e(ix,iy,izm1)  = 0.0 
         ux(ix,iy,izm1) = 0.0
         uy(ix,iy,izm1) = 0.0
         vx(ix,iy,izm1) = 0.0
         vy(ix,iy,izm1) = 0.0
         wx(ix,iy,izm1) = 0.0 
         wy(ix,iy,izm1) = 0.0 
           wbc(ix,iy,2) = 0.0
           ebc(ix,iy,2) = 0.0
           ubc(ix,iy,2) = u(ix,iy,izm1)
           vbc(ix,iy,2) = v(ix,iy,izm1)
           pbc(ix,iy,2) = 0.0
           pbc2(ix,iy,2)= 0.0
         
      enddo
      enddo


!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)
! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Tbot)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izm1) = 2.0*Tbot(iscl) - t(ix,iy,iscl,iz) 
            end if

            if (isfc(iscl)==0) then 
            flux = dzu(0)*wtsfc(iscl)/vis_s(ix,iy,iscl,izm1)
            t(ix,iy,iscl,izm1) = t(ix,iy,iscl,iz)     
     +       + flux  
            end if
         enddo
         enddo
      enddo

      end subroutine lower_dns
      subroutine upper_dns

      !Make upper BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izp1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = nnz
      izp1 = iz + 1

c
c ------------ initialize u, v, w, t and derivatives at izp1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izp1)  = 2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izp1)  = -v(ix,iy,iz) 
         w(ix,iy,izp1)  = 0.0 
         r3(ix,iy,izp1) = 0.0
         e(ix,iy,izp1)  = 0.0 
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
         wx(ix,iy,izp1) = 0.0 
         wy(ix,iy,izp1) = 0.0 
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,izp1)
           vbc(ix,iy,1) = v(ix,iy,izp1)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
      enddo
      enddo

! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Ttop)

!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izp1) = 2.0*Ttop(iscl)-t(ix,iy,iscl,iz) 
            end if
            if (isfc(iscl)==0) then 
            flux = dzu(izp1)*wtsfc(iscl)/vis_s(ix,iy,iscl,izp1)
            t(ix,iy,iscl,izp1) = t(ix,iy,iscl,iz)      
     +       - flux 
            end if
         enddo
         enddo
      enddo

      end subroutine upper_dns
      subroutine lower_free(it)
c
c --------------- setup lower boundary condition for free
c                 convection where each processor applies
c                 log-law at several (ix,iy) for iz = 1
c
c                 index f(.,.,2)  indicates lower
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
      real sbuf(2+2*nscl,mxs:mxe,iys:iye)
      real rbuf((2+2*nscl)*nnx*(iye+1-iys))
c
c -------------- broadcast level 1 data everywhere
c
      if(iss .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,1) = u(ix,iy,1)
            u_level1(ix,iy,2) = v(ix,iy,1)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
         enddo
         enddo
         enddo
      endif
      num = nnx*(iye + 1 - iys)*(2+nscl)
c
c ------ send all of root data to other processors
c
      call mpi_send_root(u_level1(1,iys,1),
     +             num,myid,numprocs,ncpu_s)
c
c --------- every task gets their own fluxes and surface scalars
c
      call suft2(u_level1,it)
c
c --------- send surface scalars and momentum fluxes
c           back to root(s)
c
      if(numprocs .eq. 1) go to 999
c
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(1,ix,iy)  = tau13m(ix,iy)
         sbuf(2,ix,iy)  = tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
         sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
      enddo
      enddo
      enddo
c
      irow_r = mod(myid,ncpu_s)
      if(myid .ge. ncpu_s) then
        num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
        call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ncpu_s,numprocs-1,ncpu_s
           num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
           call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
c          call f_suft2(rbuf,maxnx,maxny,mx_s(l),mx_e(l),iys,iye,nscl,
           call f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
        enddo
      endif
c
  999 continue
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
      if(iss .eq. 0) then
c
         buf(1) = 0.0
         buf(2) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(1) = buf(1) + tau13m(ix,iy)
            buf(2) = buf(2) + tau23m(ix,iy)
         enddo
         enddo
         do iscl=1,nscl
            buf(2+iscl)      = 0.
            buf(2+nscl+iscl) = 0.
            do iy=iys,iye
            do ix=1,nnx
               buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
               buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
            enddo
            enddo
         enddo
c
         call mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
         uwsfc = buf(1)*fnxy
         vwsfc = buf(2)*fnxy
         do iscl=1,nscl
            wtsfc(iscl) = buf(2+iscl)*fnxy
            tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
         enddo
c
         iz   = 1
         izm1 = iz - 1
         dz_i = dzu_i(iz)
c
         do iy=iys,iye
         do ix=1,nnx
            ebc(ix,iy,2)=amax1(e(ix,iy,iz),sml_eg)
            wbc(ix,iy,2)= 0.0
            pbc(ix,iy,2) = 0.0
            pbc2(ix,iy,2) = 0.0
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            dudz     = 2.*u(ix,iy,iz)*dz_i
            dvdz     = 2.*v(ix,iy,iz)*dz_i
            ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
            vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
         enddo
         enddo
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
               tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
            enddo
            enddo
         enddo
c
c ------------ initialize u, v, w, t and derivatives at izm1
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1)  = ubc(ix,iy,2)
            v(ix,iy,izm1)  = vbc(ix,iy,2)
            w(ix,iy,izm1)  = wbc(ix,iy,2)
            r3(ix,iy,izm1) =  0.0
            e(ix,iy,izm1)  = ebc(ix,iy,2)
            ux(ix,iy,izm1) = 0.0
            uy(ix,iy,izm1) = 0.0
            vx(ix,iy,izm1) = 0.0
            vy(ix,iy,izm1) = 0.0
            wx(ix,iy,izm1) = wbc(ix,iy,2)
            wy(ix,iy,izm1) = wbc(ix,iy,2)
         enddo
         enddo
c
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
            enddo
            enddo
         enddo
c
c ----- end of if block for root row
c
      endif
c
 7999 continue
c
      return
      end
      subroutine lower(it)
c
c ------ setup lower boundary condition for entire plane at (iz = 1)
c        using either businger or large formulas with wind.
c        index f(.,.,2)  indicates lower. threaded version
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      real sfc_flx(2+nscl)
c
      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(1)
c
      do iy=iys,iye
      do ix=1,nnx
         ebc(ix,iy,2)  = amax1(e(ix,iy,iz),sml_eg)
         wbc(ix,iy,2)  = 0.0
         pbc(ix,iy,2)  = 0.0
         pbc2(ix,iy,2) = 0.0
      enddo
      enddo
c
      if(iocean .eq. 1) then
         call sufto(it)
         do iy=iys,iye
         do ix=1,nnx
            tau13m(ix,iy) = -au13m
            tau23m(ix,iy) = -au23m
         enddo
         enddo
         do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              taut3m(ix,iy,iscl) = aut3m(iscl)
           enddo
           enddo
         enddo
c
      else
c
         call suft(it)
         fac = -utau**2/(windm*sqrt(u1xy**2 + v1xy**2))
         do iy=iys,iye
         do ix=1,nnx
            tau13m(ix,iy)=fac*(windm*(u(ix,iy,iz)+ugal-u1xy)+
     +                     wind(ix,iy)*u1xy)
            tau23m(ix,iy)=fac*(windm*(v(ix,iy,iz)-v1xy)+
     +                     wind(ix,iy)*v1xy)
         enddo
         enddo
         do iscl=1,nscl
            dnom3=t10xy(iscl)*windm
            if(dnom3 .ne. 0.) then
               dnom_i = 1.0/dnom3
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl)=aut3m(iscl)*
     +                 (windm*(t(ix,iy,iscl,iz)-t1xy(iscl))+
     +                  wind(ix,iy)*(t1xy(iscl)-tsfcc(iscl)))*dnom_i
               enddo
               enddo
            else
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl) = aut3m(iscl)
               enddo
               enddo
            endif
         enddo
c
      endif
c
c -------- partial sums of surface fluxes and mean scalar
c
      sfc_flx(1) = 0.0
      sfc_flx(2) = 0.0
      do iy=iys,iye
      do ix=1,nnx
         sfc_flx(1) = sfc_flx(1) + tau13m(ix,iy)
         sfc_flx(2) = sfc_flx(2) + tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         sfc_flx(2+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            sfc_flx(2+iscl) = sfc_flx(2+iscl) + taut3m(ix,iy,iscl)
         enddo
         enddo
      enddo
c
      call mpi_sum_xy(sfc_flx,myid,iss,ise,(2+nscl))
      uwsfc = sfc_flx(1)*fnxy
      vwsfc = sfc_flx(2)*fnxy
      do iscl=1,nscl
         wtsfc(iscl) = sfc_flx(2+iscl)*fnxy
      enddo
c     write(nprt,2345) uwsfc, vwsfc, wtsfc(nscl), tsfcc(nscl)
 2345 format(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,
     +       ' wtsfc = ',e15.6,' tsfcc = ',e15.6)
c
      do iy=iys,iye
      do ix=1,nnx
         dudz     = 2.*(u(ix,iy,iz) + ugal)*dz_i
         dvdz     = 2.*v(ix,iy,iz)*dz_i
         ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
         vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
            tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
         enddo
         enddo
      enddo
c
c ------------ initialize u, v, w, t and derivatives at izm1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izm1)  = ubc(ix,iy,2)
         v(ix,iy,izm1)  = vbc(ix,iy,2)
         w(ix,iy,izm1)  = wbc(ix,iy,2)
         r3(ix,iy,izm1) =  0.0
         e(ix,iy,izm1)  = ebc(ix,iy,2)
         ux(ix,iy,izm1) = 0.0
         uy(ix,iy,izm1) = 0.0
         vx(ix,iy,izm1) = 0.0
         vy(ix,iy,izm1) = 0.0
         wx(ix,iy,izm1) = wbc(ix,iy,2)
         wy(ix,iy,izm1) = wbc(ix,iy,2)
      enddo
      enddo
c
c ------------- no need to call derivatives here since
c               wbc = 0, change for more general lower bc
c
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine pbltop(itop)
c
c ---------- get estimate of pbl top
c
c            method = 0, min of wt flux
c                        (good for buoyancy cases)
c            method = 1, uw flux less than critical value
c                        (good for ekman cases)
c            method = 2, running t average exceeds criterion
c                        (good for neutral cases with capping
c                         inversions)
c            method = 3, maximum gradient in temperature field
c                        (good for finding local zi see jas paper)
c                        with minimum search height (sr. setup)
c
c ------------ if method uses average statistics then only root
c              process need find zi
c
      use pars
      use fields
      use con_data
      use con_stats
      real trun(maxnz)
      include 'mpif.h'
      real gradloc(2,nnx,nny), gradmax(2,nnx,nny)
      external get_zi
c
      if(method .le. 2 .and. l_root) then
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = 1.0
      if (method .le. 0 .or. method .gt. 2) then
         itop=1
         wttot=wtle(1,1)+wtsb(1,1)
         wtmin=wttot*sgn
         do iz=2,nnz
            wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
            if (wttot.le.wtmin) then
               itop=iz
               wtmin=wttot
            endif
         enddo
         zi=z(itop)
      else if (method .eq. 1) then
         itop = 1
         crit = 0.05
         uwsf = utau*utau
         do iz=1,nnzm1
               uwtot = (uwle(iz) + uwsb(iz))**2 +
     $                 (vwle(iz) + vwsb(iz))**2
               uwtot = sqrt(uwtot)
               if(uwtot/uwsf .gt. crit) then
                  itop=iz
               endif
         enddo
         zi=z(itop)
      else if (method .eq. 2) then
         trun(1) = txym(1,1)
         do iz=2,nnz
             weight = z(iz-1)/z(iz)
             trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
         enddo
         itop = 1
         tcrit = 0.25
         if(iocean .eq. 1) tcrit = 0.1
         do iz=2,nnz
                if(txym(iz,1) .gt. (trun(iz) + tcrit)) then
                  itop = iz
                  go to 320
                endif
         enddo
  320    continue
         zi=z(itop)
      endif
      do iy=1,nny
      do ix=1,nnx
         gradmax(2,ix,iy) = zi
      enddo
      enddo
c
c ----------- use gradient method, every process computes
c
      elseif(method .eq. 3) then
c
c ---------------- get local zi from gradient in temperaure field
c
c     dz_i = dzu_i(izs+1)
c     do iy=1,nny
c     do ix=1,nnx
c        gradloc(1,ix,iy) = (t(ix,iy,1,izs+1) - t(ix,iy,1,izs))*dz_i
c        gradloc(2,ix,iy) = z(izs)
c     enddo
c     enddo
c
c ------- similar to zeroing the stat array in sr. mean_stat
c
      do iy=1,nny
      do ix=1,nnx
         gradloc(1,ix,iy) = 0.0
         gradloc(2,ix,iy) = z(iz_min)
      enddo
      enddo
c
c ------------- now all z in this process
c
      if(iz_min .le. ize) then
      do iz=max(izs,iz_min),ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            if(grad .gt. gradloc(1,ix,iy)) then
               gradloc(1,ix,iy) = grad
               gradloc(2,ix,iy) = z(iz)
            endif
         enddo
         enddo
      enddo
      endif
c
c     call mpi_reduce(gradloc,gradmax,2*nnx*nny,mpi_real8,ziloc,
c    +                i_root,mpi_comm_world,ierror)
c
c ----------- alternate version using already defined function in mpi
c             passes 2 real8 variables
c
      call mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,
     +                mpi_maxloc,i_root,mpi_comm_world,ierror)
c
c ------------ get average on root process
c
      if(l_root) then
         zi_avg = 0.0
         do iy=1,nny
         do ix=1,nnx
            zi_avg = zi_avg + gradmax(2,ix,iy)
         enddo
         enddo
         zi = zi_avg*fnxy
c        itop = nint(zi/dz)
      endif
c
      endif
c
c -------- send average zi everywhere
c
      call mpi_bcast(zi,1,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      if(iocean .ne. 1) then
         do iz=1,nnz
            if(zi .ge. z(iz) .and.
     +         zi .lt. z(iz+1)) itop = iz
         enddo
      else
         do iz=1,nnz
            if(zi .le. z(iz) .and.
     +         zi .gt. z(iz+1)) itop = iz
         enddo
      endif
c
c     if(l_root) write(6,7001) myid,zi,itop
 7001 format(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,
     +       ' itop = ',i3)
c
      return
      end
