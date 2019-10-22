
      subroutine setup(it)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      it = iti
c
c ------------ turn on new sgs model at a particular step
c
      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif
c
      if(igrdr . eq. 3) then
         if(l_root) then
            write(6,6)iti,utau,tsfcc(1) ,qstar(1)
            write(6,510)
            write(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(6,it,1,nnz)
         endif
          if(l_debug) then
            write(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
            write(nprt,510)
            write(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(nprt,it,izs,ize)
          endif
      endif
c     if(ifilt.eq.1)call filter
      if(l_root) then
         write(6,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
      if(l_debug) then
         write(nprt,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
c
c -------------- boundary condition flags 
c
      ibcu = iradup
c     ibcu = 0
      ibcl = 0
c
c -------------------- wavenumbers, introduce a normalized
c                      set of wavenumbers to eliminate computation
c                      in derivatives , xderiv, yderiv
c
      do i=1,nnx
         xkn(i) = float(i-1)*pi2/xl
         if(i.gt.ncx)xkn(i) = -float(nnx-i+1)*pi2/xl
      enddo
      fn = 1.0/float(nnx)
      do i=1,nnx
         xk(i) = xkn(i)*fn
      enddo
      do i=1,nny
         ykn(i) = float(i-1)*pi2/yl
         if(i.gt.ncy)ykn(i) = -float(nny-i+1)*pi2/yl
      enddo
      fn = 1.0/float(nny)
      do i=1,nny
         yk(i) = ykn(i)*fn
      enddo
      ii = -1
      do i=1,ncx
         ii = ii + 2
         temp = xkn(i)**2
         do j=1,nny
            temp1       = temp + ykn(j)**2
            xks(ii,j)   = temp1
            xks(ii+1,j) = temp1
         enddo
      enddo
      xnn = abs(batag*dtdzf(1))
c
c ----------- choose correct sign so gravity waves
c             propagate out of the domain
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(ibcu.eq.1) then
         do iy=1,nny
         do ix=1,nnxp2
            if(xks(ix,iy) .le. 0.) then
              wavexy(ix,iy) = 0.0
            else
              wavexy(ix,iy) = sgn*sqrt(xnn/xks(ix,iy))
            endif
         enddo
         enddo
      endif
c
c -------------------- set length scale for SGS model
c
      if(iz_space .eq. 0) then
c
c ------------- uniform vertical spacing
c
      dx32 = dx*3./2.
      dy32 = dy*3./2.
      dsl  = (abs(dx32*dy32*dzw(1)))**(1./3.)
      dslg = dsl
      if(l_root)  write(6,2000) dsl
      if(l_debug) write(nprt,2000) dsl
c
c --------------------- create dsl array for easy indexing in comp1
c
      do iz=0,nnzp1
         dsl_z(iz) = dslg
      enddo
c
c ------------- variable vertical spacing
c
      else
c
c ----------- just estimate dsl for average spacing
c
         dx32 = dx*3./2.
         dy32 = dy*3./2.
c
         dsl_max = (abs(dx32*dy32*dzw(0)))**(1./3.)
         do iz=0,nnzp1
            dsl_z(iz) = (abs(dx32*dy32*dzw(iz)))**(1./3.)
            if(dsl_z(iz) .gt. dsl_max) dsl_max = dsl_z(iz)
         enddo
c        do iz=0,nnzp1
c           dsl_z(iz) = dsl_max
c        enddo
         dsl  = dsl_max
         dslg = dsl
      endif
c
      gridr = 1.0
      sml_eg = smal_e*gridr
c -------------------- set viscosity model parameters 
      if(ivis .ne. 1) then
        viscon = 0.0
        xksurf = 0.0
        nmatch = -1
        myid_newvis = 0
        do iz=1,nnz
           dfac(iz) = 1.0
        enddo
      endif
c ------------------- set stokes velocity for atmos/oceanic flow
      call stokesv
c
c --------- can add a time factor so as to skip into any part of
c           the specified geostrophic arrays. time factor in seconds
c
      t_factor = 7200.0
c
c ---------- for print out to get more digits
c
      t_ref = 300.0
c
c -------------------- specify cooling rate and initial
c                      temperature even for restarts
c
      c_rate   = 0.25/3600.0
      t_surf_i = 265.0
c
c -------------------- do not look for zi below zi_min
c
      zi_min = 30.0
      if(iocean .eq. 1) zi_min = -5.0
      iz_min = 1
      do iz=1,nnz-1
         if(zz(iz) .lt. zi_min .and.
     +      zz(iz+1) .ge. zi_min) iz_min = iz
      enddo
      if(l_root) then
         write(6,9000) zi_min, iz_min
      endif
c
 9998 continue
      return
c --------------------------- format statements
    6 format(///,' DATA FROM RESTART FILE AT STEP =',I5,
     +       ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510 format(' RESTART ***** CASE WITH : ******',/)
  520 format(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,
     +       ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x,
     +       '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4,
     +       ', UG = ',e12.4)
    1 format(10x,' NNX = ',i5,',  NNY = ',i5,
     + ',  NNZ = ',i5,/,10x,' SFC SMLT = ',i1,
     + ',  FILTER = ',i1,
     + ',  ITI = ',i6,',  ITMAX = ',i6,/,10x,
     + ' IUPWIND = ',i1,',  BUYNCY = ',i1,
     + ',  NO ALISNG = ',i1,',  ITCUT = ',i1,/,10x,
     + ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6,
     + ',  SUBSD = ',i1,/,
     + 10x,' BRCLICITY = ',i1,',  METHOD = ',i1,',  IOCEAN = ',i1,
     + ',  IVIS = ',i1)
 2000 format(10x,' DSL = ',e15.6)
 9000 format(' Search for zi above the height = ',e15.6,/,
     +       ' iz_min = ',i5)
      end
      subroutine set_range
c
c ---- build special x,y,z-ranges. dimensioned for 0:numprocs-1
c      indexed with myid
c
c      [ix_s:ix_e] x-range for computing y-derivatives nx-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jx_s:jx_e] x-range for computing 2d fft (nx+2)-pts/ncpu_s
c                  must be even in each x-interval for complex fft in y
c
c      [kx_s:kx_e] x-range for pressure solver transpose (nx+2)-pts/ncpu_z
c                  nx+2 fourier coefficients for xtoz and ztox transposes
c
c      [mx_s:mx_e] x-range split across z cpus as nx-pts/ncpu_z
c                  for use in surface layer routines 
c
c      [is_s:is_e] starting and ending processor id's for a
c                  particular z-level
c
c      [iy_s:iy_e] y-range for computing y-derivatives ny-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jy_s:jy_e] y-range for use in xtoz and ztox transposes 
c                  in pressure solution
c
c      [iz_s:iz_e] z-range for a particular vertical slab
c
c
      use pars
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      ii = -1
      do nn=0,ncpu_z-1
         call range(1,nnx+2,ncpu_z,nn,lx_s,lx_e)
         call range(1,nnx,ncpu_z,nn,nx_s,nx_e)
         call range(1,nny,ncpu_z,nn,ly_s,ly_e)
         call range(1,nnz,ncpu_z,nn,mz_s,mz_e)
         do mm=0,ncpu_s-1
            call range(1,nny,ncpu_s,mm,ny_s,ny_e)
            call range(1,nnx,ncpu_s,mm,nxy_s,nxy_e)
            call range(1,ncx,ncpu_s,mm,l2x_s,l2x_e)
            ii       = ii + 1
c
            ix_s(ii) = nxy_s
            ix_e(ii) = nxy_e
            jx_s(ii) = (l2x_s - 1)*2 + 1
            jx_e(ii) = l2x_e*2
            kx_s(ii) = lx_s
            kx_e(ii) = lx_e
            mx_s(ii) = nx_s
            mx_e(ii) = nx_e
c
            iy_s(ii) = ny_s
            iy_e(ii) = ny_e
            jy_s(ii) = ly_s
            jy_e(ii) = ly_e
c
            iz_s(ii) = mz_s
            iz_e(ii) = mz_e
c
            is_s(ii) = (ii/ncpu_s)*ncpu_s
            is_e(ii) = is_s(ii) + ncpu_s - 1
         enddo
      enddo
c
      iys = iy_s(myid)
      iye = iy_e(myid)
      jys = jy_s(myid)
      jye = jy_e(myid)
      ixs = ix_s(myid)
      ixe = ix_e(myid)
      jxs = jx_s(myid)
      jxe = jx_e(myid)
      kxs = kx_s(myid)
      kxe = kx_e(myid)
      mxs = mx_s(myid)
      mxe = mx_e(myid)
      izs = iz_s(myid)
      ize = iz_e(myid)
c
c ----------- get starting and  ending processor id's on each
c             vertical slab
c
      iss = is_s(myid)
      ise = is_e(myid)
c
c ------------ debug ranges 
c
      if(l_debug) then
         write(nprt,1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn),
     +                     jx_e(nn), kx_s(nn), kx_e(nn),
     +                     nn = 0,numprocs-1)
 1200    format(' myid =  ',i4,/,
     +       ' nn',5x,' ixs ',5x,' ixe ',5x,' jxs ',5x,' jxe '
     +       ,5x,' kxs ',5x,' kxe',/,(7i6))
c
         write(nprt,1213) myid, (nn, iy_s(nn), iy_e(nn),
     +                   jy_s(nn), jy_e(nn),
     +                   iz_s(nn), iz_e(nn), is_s(nn), is_e(nn),
     +                   nn=0,numprocs-1)
 1213    format(' myid = ',i3,/,
     +       ' nn ',3x,' iys ',5x,' iye ',5x,
     +       ' jys ',5x,' jye ',5x,
     +       ' izs ',5x,' ize',5x,' iss ',5x,' ise ',/,
     +       (9i6))
      endif
c
      return
      end
      subroutine set_sav(it,istart)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      data ionce /0/
      save ionce
c
      if(it .ne. istart) then
c
c ------------------- increment time if not first time through
c
         time=time+dt
      endif
c
      it=it+1
c
      dt    = dt_new
      dtg   = dt
      mnout = (mod(it,imean).eq.0).or. (it.eq.1)
      mtape = (mod(it,itape).eq.0)
      micut = (mod(it,itcut).eq.0)
      mviz  = (mod(it,i_viz).eq.0)
      if(ihst .lt. 0) then
         mhis = .false.
      else
         mhis = (mod(it,ihst).eq.0 .and. it .ge. it_his)
      endif
       mtrans = (mod(it+1,ihst).eq.0)

      if (i_viz .lt. 0) then
         msave_v = .false.
      else
         msave_v = (mod(it,i_viz).eq.0 .and. it .ge. it_viz)
      endif
c
c ---------- decide whether velocity fields are saved
c
      msave = .false.
      if(it .ge. itstr .and. mtape) then
         itn=itn+1
         msave = .true.
         call get_output_filenames
      endif
c
c ---------- decide whether viz fields are saved
c
      !msave_v = .false.
      !if(it .ge. itstr .and. mviz .and. i_viz .gt. 0) then
      !   msave_v = .true.
      !   if(ionce .eq. 0) then
      !      ionce = 1
      !      call open_viz
      !   endif
      !endif
      if((i_viz .gt. 0) .and. (it .ge. it_viz_nxt)) then
        call viz_output_filename(it)
        it_viz_nxt = it_viz_nxt + itape
      endif
c
c --------- decide whether history files are to be saved
c
      if((ihst .gt. 0) .and. (it .ge. it_his_nxt)) then
         call open_his(it)
         it_his_nxt = it_his_nxt + itape
      endif
c
      return
      end
