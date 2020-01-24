      subroutine vgrid_channel_fstrm(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      real s(0:2*nnz+1)
      integer :: zidx
      logical l_root, l_debug
c
      nnz = 2.0*nnz
      zi = 2.0*zi
      zl = 2.0*zl
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_stretop  = nnz !where strectching grids end, in free stream, set to top
c     n_stretop  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst  !center to first grid ratio
      z_fac2 = 1.0/float(n_stretop/2) !float strech increments
      z_fac  = 1.1 !stretching factor
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
      s(1) = z_frst
      do iz=2,n_stretop/2-1
         s(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      s(n_stretop/2) = z_cntr
      do iz=1,n_stretop/2 - 1
         s(n_stretop-iz) = zi - s(iz)
      enddo
      s(n_stretop) = zi
      s(0)   = 0.0
c
      if(l_root) write(6,5300) n_stretop
 5300 format(' n_stretop = ',i4)
c
c -------------- build grid from zi on up
!     For the channel, zi represents the channel centerline
!     Want the mesh to be a mirror image across this:
c
!      zidx = 1
!      do iz=n_stretop+1,nnz
!         z(iz) = zi + (zi - z(n_stretop-zidx))
!         zidx = zidx + 1
!      enddo
      s(nnz+1) = s(nnz) + (s(nnz) - s(nnz-1))

           nnz=nnz/2.0
           zi = zi/2.0
           zl = zl/2.0

           do iz=0,nnz
                   z(iz)=s(iz)
           enddo
            z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c
!          do iz=0,nnz+1
!          print*,'n = ',iz,'z = ',z(iz)
!          enddo
c
      return
      end subroutine vgrid_channel_fstrm
