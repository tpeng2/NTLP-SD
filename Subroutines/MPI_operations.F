      subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
      real g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(j,i,k) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
      real g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(k,j,i) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         f(i,j,k) = ft(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_yz_var(temp_x,nvar,nny,iys,iye,izs,ize,ir)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real buf(nvar,iys:iye,izs:ize)
      real(kind=4), dimension(nvar,nny,izs:ize) :: temp_x
c
      num = nvar*(ize+1-izs)*(iye+1-iys)
      call mpi_recv(buf(1,iys,izs),num,mpi_real8,ir,1,
     +             mpi_comm_world,istatus,ierr)
      do k=izs,ize
      do j=iys,iye
      do ii=1,nvar
         temp_x(ii,j,k) = buf(ii,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do i=ixs,ixe
      do j=iys,iye
      do k=izs-1,ize+1
         f(i,j,k) = ft(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------- grab correct chunk of array to be sent and skip ghost points
c
      real f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         gt(j,i,k) = g(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent,
c               account for ghost points
c
      real g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do j=iys,iye
      do i=ixs,ixe
      do k=izs-1,ize+1
         gt(k,j,i) = g(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine mpi_gath_root(fs,fr,iz_s,iz_e,izs,ize,nz,myid,np,ns)
c
c ---------- gather results on root processors
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      integer iz_s(0:np-1), iz_e(0:np-1)
      real fs(izs:ize), fr(1:nz)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .gt. ns) then
        call mpi_send(fs(izs),ize+1-izs,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ns,np-1,ns
           ind = iz_s(l) + 1
           num = iz_e(l) + 1 - iz_s(l)
           call mpi_recv(fr(ind),num,mpi_real8,l,1,
     +       mpi_comm_world,istatus,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine mpi_send_root(fs,num,myid,np,ns)
c
c ---------- send root results to other processors above it
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real fs(num)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .ge. ns) then
        call mpi_recv(fs(1),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
        do l=irow_r+ns,np-1,ns
           call mpi_send(fs(1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine mpi_sum_xy(f,myid,iss,ise,nsend)
c
c --------- get horizontal x-y sum over a set of proccessors [iss:ise]
c           for vector f(i). f(i) is overwritten. skip if single processor
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real work(nsend,iss:ise), f(nsend)
c
      if(iss .eq. ise) go to 999
c
      do j=1,nsend
         work(j,myid) = f(j)
         f(j)         = 0.0
      enddo
      do i=iss,ise
         if(i .ne. myid) then
            call mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,
     +               work(1,i),nsend,mpi_real8,i,1,
     +           mpi_comm_world,istatus,ierr)
         endif
      enddo
      do i=iss,ise
      do j=1,nsend
         f(j) = f(j) + work(j,i)
      enddo
      enddo
c
  999 continue
c
      return
      end
      subroutine mpi_sum_z(f,i_root,myid,nsend,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend), f(nsend)
c
c -------- just root gets the result
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,
     +                  mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do i=1,nsend
               f(i) = recv_b(i)
            enddo
         endif
      else
c
c -------- everyone gets the result
c
         call mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,
     +                  mpi_comm_world,ierr)
         do i=1,nsend
            f(i) = recv_b(i)
         enddo
      endif
c
      return
      end
      subroutine mpi_sum_z_s(f,i_root,myid,nsend,nscl,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i,nscl)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend,nscl), f(nsend,nscl)
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum,i_root,mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do iscl=1,nscl
            do i=1,nsend
               f(i,iscl) = recv_b(i,iscl)
            enddo
            enddo
         endif
      else
         call mpi_allreduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum, mpi_comm_world,ierr)
         do iscl=1,nscl
         do i=1,nsend
            f(i,iscl) = recv_b(i,iscl)
         enddo
         enddo
      endif
c
      return
      end
      subroutine exchange
c
c ------------- exchange ghost points with mpi,
c               nb and nt are the destination and
c               source nodes. Allows for 1z per cpu
c
      use pars
      use fields
c     use fftwk
      include 'mpif.h'
c
      real fs(nnx,iys:iye,(4+nscl)),fr(nnx,iys:iye,(4+nscl))
      integer istatus(mpi_status_size)
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ account for endpoints
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)*(4+nscl)
      nrecv = nsend
c
c --------- send top of myid, receive bottom from myid - ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,ize)
         fs(ix,iy,2) = v(ix,iy,ize)
         fs(ix,iy,3) = w(ix,iy,ize)
         fs(ix,iy,4) = e(ix,iy,ize)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm1 = izs-1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1) = fr(ix,iy,1)
            v(ix,iy,izm1) = fr(ix,iy,2)
            w(ix,iy,izm1) = fr(ix,iy,3)
            e(ix,iy,izm1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c -------- send bottom of myid, receive bottom from myid + ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,izs)
         fs(ix,iy,2) = v(ix,iy,izs)
         fs(ix,iy,3) = w(ix,iy,izs)
         fs(ix,iy,4) = e(ix,iy,izs)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp1 = ize+1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izp1) = fr(ix,iy,1)
            v(ix,iy,izp1) = fr(ix,iy,2)
            w(ix,iy,izp1) = fr(ix,iy,3)
            e(ix,iy,izp1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c --------------- send extra scalar points 
c
      nsend = nnx*(iye + 1 - iys)*nscl
      nrecv = nsend
c
c -------------- send top of myid, receive bottom from myid - ncpu_s
c
      izm1 = ize-1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm2 = izs-2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
c -------------- send bottom of myid, receive bottom from myid + ncpu_s
c
      izp1 = izs+1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp2 = ize+2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array  f(nx,iys:iye,iz1:iz2) ---> g(ny,ixs:ixe,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
c
      do i=iss,ise
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iy_e(i) - iy_s(i) + 1)*ik
         if(i .eq. myid) then
            call send_xtoy(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
         else
            call send_xtoy(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
      enddo
c
      return
      end
      subroutine xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array  f(nx,iys:iye,izs-1:ize+1) 
c                     ---> g(0:nz+1,iys:iye,ixs:ixe)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+1)),
     +     gt(nz*(ixe+1-ixs)*(iye-iys+1))
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 1)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ----------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iz_e(i) - iz_s(i) + 1)*ij
         if(i .eq. myid) then
            call send_xtoz(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
         else
            call send_xtoz(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
      enddo
c
      return
      end
      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array g(ny,ixs:ixe,iz1:iz2) ---> f(nx,iys:iye,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
      do i=iss,ise
         nsend = (iy_e(i) - iy_s(i) + 1)*ik
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
         else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ytox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
      enddo
c
      return
      end
      subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array g(0:nz+1,iys:iye,ixs:ixe) 
c                    ---> f(nx,iys:iye,izs-1:ize+1)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+3)),
     +     gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))
c
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 3)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ------------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (iz_e(i) - iz_s(i) + 3)*ij
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
         else
            call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ztox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
      enddo
c
      return
      end
      subroutine bcast_pbc
c
c ---- send upper boundary conditions to other processors 
c      for fft solution of pressure
c
      use pars
      use fields
      include 'mpif.h'
      integer istatus(mpi_status_size),ierr
c
      if(numprocs .eq. 1) go to 999
c
      irow_r = mod(myid,ncpu_s)
      irow_t = is_s(numprocs-1) + irow_r
      num = nnx*(iye+1-iys)
c
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
