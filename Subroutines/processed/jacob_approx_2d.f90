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
