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
