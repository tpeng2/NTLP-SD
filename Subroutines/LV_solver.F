        subroutine LV_solver(vnext,h,vec1,vec2,flag)
        use particles
        implicit none

        real, intent(in) :: vnext(3),h, vec1(2)
        real, intent(out) :: vec2(2)
        integer, intent(out) :: flag
        real :: error = 1E-8, fv1(2), fv2(2), v1(2), v_output(3), rel
        real :: diff, lambda,lup,ldown
        real :: C(2), newC(2), gradC(2), correct(2)
        real, dimension(1:2, 1:2) :: J,I,g,invg
        integer :: iterations, neg

        I = reshape((/1, 0, 0, 1/),shape(I))
        iterations = 0
        flag = 0
        v1 = vec1
        fv2 = (/1., 1./)

        lambda = 0.001
        lup = 2.0
        ldown = 2.0

        do while ((sqrt(dot_product(fv2, fv2)) > error).AND.
     +              (iterations<1000))

        iterations = iterations + 1
        call jacob_approx_2d(vnext, v1(1), v1(2), h,J)

        call ie_vrt_nd(vnext, v1(1), v1(2),v_output,fv1,h)
        g = matmul(transpose(J),J)+lambda*I
        gradC = matmul(transpose(J),fv1)
        C = 0.5*fv1*fv1

        call inverse_finder_2d(g, invg)
        correct = matmul(invg, gradC)
        if (sqrt(dot_product(correct,correct)) < 1E-12) then
                !flag = 1
                EXIT
        end if

        vec2 = v1 - correct
        call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,fv2,h)
        newC = 0.5*fv2*fv2

        if (sqrt(dot_product(newC,newC))<sqrt(dot_product(C,C))) then
                v1 = vec2
                lambda = lambda/ldown
        else
                lambda = lambda*lup
        end if

        end do

        if (iterations==1000) then
                flag = 1
                !print*,'not converged ', sqrt(dot_product(fv2, fv2))
        end if

        if (vec2(1) < 0 .OR. vec2(2) < 0) then
                flag = 1
                !print*, 'negative', sqrt(dot_product(fv2,fv2))
        end if


      end subroutine LV_solver
