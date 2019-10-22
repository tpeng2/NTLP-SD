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
