subroutine opticaldepth
use RTsynth_module

!Calculate optical depth
write(6,*) 'Calculating optical depth'
#ifdef VELOCITY
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(vel,ds_freq)&
!$OMP PRIVATE(ci,cj,ck,sigma,frac,phi)&
!$OMP PRIVATE(alpha,alpha_dust,Bnu_dust)
#endif
do vel = velmin, velmax
  ds_freq = coolant(jr)%frequencies(lr+1,lr) / (1. + real(vel)*1d5/C/real(velstep)) !Doppler-shifted frequency
#endif
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      do ck=1+ll,ctot-ll-1
        sigma = (coolant(jr)%FREQUENCIES(lr+1,lr)/C) * SQRT(KB*pdr(ci,cj,ck)%Tgas/MH/coolant(jr)%molweight + vturb**2/2.)
        phi = 1. / sigma / sqrt(2.*PI)
        frac = (pdr(ci,cj,ck)%tpop(1)*coolant(jr)%WEIGHTS(lr+1)/coolant(jr)%WEIGHTS(lr)) - pdr(ci,cj,ck)%tpop(2)
#ifdef VELOCITY
        phi = phi * dexp(-((1+pdr(ci,cj,ck)%vz*1d5/C)*ds_freq - coolant(jr)%frequencies(lr+1,lr))**2/2./sigma**2)
        alpha = phi*(coolant(jr)%A_COEFFS(lr+1,lr)*C**2/8./pi/coolant(jr)%frequencies(lr+1,lr)**2)*frac*step
#ifdef DUST
        alpha_dust = 0.1d0*d2g*(ds_freq/1000d9)**2*pdr(ci,cj,ck)%rho*MH*1.402 !Arzoumanian+11
        Bnu_dust = (2.*HP*ds_freq**3/C**2)/(dexp(HP*ds_freq/KB/pdr(ci,cj,ck)%Tdust)-1)
        pdr(ci,cj,ck)%tSource(vel) = (alpha*pdr(ci,cj,ck)%tBnu + alpha_dust*Bnu_dust)/(alpha + alpha_dust)
        pdr(ci,cj,ck+1)%ttau(vel) = pdr(ci,cj,ck)%ttau(vel) + alpha + alpha_dust

#else 
        pdr(ci,cj,ck+1)%ttau(vel) = pdr(ci,cj,ck)%ttau(vel) + alpha!*step
#endif

!else VELOCITY
#else 
        alpha = phi*(coolant(jr)%A_COEFFS(lr+1,lr)*C**2/8./pi/coolant(jr)%frequencies(lr+1,lr)**2)*frac*step
#ifdef DUST 
        alpha_dust = 0.1d0*d2g*(coolant(jr)%FREQUENCIES(lr+1,lr)/1000d9)**2*pdr(ci,cj,ck)%rho*MH*1.402 !Arzoumanian+11
        Bnu_dust = (2.*HP*coolant(jr)%FREQUENCIES(lr+1,lr)**3/C**2)/&
                (dexp(HP*coolant(jr)%FREQUENCIES(lr+1,lr)/KB/pdr(ci,cj,ck)%Tdust)-1)
        pdr(ci,cj,ck)%tSource = (alpha*pdr(ci,cj,ck)%tBnu + alpha_dust*Bnu_dust)/(alpha + alpha_dust)
        pdr(ci,cj,ck+1)%ttau = pdr(ci,cj,ck)%ttau + (alpha + alpha_dust)*step
#else
        pdr(ci,cj,ck+1)%ttau = pdr(ci,cj,ck)%ttau + phi*(coolant(jr)%A_COEFFS(lr+1,lr)*C**2/8./pi/&
                coolant(jr)%frequencies(lr+1,lr)**2)*frac*step
#endif
#endif
      enddo
    enddo
  enddo
#ifdef VELOCITY
enddo !vel
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
#endif

return
end subroutine
