subroutine hi21cm
use RTsynth_module
implicit none
double precision::A21=2.85d-15 !s-1
double precision::nu21=1420.405d6 !Hz
double precision::tau_factor,NHI,transmittance,T_HI_k,tau_k

write(6,*) 'Calculating HI optical depth and 21cm line intensity'
tau_factor = (3./32./PI)*(HK*C**2)*(A21/nu21)

do vel = velmin,velmax
  ds_freq = nu21 / (1. + real(vel)*1d5/C/real(velstep)) !Doppler-shifted frequency
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      !NHI=0.0D0
      transmittance = 1.0D0
      do ck=1+ll,ctot-ll-1
        step = abs(pdr(ci,cj,ck+1)%z-pdr(ci,cj,ck)%z)*PC
        sigma = (nu21/C) * SQRT(KB*pdr(ci,cj,ck)%Tgas/MH + vturb**2/2.)
        phi = 1. / sigma / sqrt(2.*PI)
        phi = phi * dexp(-((1+pdr(ci,cj,ck)%vz*1d5/C)*ds_freq - nu21)**2/2./sigma**2)
!        NHI = NHI + 0.5*(pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(32)+ &
!                pdr(ci,cj,ck+1)%rho*pdr(ci,cj,ck+1)%abun(32))*step
!        pdr(ci,cj,ck+1)%tau_HI(vel) = pdr(ci,cj,ck)%tau_HI(vel) + &
!                tau_factor*(NHI/pdr(ci,cj,ck)%Tgas)*phi
        NHI = pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(32)*step
        tau_k = tau_factor * (NHI/pdr(ci,cj,ck)%Tgas) * phi
        T_HI_k = (Tcmb*dexp(-tau_k) + pdr(ci,cj,ck)%Tgas * (1.0D0 - dexp(-tau_k)))*transmittance
        pdr(ci,cj,1)%T_HI(vel) = pdr(ci,cj,1)%T_HI(vel) + T_HI_k
        transmittance = transmittance * dexp(-tau_k)
      enddo
    enddo
  enddo
!  do ci=1+ll,ctot-ll
!    do cj=1+ll,ctot-ll
!      do ck=1+ll,ctot-ll-1
!        pdr(ci,cj,1)%T_HI(vel) = pdr(ci,cj,1)%T_HI(vel) + & 
!                pdr(ci,cj,ck+1)%Tgas*(1.-dexp(-pdr(ci,cj,ck+1)%tau_HI(vel)))
!      enddo
!    enddo
!  enddo
enddo

!write(6,*) 'Calculating velocity integrated 21cm intensity'
open(unit=22,file='HI21cm_W.dat',status='replace')
do ci=1+ll,ctot-ll
  do cj=1+ll,ctot-ll
    pdr(ci,cj,1)%W_HI = 0
    do vel=velmin,velmax-1
      pdr(ci,cj,1)%W_HI = pdr(ci,cj,1)%W_HI + 0.5*(pdr(ci,cj,1)%T_HI(vel) + &
              pdr(ci,cj,1)%T_HI(vel+1))*real(velstep)/real(velmax-velmin)
    enddo
    write(22,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%W_HI
  enddo
  write(22,*) ''
enddo

write(6,*) 'Writing output'
open(unit=21,file='HI21cm.dat',status='replace')
do vel=velmin,velmax
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      write(21,'(100ES15.7)') real(vel)/real(velstep),pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%T_HI(vel)
    enddo
    write(21,*) '' 
  enddo
  write(21,*) ''
enddo


write(6,*) 'Finished!'
return
end subroutine
