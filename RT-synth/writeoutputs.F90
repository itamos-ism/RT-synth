subroutine writeoutputs
use RTsynth_module
double precision :: NH2,NHI,NCp,NC,NCO,NHCOp,Ntgas

write(6,*) 'Writing outputs'
if (jr.eq.1.and.lr.eq.1) suffix="_CO10.dat"
if (jr.eq.1.and.lr.eq.2) suffix="_CO21.dat"
if (jr.eq.1.and.lr.eq.3) suffix="_CO32.dat"
if (jr.eq.1.and.lr.eq.4) suffix="_CO43.dat"
if (jr.eq.1.and.lr.eq.5) suffix="_CO54.dat"
if (jr.eq.1.and.lr.eq.6) suffix="_CO65.dat"
if (jr.eq.1.and.lr.eq.7) suffix="_CO76.dat"
if (jr.eq.2) suffix="_CII.dat"
if (jr.eq.3.and.lr.eq.1) suffix="_CI10.dat"
if (jr.eq.3.and.lr.eq.2) suffix="_CI21.dat"
if (jr.eq.4.and.lr.eq.1) suffix="_OI63.dat"
if (jr.eq.4.and.lr.eq.2) suffix="_OI146.dat"
if (jr.eq.5.and.lr.eq.1) suffix="_HCO+10.dat"
if (jr.eq.5.and.lr.eq.2) suffix="_HCO+21.dat"
if (jr.eq.5.and.lr.eq.3) suffix="_HCO+32.dat"
if (jr.eq.5.and.lr.eq.4) suffix="_HCO+43.dat"
if (jr.eq.5.and.lr.eq.5) suffix="_HCO+54.dat"
if (jr.eq.5.and.lr.eq.6) suffix="_HCO+65.dat"
if (jr.eq.5.and.lr.eq.7) suffix="_HCO+76.dat"
write(6,*) 'writing file: ','RT_'//adjustl(trim(prefix))//'_'//adjustl(trim(los_direction))//adjustl(trim(suffix))

!suf = '_'//trim(adjustl(coolant(jr)%cname))
!fix = 
!suffix = trim(adjustl(suf))//trim(adjustl(fix))

!close(1);open(unit=1,file='RT_'//adjustl(trim(prefix))//'.dat',status='replace')
close(1);open(unit=1,file='RT_'//adjustl(trim(prefix))//'_'//adjustl(trim(los_direction))//adjustl(trim(suffix)),status='replace')
do ci=1+ll,ctot-ll
  do cj=1+ll,ctot-ll
#ifdef VELOCITY
    write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%W*coolant(jr)%Vunits(lr+1,lr)
#else
    write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Ntot,pdr(ci,cj,1)%tTr*coolant(jr)%units(lr+1,lr)
#endif
  enddo
  write(1,*) ''
enddo

#ifdef VELOCITY
close(1);open(unit=1,file='RT_vel_'//adjustl(trim(prefix))//'_'&
        //adjustl(trim(los_direction))//adjustl(trim(suffix)),status='replace')
do vel=velmin,velmax-1
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      write(1,'(100ES15.7)') real(vel)/real(velstep),pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%tTr(vel)*coolant(jr)%units(lr+1,lr)
    enddo
    write(1,*) ''
  enddo
  write(1,*) ''
enddo

close(1);open(unit=1,file='RT_vel_tau_'//adjustl(trim(prefix))//'_'&
        //adjustl(trim(los_direction))//adjustl(trim(suffix)),status='replace')
do vel=velmin,velmax-1
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      write(1,'(100ES15.7)') real(vel)/real(velstep),pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,ctot-ll)%ttau(vel)
    enddo
    write(1,*) ''
  enddo
  write(1,*) ''
enddo

#endif


if (jr.eq.1.and.lr.eq.1) then
  close(1);open(unit=1,file='RT_'//adjustl(trim(prefix))//'_'//adjustl(trim(los_direction))//'_cds.dat',status='replace')
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      NH2=0; NHI=0; NCp=0; NC=0; NCO=0; NHCOp=0; Ntgas=0
      do ck=1+ll,ctot-ll
        NH2 = NH2 + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(31)*step
        NHI = NHI + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(32)*step
        NCp = NCp + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(11)*step
        NC = NC + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(25)*step
        NCO = NCO + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(28)*step
        NHCOp = NHCOp + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%abun(23)*step
        Ntgas = Ntgas + pdr(ci,ci,ck)%rho*pdr(ci,ci,ck)%Tgas*step
      enddo
      write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,NH2,NHI,NCp,NC,NCO,NHCOp,sum(pdr(ci,cj,:)%rho)*step, & 
        Ntgas / sum(pdr(ci,ci,:)%rho)/step
!      if (network.eq.'REDUCED') then
!         !write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Nabn(11),pdr(ci,cj,1)%Nabn(25),&
!         !     pdr(ci,cj,1)%Nabn(28),pdr(ci,cj,1)%Nabn(30),pdr(ci,cj,1)%Nabn(32),pdr(ci,cj,1)%Nabn(31),pdr(ci,cj,1)%Nabn(23)
!         write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,sum(pdr(ci,cj,:)%rho)*real(ctot)*step!*pdr(ci,cj,:)%abun(11))
!      else if (network.eq.'MEDIUM') then
!         write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Nabn(14),pdr(ci,cj,1)%Nabn(8),&
!              pdr(ci,cj,1)%Nabn(28),pdr(ci,cj,1)%Nabn(7)
!      else if (network.eq.'FULL') then
!         write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Nabn(172),pdr(ci,cj,1)%Nabn(210),&
!              pdr(ci,cj,1)%Nabn(211),pdr(ci,cj,1)%Nabn(213)
!      endif
    enddo
  write(1,*) ''
  enddo

#ifdef CRATTENUATION
  close(1);open(unit=1,file='RT_'//adjustl(trim(prefix))//'_'//adjustl(trim(los_direction))//'_Ncr.dat',status='replace')
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      write(1,'(3ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Ncr/pdr(ci,cj,1)%Ntot
    enddo
  write(1,*) ''
  enddo
#endif

endif

write(6,*) 'Finished!'
return
end subroutine
