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
close(1);open(unit=1,file=adjustl(trim(outdir))//'/'//'RT_'//adjustl(trim(prefix))//'_'//&
      &adjustl(trim(los_direction))//adjustl(trim(suffix)),status='replace')
do ci=1+ll,ctot-ll
  do cj=1+ll,ctot-ll
#ifdef VELOCITY
    write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%W*coolant(jr)%Vunits(lr+1,lr)
#else
    write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%tTr*coolant(jr)%units(lr+1,lr)
#endif
  enddo
  write(1,*) ''
enddo

#ifdef VELOCITY
close(1);open(unit=1,file=adjustl(trim(outdir))//'/'//'RT_vel_'//adjustl(trim(prefix))//'_'&
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

close(1);open(unit=1,file=adjustl(trim(outdir))//'/'//'RT_tau_'//adjustl(trim(prefix))//'_'&
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


#ifdef CRATTENUATION
if (jr.eq.1.and.lr.eq.1) then
  close(1);open(unit=1,file=adjustl(trim(outdir))//'/'//'RT_'//adjustl(trim(prefix))//&
       '_'//adjustl(trim(los_direction))//'_Ncr.dat',status='replace')
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      write(1,'(3ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Ncr/pdr(ci,cj,1)%Ntot
    enddo
  write(1,*) ''
  enddo
endif
#endif

write(6,*) 'Finished!'
return
end subroutine
