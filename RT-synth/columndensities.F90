subroutine columndensities
use RTsynth_module
double precision :: NH2,NHI,NCp,NC,NCO,NHCOp,Ntgas

write(6,*) 'Calculating column density maps'
close(1);open(unit=1,file=adjustl(trim(outdir))//'/'//'RT_'//adjustl(trim(prefix))//&
      &'_'//adjustl(trim(los_direction))//'_cds.dat',status='replace')
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
      Ntgas = Ntgas + pdr(ci,cj,ck)%rho*pdr(ci,cj,ck)%Tgas*step
    enddo
    write(1,'(100ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,NH2,NHI,NCp,NC,NCO,NHCOp,sum(pdr(ci,cj,:)%rho)*step, & 
      Ntgas / sum(pdr(ci,cj,:)%rho)/step
  enddo
write(1,*) ''
enddo

#ifdef CRATTENUATION
  close(1);open(unit=1,file=adjustl(trim(outdir))//'/'//'RT_'//adjustl(trim(prefix))//&
       '_'//adjustl(trim(los_direction))//'_Ncr.dat',status='replace')
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      write(1,'(3ES15.7)') pdr(ci,cj,1)%x,pdr(ci,cj,1)%y,pdr(ci,cj,1)%Ncr/pdr(ci,cj,1)%Ntot
    enddo
  write(1,*) ''
  enddo
#endif

#ifdef CDONLY
stop
#endif

return
end subroutine
