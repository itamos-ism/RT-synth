subroutine excitation
use RTsynth_module

!Calculate excitation temperature & black body emission
write(6,*) 'Calculating excitation temperature & black-body emission'
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ci,cj,ck)
#endif
do ci=1+ll,ctot-ll
  do cj=1+ll,ctot-ll
    do ck=1+ll,ctot-ll
       if (pdr(ci,cj,ck)%tpop(1).eq.0) then
           pdr(ci,cj,ck)%tTex = 0.
           pdr(ci,cj,ck)%tBnu = 0.
       else
           pdr(ci,cj,ck)%tTex = (HP*coolant(jr)%FREQUENCIES(lr+1,lr)/KB) / log(coolant(jr)%WEIGHTS(lr+1)*&
               pdr(ci,cj,ck)%tpop(1)/pdr(ci,cj,ck)%tpop(2)/coolant(jr)%WEIGHTS(lr))
#ifdef BACKGROUND
           pdr(ci,cj,ck)%tBnu = (2.*HP*coolant(jr)%FREQUENCIES(lr+1,lr)**3/C**2) / (dexp(HP*&
               coolant(jr)%FREQUENCIES(lr+1,lr)/KB/pdr(ci,cj,ck)%tTex)-1.0) - &
               (2.*HP*coolant(jr)%FREQUENCIES(lr+1,lr)**3/C**2) / (dexp(HP*&
               coolant(jr)%FREQUENCIES(lr+1,lr)/KB/Tcmb)-1.0)
#else
           pdr(ci,cj,ck)%tBnu = (2.*HP*coolant(jr)%FREQUENCIES(lr+1,lr)**3/C**2) / (dexp(HP*&
               coolant(jr)%FREQUENCIES(lr+1,lr)/KB/pdr(ci,cj,ck)%tTex)-1.0)
#endif
       endif
    enddo
  enddo
enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

return
end subroutine
