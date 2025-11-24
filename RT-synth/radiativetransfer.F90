subroutine radiativetransfer
use RTsynth_module

!Solve radiative transfer equation
write(6,*) 'Solving radiative transfer equation'
#ifdef VELOCITY
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(vel,ci,cj,ck)&
!$OMP PRIVATE(dtau,Tr_incr1,Tr_incr2)
#endif
do vel = velmin, velmax
#endif
  do ci=1+ll,ctot-ll
    do cj=1+ll,ctot-ll
      do ck=1+ll,ctot-ll-1
#ifdef VELOCITY
         dtau=abs(pdr(ci,cj,ck)%ttau(vel)-pdr(ci,cj,ck+1)%ttau(vel))
#else
         dtau=abs(pdr(ci,cj,ck)%ttau-pdr(ci,cj,ck+1)%ttau)
#endif

         if (dtau.gt.1d10) then
#ifdef VELOCITY
#ifdef DUST
            pdr(ci,cj,1)%tTr(vel) = pdr(ci,cj,ck)%tSource(vel)
#else
            pdr(ci,cj,1)%tTr(vel) = pdr(ci,cj,ck)%tBnu
#endif

!else VELOCITY
#else

#ifdef DUST
            pdr(ci,cj,1)%tTr = pdr(ci,cj,ck)%tSource
#else
            pdr(ci,cj,1)%tTr = pdr(ci,cj,ck)%tBnu
#endif
#endif
         else if (dtau.gt.1d-6) then
#ifdef VELOCITY
#ifdef DUST
            Tr_incr1 = pdr(ci,cj,ck)%tSource(vel) * ((1.-dexp(-dtau))/dtau-dexp(-dtau))
            Tr_incr2 = pdr(ci,cj,ck+1)%tSource(vel) * (1.-(1.-dexp(-dtau))/dtau)
#else
            Tr_incr1 = pdr(ci,cj,ck)%tBnu * ((1.-dexp(-dtau))/dtau-dexp(-dtau))
            Tr_incr2 = pdr(ci,cj,ck+1)%tBnu * (1.-(1.-dexp(-dtau))/dtau)
#endif
            pdr(ci,cj,1)%tTr(vel)  = pdr(ci,cj,1)%tTr(vel) * dexp(-dtau) + Tr_incr1 + Tr_incr2
#else
#ifdef DUST
            Tr_incr1 = pdr(ci,cj,ck)%tSource * ((1.-dexp(-dtau))/dtau-dexp(-dtau))
            Tr_incr2 = pdr(ci,cj,ck+1)%tSource * (1.-(1.-dexp(-dtau))/dtau)
#else
            Tr_incr1 = pdr(ci,cj,ck)%tBnu * ((1.-dexp(-dtau))/dtau-dexp(-dtau))
            Tr_incr2 = pdr(ci,cj,ck+1)%tBnu * (1.-(1.-dexp(-dtau))/dtau)
#endif

            pdr(ci,cj,1)%tTr = pdr(ci,cj,1)%tTr * dexp(-dtau) + Tr_incr1 + Tr_incr2
#endif
         else
#ifdef VELOCITY
#ifdef DUST
            pdr(ci,cj,1)%tTr(vel) = pdr(ci,cj,1)%tTr(vel) * (1.-dtau) + (pdr(ci,cj,ck)%tSource(vel) + &
                    pdr(ci,cj,ck+1)%tSource(vel))*dtau/2.
#else            
            pdr(ci,cj,1)%tTr(vel) = pdr(ci,cj,1)%tTr(vel) * (1.-dtau) + (pdr(ci,cj,ck)%tBnu + &
                    pdr(ci,cj,ck+1)%tBnu)*dtau/2.
#endif

!else VELOCITY
#else

#ifdef DUST
            pdr(ci,cj,1)%tTr = pdr(ci,cj,1)%tTr * (1.-dtau) + (pdr(ci,cj,ck)%tSource + &
                    pdr(ci,cj,ck+1)%tSource)*dtau/2.
#else            
            pdr(ci,cj,1)%tTr = pdr(ci,cj,1)%tTr * (1.-dtau) + (pdr(ci,cj,ck)%tBnu + &
                    pdr(ci,cj,ck+1)%tBnu)*dtau/2.
#endif
#endif

         endif
      enddo
    enddo
  enddo
#ifdef VELOCITY
enddo !vel
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

!Velocity integrated emission
write(6,*) 'Calculating velocity integrated emission'
do ci=1+ll,ctot-ll
  do cj=1+ll,ctot-ll
    pdr(ci,cj,1)%W = 0
    do vel=velmin,velmax-1
      pdr(ci,cj,1)%W = pdr(ci,cj,1)%W + 0.5*(pdr(ci,cj,1)%tTr(vel) + pdr(ci,cj,1)%tTr(vel+1))*&
              coolant(jr)%frequencies(lr+1,lr)*abs(1./real(velstep))*1d5/C
    enddo
  enddo
enddo
#ifdef DUST
do ci=1+ll,ctot-ll
  do cj=1+ll,ctot-ll
      pdr(ci,cj,1)%W = pdr(ci,cj,1)%W -0.5*(pdr(ci,cj,1)%tTr(velmin)+pdr(ci,cj,1)%tTr(velmax))*(velmax-velmin)
  enddo
enddo
#endif
#endif

return
end subroutine
