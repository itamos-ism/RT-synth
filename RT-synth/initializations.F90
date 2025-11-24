subroutine initializations
use RTsynth_module

allocate(Nabn(1:nspec))
write(6,*) 'Initializing'
#ifdef VELOCITY
velmin = velmin * velstep
velmax = velmax * velstep
#endif
allocate(pdr(1:ctot,1:ctot,1:ctot))
do ci=1,ctot
  do cj=1,ctot
    do ck=1,ctot
      allocate(pdr(ci,cj,ck)%abun(1:nspec))
      allocate(pdr(ci,cj,ck)%Nabn(1:nspec))
#ifdef VELOCITY
      allocate(pdr(ci,cj,ck)%ttau(velmin:velmax))
#ifdef DUST
      allocate(pdr(ci,cj,ck)%tSource(velmin:velmax))
#endif
      allocate(pdr(ci,cj,ck)%tTr(velmin:velmax))   
#ifdef HI21CM
      allocate(pdr(ci,cj,ck)%tau_HI(velmin:velmax))
      allocate(pdr(ci,cj,ck)%T_HI(velmin:velmax))
      pdr(ci,cj,ck)%tau_HI=0.0D0
      pdr(ci,cj,ck)%T_HI=0.0D0
#endif
#endif              
    enddo                                     
  enddo                                       
enddo                                         

do j=1,coo                                    
  allocate(coolant(j)%Tr(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
  allocate(coolant(j)%units(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
#ifdef VELOCITY
  allocate(coolant(j)%Vunits(1:coolant(j)%cnlev,1:coolant(j)%cnlev))
#endif
!  coolant(j)%Tr=0.0D0
!  coolant(j)%N=0.0D0
  do l=1,coolant(j)%cnlev-1
    coolant(j)%units(l+1,l) = C**2/2./KB/coolant(j)%frequencies(l+1,l)**2
#ifdef VELOCITY
    coolant(j)%Vunits(l+1,l) = C**3/2./KB/1D5/coolant(j)%frequencies(l+1,l)**3
#endif
  enddo
enddo

return
end subroutine
