subroutine readpdrfiles
use RTsynth_module

close(1);open(unit=1,file='chemfiles/species_'//trim(adjustl(chemsuf)),status='old')
do i=1,nspec
  read(1,*) index,species(i)
  do j=1,coo
    if (coolant(j)%cname.eq.species(i)) coolant(j)%cspec=i
  enddo
enddo

pdrfile = trim(adjustl(filein))//".pdr.fin"
close(1);open(unit=1,file=pdrfile,status='old')

write(6,*) 'Reading abundances'
do ci=1,ctot
  do cj=1,ctot
    do ck=1,ctot
      allocate(pdr(ci,cj,ck)%abun(1:nspec))
      read(1,*) id,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%Tgas,&
              pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
    enddo
  enddo
enddo



#ifdef CRATTENUATION
write(6,*) 'Reading CRIR'
crfile = trim(adjustl(filein))//".cr.fin"
close(1);open(unit=1,file=crfile,status='old')
do ci=1,ctot
  do cj=1,ctot
    do ck=1,ctot
      read(1,*) pdr(ci,cj,ck)%zetalocal
    enddo
  enddo
enddo
#endif

#ifdef VELOCITY
write(6,*) 'Reading velocities'
close(1);open(unit=1,file=filevel,status='old')
do ci=1,ctot
  do cj=1,ctot
    do ck=1,ctot
      read(1,*) pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vy,pdr(ci,cj,ck)%vz
    enddo
  enddo
enddo
#endif

return
end subroutine

