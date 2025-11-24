subroutine readpdrfiles
use RTsynth_module
double precision::Av(1:12)

close(1);open(unit=1,file='chemfiles/species_'//trim(adjustl(chemsuf)),status='old')
do i=1,nspec
  read(1,*) index,species(i)
  do j=1,coo
    if (coolant(j)%cname.eq.species(i)) coolant(j)%cspec=i
  enddo
enddo

write(6,*) 'Reading input distribution'

pdrfile = trim(adjustl(filein))//".pdr.fin"
close(1);open(unit=1,file=pdrfile,status='old')

#ifdef CRATTENUATION
crfile = trim(adjustl(filein))//".cr.fin"
close(2);open(unit=2,file=crfile,status='old')
#endif

#ifdef VELOCITY
close(3);open(unit=3,file=filevel,status='old')
#endif

select case(los_direction)
  case('+z')  ! Original case (z is the line-of-sight)
    do ci=1,ctot                     
      do cj=1,ctot                 
        do ck=1,ctot
!          allocate(pdr(ci,cj,ck)%abun(1:nspec))
          read(1,*) id,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%Tgas,&
                  pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
#ifdef CRATTENUATION
          read(2,*) pdr(ci,cj,ck)%zetalocal
#endif
#ifdef VELOCITY
          read(3,*) pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vy,pdr(ci,cj,ck)%vz
#endif
        enddo       
      enddo          
    enddo
!    do ck=1,ctot
!      write(79,*) pdr(64,64,ck)%vz
!    enddo
!    stop


  case('-z')  ! Reverse z direction
    do ci=1,ctot                     
      do cj=1,ctot                 
        do ck=ctot,1,-1
!          allocate(pdr(ci,cj,ck)%abun(1:nspec))
          read(1,*) id,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%Tgas,&
                  pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
          pdr(ci,cj,ck)%z = -pdr(ci,cj,ck)%z
#ifdef CRATTENUATION
          read(2,*) pdr(ci,cj,ck)%zetalocal
#endif
#ifdef VELOCITY
          read(3,*) pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vy,pdr(ci,cj,ck)%vz
          pdr(ci,cj,ck)%vz = -pdr(ci,cj,ck)%vz
#endif
        enddo       
      enddo          
    enddo
!    do ck=1,ctot
!      write(78,*) pdr(64,64,ck)%vz
!    enddo
!    stop


  case('+x')  ! x is the line-of-sight
    do ck=1,ctot                     
      do cj=1,ctot                 
        do ci=1,ctot
!          allocate(pdr(ci,cj,ck)%abun(1:nspec))
          read(1,*) id,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%Tgas,&
                  pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
#ifdef CRATTENUATION
          read(2,*) pdr(ci,cj,ck)%zetalocal
#endif
#ifdef VELOCITY
          read(3,*) pdr(ci,cj,ck)%vz,pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vy
#endif
        enddo       
      enddo          
    enddo
!    do ck=1,ctot
!      write(81,*) pdr(64,64,ck)%z,pdr(64,64,ck)%vz
!    enddo
!    stop


  case('-x')  ! Reverse x direction
    do ck=1,ctot                     
      do cj=1,ctot                 
        do ci=ctot,1,-1
!          allocate(pdr(ci,cj,ck)%abun(1:nspec))
          read(1,*) id,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%Tgas,&
                  pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
!          pdr(ci,cj,ck)%z = -pdr(ci,cj,ck)%z
#ifdef CRATTENUATION
          read(2,*) pdr(ci,cj,ck)%zetalocal
#endif
#ifdef VELOCITY
          read(3,*) pdr(ci,cj,ck)%vz,pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vy
          pdr(ci,cj,ck)%vz = -pdr(ci,cj,ck)%vz
#endif
        enddo       
      enddo          
    enddo
!    do ck=1,ctot
!      write(80,*) pdr(64,64,ck)%z,pdr(64,64,ck)%vz
!    enddo
!    stop

  case('+y')  ! y is the line-of-sight
    do ci=1,ctot                     
      do ck=1,ctot                 
        do cj=1,ctot
!          allocate(pdr(ci,cj,ck)%abun(1:nspec))
          read(1,*) id,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%Tgas,&
                  pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
#ifdef CRATTENUATION
          read(2,*) pdr(ci,cj,ck)%zetalocal
#endif
#ifdef VELOCITY
          read(3,*) pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vz,pdr(ci,cj,ck)%vy
#endif
        enddo       
      enddo          
    enddo

  case('-y')  ! Reverse y direction
    do ci=1,ctot                     
      do ck=1,ctot                 
        do cj=ctot,1,-1
!          allocate(pdr(ci,cj,ck)%abun(1:nspec))
          read(1,*) id,pdr(ci,cj,ck)%x,pdr(ci,cj,ck)%z,pdr(ci,cj,ck)%y,pdr(ci,cj,ck)%Tgas,&
                  pdr(ci,cj,ck)%Tdust,etype,pdr(ci,cj,ck)%rho,uv,pdr(ci,cj,ck)%abun
          pdr(ci,cj,ck)%z = -pdr(ci,cj,ck)%z
#ifdef CRATTENUATION
          read(2,*) pdr(ci,cj,ck)%zetalocal
#endif
#ifdef VELOCITY
          read(3,*) pdr(ci,cj,ck)%vx,pdr(ci,cj,ck)%vz,pdr(ci,cj,ck)%vy
          pdr(ci,cj,ck)%vz = -pdr(ci,cj,ck)%vz
#endif
        enddo       
      enddo          
    enddo
end select

step = abs(pdr(1,1,2)%z-pdr(1,1,1)%z)*PC
!write(6,*) 'plength = ',step/PC

return
end subroutine

