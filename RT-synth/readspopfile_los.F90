subroutine readspopfile
use RTsynth_module

write(6,*) 'Reading level populations'
close(2);open(unit=2,file=adjustl(trim(filein))//trim(adjustl('.RTspop.fin')),status='old')
do i=1,nheader
  read(2,*) skipline
enddo

allocate(temp_coolant(1:coo))
do j=1,coo
  allocate(temp_coolant(j)%pop(1:coolant(j)%cnlev))
enddo
!continue reading RTspop.fin

select case(los_direction)
  case('+z')  ! Original case (z is the line-of-sight)
    do ci=1,ctot                     
      do cj=1,ctot                 
        do ck=1,ctot
          do j=1,coo
            read(2,'(100ES15.7)',advance='no')  temp_coolant(j)%pop(1:coolant(j)%cnlev)
          enddo
          read(2,*)
          pdr(ci,cj,ck)%tpop(1) = temp_coolant(jr)%pop(lr)
          pdr(ci,cj,ck)%tpop(2) = temp_coolant(jr)%pop(lr+1)
        enddo       
      enddo          
    enddo

  case('-z')  ! Reverse z direction
    do ci=1,ctot                     
      do cj=1,ctot                 
        do ck=ctot,1,-1
          do j=1,coo
            read(2,'(100ES15.7)',advance='no')  temp_coolant(j)%pop(1:coolant(j)%cnlev)
          enddo
          read(2,*)
          pdr(ci,cj,ck)%tpop(1) = temp_coolant(jr)%pop(lr)
          pdr(ci,cj,ck)%tpop(2) = temp_coolant(jr)%pop(lr+1)
        enddo       
      enddo          
    enddo

  case('+x')  ! x is the line-of-sight
    do ck=1,ctot                     
      do cj=1,ctot                 
        do ci=1,ctot
          do j=1,coo
            read(2,'(100ES15.7)',advance='no')  temp_coolant(j)%pop(1:coolant(j)%cnlev)
          enddo
          read(2,*)
          pdr(ci,cj,ck)%tpop(1) = temp_coolant(jr)%pop(lr)
          pdr(ci,cj,ck)%tpop(2) = temp_coolant(jr)%pop(lr+1)
        enddo       
      enddo          
    enddo

  case('-x')  ! Reverse x direction
    do ck=1,ctot                     
      do cj=1,ctot                 
        do ci=ctot,1,-1
          do j=1,coo
            read(2,'(100ES15.7)',advance='no')  temp_coolant(j)%pop(1:coolant(j)%cnlev)
          enddo
          read(2,*)
          pdr(ci,cj,ck)%tpop(1) = temp_coolant(jr)%pop(lr)
          pdr(ci,cj,ck)%tpop(2) = temp_coolant(jr)%pop(lr+1)
        enddo       
      enddo          
    enddo

  case('+y')  ! y is the line-of-sight
    do ci=1,ctot                     
      do ck=1,ctot                 
        do cj=1,ctot
          do j=1,coo
             read(2,'(100ES15.7)',advance='no')  temp_coolant(j)%pop(1:coolant(j)%cnlev)
          enddo
          read(2,*)
          pdr(ci,cj,ck)%tpop(1) = temp_coolant(jr)%pop(lr)
          pdr(ci,cj,ck)%tpop(2) = temp_coolant(jr)%pop(lr+1)
        enddo       
      enddo          
    enddo

  case('-y')  ! Reverse y direction
    do ci=1,ctot                     
      do ck=1,ctot                 
        do cj=ctot,1,-1
          do j=1,coo
            read(2,'(100ES15.7)',advance='no')  temp_coolant(j)%pop(1:coolant(j)%cnlev)
          enddo
          read(2,*)
          pdr(ci,cj,ck)%tpop(1) = temp_coolant(jr)%pop(lr)
          pdr(ci,cj,ck)%tpop(2) = temp_coolant(jr)%pop(lr+1)
        enddo       
      enddo          
    enddo
end select

close(2)
deallocate(temp_coolant)

return
end subroutine

