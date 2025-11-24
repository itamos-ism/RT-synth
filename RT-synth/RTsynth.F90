!Written by Thomas G. Bisbas
program RTsynth
use RTsynth_module
#ifdef HDF5
use m_readh5
#endif
#ifdef MULTILINES
integer::lrtot
#endif
#ifdef MULTIDIRECTION
integer::dir
character(len=2)::direction(1:6)
#endif

write(6,*) '----------------'
write(6,*) '-----RTsynth----'
write(6,*) '----------------'

#ifdef OPENMP
write(6,*) 'PARALLEL'
#else
write(6,*) 'SERIAL'
#endif
#ifdef HDF5
write(6,*) 'HDF5'
#endif
#ifdef VELOCITY
write(6,*) 'VELOCITY'
#endif
#ifdef DUST
write(6,*) 'DUST'
#endif
#ifdef BACKGROUND
write(6,*) 'BACKGROUND'
#endif
#ifdef CRATTENUATION
write(6,*) 'CRATTENUATION'
#endif
#ifdef HI21CM
write(6,*) 'HI21CM'
#endif
#ifdef MULTILINES
write(6,*) 'MULTILINES'
#endif
#ifdef MULTIDIRECTION
write(6,*) 'MULTIDIRECTION'
#endif

call readparams
call readcoolants
call initializations


#ifdef MULTIDIRECTION
direction(1) = '+z'
direction(2) = '+x'
direction(3) = '+y'
direction(4) = '-z'
direction(5) = '-x'
direction(6) = '-y'
do dir = 1,3
   los_direction = direction(dir)
#endif
   write(6,*) '--------------'
   write(6,*) 'LOS direction = ',los_direction
   write(6,*) '--------------'
#ifdef HDF5
   call readpdrfiles_h5
#else
   call readpdrfiles
#endif

   call columndensities

#ifdef HI21CM
   call hi21cm
#else
#ifdef MULTILINES
   do jr=1,5
     if (jr.eq.1) lrtot=7
     if (jr.eq.2) lrtot=1
     if (jr.eq.3) lrtot=2
     if (jr.eq.4) lrtot=2
     if (jr.eq.5) lrtot=7
     do lr=1,lrtot
#endif
       if (jr.eq.1.and.lr.eq.1) write(6,*) 'Calculating line CO (1-0)'
       if (jr.eq.1.and.lr.eq.2) write(6,*) 'Calculating line CO (2-1)'
       if (jr.eq.1.and.lr.eq.3) write(6,*) 'Calculating line CO (3-2)'
       if (jr.eq.1.and.lr.eq.4) write(6,*) 'Calculating line CO (4-3)'
       if (jr.eq.1.and.lr.eq.5) write(6,*) 'Calculating line CO (5-4)'
       if (jr.eq.1.and.lr.eq.6) write(6,*) 'Calculating line CO (6-5)'
       if (jr.eq.1.and.lr.eq.7) write(6,*) 'Calculating line CO (7-6)'
       if (jr.eq.2) write(6,*) 'Calculating line [CII] 158um'
       if (jr.eq.3.and.lr.eq.1) write(6,*) 'Calculating line [CI] (1-0)'
       if (jr.eq.3.and.lr.eq.2) write(6,*) 'Calculating line [CI[ (2-1)'
       if (jr.eq.4.and.lr.eq.1) write(6,*) 'Calculating line [OI] 63um'
       if (jr.eq.4.and.lr.eq.2) write(6,*) 'Calculating line [OI] 147um'
       if (jr.eq.5.and.lr.eq.1) write(6,*) 'Calculating line HCO+ (1-0)'
       if (jr.eq.5.and.lr.eq.2) write(6,*) 'Calculating line HCO+ (2-1)'
       if (jr.eq.5.and.lr.eq.3) write(6,*) 'Calculating line HCO+ (3-2)'
       if (jr.eq.5.and.lr.eq.4) write(6,*) 'Calculating line HCO+ (4-3)'
       if (jr.eq.5.and.lr.eq.5) write(6,*) 'Calculating line HCO+ (5-4)'
       if (jr.eq.5.and.lr.eq.6) write(6,*) 'Calculating line HCO+ (6-5)'
       if (jr.eq.5.and.lr.eq.7) write(6,*) 'Calculating line HCO+ (7-6)'

       do ci=1,ctot
         do cj=1,ctot
           do ck=1,ctot
             pdr(ci,cj,ck)%ttau=0.0D0                
             pdr(ci,cj,ck)%tTr=0.0D0
           enddo
         enddo
       enddo
       do j=1,coo
         coolant(j)%Tr=0.0D0
         coolant(j)%N=0.0D0
       enddo
#ifdef HDF5
       call readspopfile_h5_lr
#else
       call readspopfile
#endif
       call excitation
       call opticaldepth
       call radiativetransfer
       call writeoutputs
#ifdef MULTILINES
     enddo
   enddo
#endif
#endif
#ifdef MULTIDIRECTION
enddo
#endif

end program
