subroutine readparams
use RTsynth_module
#ifdef HDF5
use m_readh5
#endif
character(len=20)::cfile

close(1);open(unit=1,file=paramFile,status='old')
read(1,'(//)')
read(1,*) directory
read(1,*) outdir
read(1,*) prefix
filein = adjustl(trim(directory))//'/'//adjustl(trim(prefix))

! read header
#ifdef HDF5
  call readspopheader_h5
#else
  open(unit=2,file=adjustl(trim(filein))//trim(adjustl('.RTspop.fin')),status='old')
  nheader = 0
  read(2,*) network
  nheader = nheader + 1
  coo=0
  do 
    read(2,'(A)') cfile
    cfile=trim(adjustl(cfile))
    if (cfile=='ENDCOOLFILES') exit
    coo=coo+1
    coolfile(coo)=cfile
    write(6,*) coolfile(coo)
  enddo
  100 continue
  close(2)
  nheader = nheader + coo + 1 !1 for ENDCOOLFILES
#endif

if (network.eq.'REDUCED') then 
  nspec=33
  chemsuf = 'reduced.d'
endif
if (network.eq.'MEDIUM') then
  nspec=77
  chemsuf = 'medium.d'
endif
if (network.eq.'FULL') then 
  nspec=215
  chemsuf = 'full.d'
endif
allocate(species(1:nspec))
write(6,*) 'Chemical network = ',network
write(6,*) 'Coolants found = ',coo

read(1,*) velfile
filevel = adjustl(trim(directory))//"/"//adjustl(trim(velfile))
read(1,*) vturb
vturb = vturb*1d5 !km/s to cm/s
read(1,*) d2g !dust-to-gas normalized to 1e-2
read(1,*) redshift
Tcmb = 2.725*(1.+redshift)
read(1,'(//)')
read(1,*) los_direction
read(1,*) ctot
read(1,*) jr
read(1,*) lr
read(1,*) ll
#ifdef VELOCITY
read(1,*) velmin
read(1,*) velmax
read(1,*) velstep
#endif

return
end subroutine
