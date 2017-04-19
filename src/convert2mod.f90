
  Program convert2mod

! *** read a grid file and convert the bottom elevation (z) to agree with the lab experiment
  
  implicit none
  
  integer :: i,j,ni,np,nngh,nadj(10,50000),ncode(50000)
  real*8 :: tp,x(50000),y(50000),z(50000)
  real*8,parameter :: pi = 3.141592653589793D0
  character(80) :: line, line2
 
! *** open an input and output file
  open(22,file='UCF-15mod2.ngh')  !input
  open(25,file='UCF-15mod3.ngh')  !corrected output
  
! *** change read format as required. This is for a *.ngh format file
  read(22,'(a)') line
  read(22,'(a)') line2
  read(22,*) np
  read(22,*) nngh
  
  tp = tan(pi*15.D0/180D0)
  write(*,*) ' tp=',tp

  do i=1,np
    read(22,*) ni,x(i),y(i),ncode(i),z(i),(nadj(j,i),j=1,nngh)
    if(x(i).lt.1.297D0) then
      z(i) = -x(i)*tp
    elseif(x(i).lt.1.807D0) then
      z(i)= 0.19D0*(x(i)-1.297D0)**3+0.1024D0*(x(i)-1.297D0)**2-0.2728D0*(x(i)-1.297D0)-0.3475D0 
    else
      z(i) = -0.435D0
    endif
  enddo
  
! *** change output format as required. This is for a *.ngh format file
  write(25,'(a)') line
  write(25,'(a)') line2
  write(25,*) np
  write(25,*) nngh
  do i=1,np
    write(25,'(I10,2(1x,F5.2),I4,1x,f8.5,10(i8))') i,x(i),y(i),ncode(i),z(i),(nadj(j,i),j=1,nngh)
  enddo
  
  end
  
  