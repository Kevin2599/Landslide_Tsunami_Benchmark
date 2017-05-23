!**************************************************************************

      MODULE MovingBottom

!**************************************************************************
! ***  Routines to accomodate moving bottom
! ***  Solid slider testcases from Langford Sue
!**************************************************************************

!    logical, save :: OutputDetails
!      integer, save :: neMB, npMB, nsMB 
      real*8, save :: tandip,cosdip,sindip
      real*8, save :: Az, Bz, Cz, Xslide, xmaj2, ymaj2
      real*8, save :: xdif1,xdif2,xdif0
      real*8, save :: timein(60001),Txin(60001),deg2rad
      real*8, parameter :: SlopeAngle = 15.D0  !initial slope
      real*8, parameter :: zmin = -0.435D0     !float bottom elevation
      real*8, parameter :: xmaj = 0.25D0       !half of the slider length
      real*8, parameter :: ymaj = 0.026D0      !half of the slider height
      real*8, parameter :: pi = dacos(-1.D0)

      END MODULE

!**************************************************************************

      program driver
      
      implicit none
      
      integer :: np,j,jj
      real*8 :: time, xdist(1021),bottom(1021)
      
      np = 1021
      do j=1,np
        xdist(j)= -0.2 + 0.01D0*dble(j-1)
      enddo
      
      open(23,file='testbot83.dat')
      
      do jj=1,601  !test time loop
      
        time = 0.D0 + 0.01D0*dble(jj-1)
      
        call CalcBottomElevation(time,xdist,bottom,np)
      
!  *** write data in your format here. This is for Tecplot.
        write(23,*) 'ZONE'
        write(23,*) 'SOLUTIONTIME=' , time
        do j=1,np
          write(23,'(f10.2,1x,f10.5)') xdist(j),bottom(j)
        enddo
        
      enddo
      
      end

!**************************************************************************

      SUBROUTINE CalcBottomElevation(time,xdist,bottom,np)
      
      USE MovingBottom

      implicit none
      
! *** passed variables
      integer, intent(in) :: np      !number of points in x
      real*8 ,intent(in)  :: time    !time in seconds from start
      real*8 ,intent(in)  :: xdist(np)   !distance(m) from x=0 at intersection of slope and initial water surface
      real*8 ,intent(out) :: bottom(np)  !bottom elevation(m) measured from free surface

      integer :: jj,k1,jn
      integer, save :: jlast
      real*8 :: xcb,zbot,zbot2,xdif,x2new,arg,zbotleft,zbotright
      real*8 :: BzonAz,CzonAz,xleft,xright
      real*8 :: tx,ty,lx,ly,ox,oy
      real*8 :: stheta,ctheta,a2,b2
      real*8 :: zmin2,xdif0b,xdif1b,xdif2b
      logical :: start=.true.

! *** Case C1 UCant Flume (15 deg) with elliptical slider

!      pi = dacos(-1.D0)
      deg2rad = dabs(pi)/180.D0
      cosdip = cos(deg2rad*SlopeAngle)
      sindip = sin(deg2rad*SlopeAngle)
      tandip = tan(deg2rad*SlopeAngle)
      xmaj2 = 1./(xmaj*xmaj)
      ymaj2 = 1./(ymaj*ymaj)
      xdif0 = (xmaj+0.0625D0)*cosdip
      xdif1 = -0.8D0*xmaj*cosdip + 0.6D0*ymaj*sindip
      xdif2 = 0.8D0*xmaj*cosdip + 0.6D0*ymaj*sindip

      Az = xmaj2*sindip*sindip + ymaj2*cosdip*cosdip

      if(start) then
        open(51, file='TimeTx.dat')
        start = .false.
! *** read input location of trailing edge of slider
        do jj=1,60001
          READ(51,*) timein(jj), Txin(jj)
        enddo
        close(51)
        jlast = 1
      endif

! *** find Tx for time
      k1=jlast
      do jj=k1,60001
        jlast=jj
        Tx = Txin(jj)
        if(Tx.ge.time) exit
      enddo


        
      xslide = Tx + xmaj*cosdip  !slide cg
      xcb = xslide/cosdip

      do jn=1,np
          
        if (xdist(jn).le.1.297D0) then   ! at each time step, calculates Trailing Edge y-pos
          zmin2=-(xdist(jn)*tandip)   ! TE on planar slope
        elseif (xdist(jn).ge.1.807D0) then
          zmin2=-0.435D0                    ! TE on horizontal flume floor
        else                                ! TE on transition curve
          zmin2=(0.19D0*((xdist(jn)-1.297D0)**3))+ &
              (0.1024D0*((xdist(jn)-1.297D0)**2))-(.2728D0*(xdist(jn)-1.297D0))-0.34753D0 
        endif

        zbot = -xdist(jn)*tandip
        zbot2 = zmin
        
! *** slide down ramp if Lx<= 1.297 or Tx<=1.297-2.*xmaj*cosdip
        if(Tx.lt.1.297-2.*xmaj*cosdip) then
          xdif = xdist(jn) - xslide
          if(abs(xdif).lt.xdif0) then
            if(xdif.lt.xdif1) then  !linear end
              arg = (xdif+xdif0)/(xdif1+xdif0)
              zbotleft = -(xslide - xdif0)*tandip
              zbotright = -(xslide + xdif1)*tandip + 0.6D0*ymaj/cosdip
              zbot = zbotleft*(1.D0 - arg) + zbotright*arg
            elseif(xdif.lt.xdif2) then
!          if(abs(xdif).le.(xmaj*cosdip)) then
              xdif = xdist(jn)*cosdip - xcb
              x2new = xdist(jn)*sindip
              Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xdist(jn)*sindip*cosdip
              Cz = xmaj2*xdif*xdif + ymaj2*x2new*x2new - 1.D0
                BzonAz = 0.5*Bz/Az
                CzonAz = Cz/Az
                arg = BzonAz*BzonAz-CzonAz
              if(arg.ge.0.) then
                zbot = -BzonAz + sqrt(arg)
              else
                write(*,*) ' sqrt out of range'
                read(*,*)
              endif
            elseif(xdif.lt.xdif0) then !linear end
              arg = (xdif-xdif0)/(xdif2-xdif0)
              zbotright = -(xslide + xdif0)*tandip
              zbotleft =  -(xslide + xdif2)*tandip + 0.6D0*ymaj/cosdip
              zbot = zbotleft*arg + zbotright*(1.D0 - arg)
            endif
          endif
! *** slide across bottom if Tx>=1.807
        elseif(Tx.gt.1.807) then
          xdif = xdist(jn) - Tx - xmaj
          if(abs(xdif).lt.(xmaj+0.0625D0)) then
            if(xdif.lt.-0.8D0*xmaj) then  !linear end
              arg = (xdif+xmaj+0.0625D0)/(0.2D0*xmaj+0.0625D0)
              xleft = Tx - 0.0625D0
              if(xleft.lt.1.807) then
                zbotleft = (0.19D0*((xleft-1.297D0)**3))+(0.1024D0*((xleft-1.297D0)**2))&
                            -(.2728D0*(xleft-1.297D0))-0.34753D0
              else
                zbotleft = zmin
              endif
              zbotright = zmin + 0.6D0*ymaj
              zbot2 = zbotleft*(1.D0 - arg) + zbotright*arg
            elseif(xdif.lt.0.8D0*xmaj) then
              arg = 1.D0 - xmaj2*xdif*xdif
              if(arg.ge.0.D0) then
                zbot2 = zmin + ymaj*sqrt(arg)
              else
                write(*,*) ' sqrt out of range'
                read(*,*)
              endif
            elseif(xdif.lt.xmaj+0.0625D0) then !linear end
              arg = (xdif-xmaj-0.0625D0)/(-0.2D0*xmaj-0.0625D0)
              zbotright = zmin
              zbotleft =  zmin + 0.6D0*ymaj
              zbot2 = zbotleft*arg + zbotright*(1.D0 - arg)
            endif
          endif
        else !curved part
          Ty=(0.19D0*((Tx-1.297D0)**3))+(0.1024D0*((Tx-1.297D0)**2))&
              -(.2728D0*(Tx-1.297D0))-0.34753D0 
          Lx = Tx + sqrt(4.D0*xmaj**2 - (-0.435-Ty)**2)
          if(Lx.gt.1.807) then
            Ly = -0.435D0
            Ox = 0.5D0*(Tx+Lx)
            Oy = 0.5D0*(Ty+Ly)
          else
            ctheta = cosdip
            do k1=1,10  !iterate
              Lx = Tx + 2.D0*xmaj*ctheta
              Ly = (0.19D0*((Lx-1.297D0)**3))+(0.1024D0*((Lx-1.297D0)**2))&
                     -(.2728D0*(Lx-1.297D0))-0.34753D0 
              Ox = 0.5D0*(Tx+Lx)
              Oy = 0.5D0*(Ty+Ly)
              ctheta = (Ox - Tx)/xmaj
            enddo
          endif
          ctheta = (Ox - Tx)/xmaj
          stheta = (Ty - Oy)/xmaj
          xdif = xdist(jn) - Ox
          xdif0b = (xmaj+0.0625D0)*ctheta
          xdif1b = -0.8D0*xmaj*ctheta + 0.6D0*ymaj*stheta
          xdif2b = 0.8D0*xmaj*ctheta + 0.6D0*ymaj*stheta
!               Oy =  -Ox*tandip
          if(abs(xdif).lt.xdif0b) then
            if(xdif.lt.xdif1b) then  !linear end
              arg = (xdif+xdif0b)/(xdif1b+xdif0b)
              xleft = Tx - 0.0625D0*ctheta
              if(xleft.lt.1.297D0) then
                zbotleft = -xleft*tandip
              else
                zbotleft =  (0.19D0*((xleft-1.297D0)**3))+(0.1024D0*((xleft-1.297D0)**2))&
                            -(.2728D0*(xleft-1.297D0))-0.34753D0
              endif
              zbotright =   Oy + 0.8D0*xmaj*stheta + 0.6D0*ymaj*ctheta
              zbot2 = zbotleft*(1.D0 - arg) + zbotright*arg
            elseif(xdif.lt.xdif2b) then
              a2 = (stheta/xmaj)**2 + (ctheta/ymaj)**2
              b2 = -2.D0*xdif*(stheta/ctheta)/(xmaj**2)
              Cz = (xdif/ctheta/xmaj)**2 - 1.D0
              BzonAz = 0.5*B2/A2
              CzonAz = Cz/A2
              arg = B2**2 - 4.D0*a2*Cz
              if(arg.ge.0.) then
                zbot2 = Oy - xdif*stheta/ctheta - 0.5D0*b2/a2 + 0.5D0*sqrt(arg)/a2
              else
                write(*,*) ' sqrt out of range'
                read(*,*)
              endif
            elseif(xdif.lt.xdif0b) then !linear end
              arg = (xdif-xdif0b)/(xdif2b-xdif0b)
              xright = Ox + xdif0b
              xleft =  Ox + xdif2b
              if(xright.gt.1.807D0) then
                zbotleft =   Oy - 0.8D0*xmaj*stheta + 0.6D0*ymaj*ctheta
                zbotright =  zmin
              else
                zbotleft =   Oy - 0.8D0*xmaj*stheta + 0.6D0*ymaj*ctheta
                zbotright =  (0.19D0*((xright-1.297D0)**3))+(0.1024D0*((xright-1.297D0)**2))&
                            -(.2728D0*(xright-1.297D0))-0.34753D0                       
              endif
              zbot2 = zbotleft*arg + zbotright*(1.D0 - arg)
            endif
          endif
        endif
        bottom(jn) = max(zmin,zmin2,zbot,zbot2)
      enddo
                                
      END SUBROUTINE

!******************************************************************************
!******************************************************************************
