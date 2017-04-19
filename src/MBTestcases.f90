!**************************************************************************

      MODULE MovingBottom

!**************************************************************************
! ***  Routines to accomodate moving bottom
! ***  Solid slider testcases from Langford Sue
!**************************************************************************

!    logical, save :: OutputDetails
!      integer, save :: neMB, npMB, nsMB 
      real*8, save :: xmaj,ymaj,zmin,tandip,cosdip,sindip
      real*8, save :: Az, Bz, Cz, Xslide, xmaj2, ymaj2
      real*8, save :: SlopeAngle,xdif1,xdif2,xdif0
      real*8 :: deg2rad

      END MODULE

!**************************************************************************

      SUBROUTINE InitSlide()
      
      USE SISL3DMOD
      USE MovingBottom

      implicit none
      
      integer :: nn,jj,k,k1,k6,jn,ncn2,mr,n1,n2,istat
      real*8 :: xcb,zc,zncn,zbot,xdif,hc0,x2new,arg
      real*8 :: BzonAz,CzonAz,dhonb,dLSonb,dLSonhc0,tmp,xleft,xright
      real*8 :: zbotleft,zbotright
      real*8 :: depc,xmin,radius,tx,ty,lx,ly,ox,oy,timein,diff,newdiff
      real*8 :: stheta,ctheta,s2theta,c2theta,det,a2,b2,ze,zmin2

! *** initialize landslide
! *** Case 80
! *** Case 81
! *** Case 82
! *** Case 83

      if(icase.ge.80.and.icase.le.83) then
! *** Case C1 UCant Flume (15 deg) with elliptical slider
          SlopeAngle = 15.D0
          zmin = -0.435D0
          xmaj = 0.25D0
          ymaj = 0.026D0

          pi = dacos(-1.D0)
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

          open(51, file='TimeTx.dat')
! *** read input
          READ(51,*) timein, Tx

          xslide = Tx + xmaj*cosdip
          xcb = xslide/cosdip

          do nn=1,ne
            ncn2 = ncn -1 + min0(1,nen(nn,ncn))
            zncn = 1./float(ncn2)
            zc = 0.
            DO k=1,ncn2
              jn=nen(nn,k)
              if(icase.gt.81) then
                if (xyz(jn,1).le.1.297D0) then   ! at each time step, calculates Trailing Edge y-pos
                  zmin2=-(xyz(jn,1)*tandip)   ! TE on planar slope
                elseif (xyz(jn,1).ge.1.807D0) then
                  zmin2=-0.435D0                    ! TE on horizontal flume floor
                else                                ! TE on transition curve
                  zmin2=(0.19D0*((xyz(jn,1)-1.297D0)**3))+ &
                      (0.1024D0*((xyz(jn,1)-1.297D0)**2))-(.2728D0*(xyz(jn,1)-1.297D0))-0.34753D0 
                endif
              else
                zmin2 = zmin
              endif
              zbot = -xyz(jn,1)*tandip
              xdif = xyz(jn,1) - xslide
              if(icase.eq.80.or.icase.eq.82) then
                if(abs(xdif).lt.(xmaj*cosdip)) then
                  xdif = xyz(jn,1)*cosdip - xcb
                  x2new = xyz(jn,1)*sindip
                  Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
                  Cz = xmaj2*xdif*xdif + ymaj2*x2new*x2new - 1.
                  BzonAz = 0.5*Bz/Az
                  CzonAz = Cz/Az
                  arg = BzonAz*BzonAz-CzonAz
                  if(arg.ge.0.) then
                    zbot = -BzonAz + sqrt(arg)
                  else
                    write(*,*) 'Serious problem in calculating ellipse- sqrt(<0)'
                    read(*,*)
                  endif
                endif                    
              else
                if(abs(xdif).lt.xdif0) then
                  if(xdif.lt.xdif1) then  !linear end
                    arg = (xdif+xdif0)/(xdif1+xdif0)
                    zbotleft = -(xslide - xdif0)*tandip
                    zbotright = -(xslide + xdif1)*tandip + 0.6D0*ymaj/cosdip
                    zbot = zbotleft*(1.D0 - arg) + zbotright*arg
                  elseif(xdif.lt.xdif2) then
                    xdif = xyz(jn,1)*cosdip - xcb
                    x2new = xyz(jn,1)*sindip
                    Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
                    Cz = xmaj2*xdif*xdif + ymaj2*x2new*x2new - 1.
                    BzonAz = 0.5*Bz/Az
                    CzonAz = Cz/Az
                    arg = BzonAz*BzonAz-CzonAz
                    if(arg.ge.0.) then
                      zbot = -BzonAz + sqrt(arg)
                    else
                      write(*,*) 'Serious problem in calculating ellipse- sqrt(<0)'
                      read(*,*)
                    endif
                  elseif(xdif.lt.xdif0) then !linear end
                    arg = (xdif-xdif0)/(xdif2-xdif0)
                    zbotright = -(xslide + xdif0)*tandip
                    zbotleft =  -(xslide + xdif2)*tandip + 0.6D0*ymaj/cosdip
                    zbot = zbotleft*arg + zbotright*(1.D0 - arg)
                  endif
                endif
              endif
              xyz(jn,3) = max(zmin,zmin2,zbot)
              zc = zc + xyz(jn,3)*zncn
            enddo
          enddo
      else
        write(*,*) 'WRONG slide option'
      endif
      
      END SUBROUTINE

!**************************************************************************

      SUBROUTINE Slide()
      
      USE SISL3DMOD
      USE MovingBottom

      implicit none
      
      integer :: nn,j,jj,jn,k,k1,k6,ncn2,mr,n1,n2
      integer :: istat,maxetmp
      real*8 :: zc,zncn,zbot,xdif,x2new,arg,zbotleft,zbotright
      real*8 :: xc,xcb,BzonAz,CzonAz,xleft,xright
!      real*8 :: etatmp(maxe2)
      real*8 :: depc,xmin,radius,tx,ty,lx,ly,ox,oy,timein
      real*8 :: stheta,ctheta,s2theta,c2theta,det,a2,b2,ze,fczc
      real*8 :: oldTime,oldTx
      real*8 :: diff, newdiff,xslide2,zbot2,zmin2,xdif0b,xdif1b,xdif2b


      if(icase.eq.80) then  !elliptical solid slider to mimic Langford's exp.

! *** read input
        do j=1,100
          READ(51,*) timein, Tx
        enddo
        
        xslide = Tx + xmaj*cosdip
        xcb = xslide/cosdip
        xslide2 = xslide+(zmin-xslide*tandip)*tan(0.5D0*deg2rad*SlopeAngle)

        do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              zncn = 1./float(ncn2)
              zc = 0.
              DO k=1,ncn2
                jn=nen(nn,k)
                zbot = -xyz(jn,1)*tandip
! *** slide down ramp
                xdif = xyz(jn,1) - xslide
                if(abs(xdif).lt.(xmaj*cosdip)) then
                  xdif = xyz(jn,1)*cosdip - xcb
                  x2new = xyz(jn,1)*sindip
                  Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
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
                endif
! *** slide across bottom
                xdif = xyz(jn,1) - xslide2
                zbot2 = zmin
                if(abs(xdif).lt.(xmaj)) then
!                  xdif = xyz(jn,1)*cosdip - xcb
!                  x2new = xyz(jn,1)*sindip
!                  Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
                  Cz = 1.D0 - xmaj2*xdif*xdif
!                  BzonAz = 0.5*Bz/Az
!                  CzonAz = Cz/Az
                  arg = Cz !BzonAz*BzonAz-CzonAz
                  if(arg.ge.0.) then
                    zbot2 = zmin + 0.25D0*ymaj*sqrt(arg)
                  else
                    write(*,*) ' sqrt out of range'
                    read(*,*)
                  endif
                endif
                xyz(jn,3) = max(zmin,zbot,zbot2)
                zc = zc + xyz(jn,3)*zncn
              enddo
        enddo

      elseif(icase.eq.81) then  !elliptical solid slider linear ends

! *** read input
        do j=1,100
          READ(51,*) timein, Tx
        enddo
        
        xslide = Tx + xmaj*cosdip
        xcb = xslide/cosdip
        xslide2 = xslide+(zmin-xslide*tandip)*tan(0.5D0*deg2rad*SlopeAngle)

        do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              zncn = 1./float(ncn2)
              zc = 0.
              DO k=1,ncn2
                jn=nen(nn,k)
                zbot = -xyz(jn,1)*tandip
! *** slide down ramp
                xdif = xyz(jn,1) - xslide
                if(abs(xdif).lt.xdif0) then
                  if(xdif.lt.xdif1) then  !linear end
                    arg = (xdif+xdif0)/(xdif1+xdif0)
                    zbotleft = -(xslide - xdif0)*tandip
                    zbotright = -(xslide + xdif1)*tandip + 0.6D0*ymaj/cosdip
                    zbot = zbotleft*(1.D0 - arg) + zbotright*arg
                  elseif(xdif.lt.xdif2) then
                    xdif = xyz(jn,1)*cosdip - xcb
                    x2new = xyz(jn,1)*sindip
                   Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
                    Cz = xmaj2*xdif*xdif + ymaj2*x2new*x2new - 1.
                    BzonAz = 0.5*Bz/Az
                    CzonAz = Cz/Az
                    arg = BzonAz*BzonAz-CzonAz
                    if(arg.ge.0.) then
                      zbot = -BzonAz + sqrt(arg)
                    else
                      write(*,*) 'Serious problem in calculating ellipse- sqrt(<0)'
                      read(*,*)
                    endif
                  elseif(xdif.lt.xdif0) then !linear end
                    arg = (xdif-xdif0)/(xdif2-xdif0)
                    zbotright = -(xslide + xdif0)*tandip
                    zbotleft =  -(xslide + xdif2)*tandip + 0.6D0*ymaj/cosdip
                    zbot = zbotleft*arg + zbotright*(1.D0 - arg)
                  endif
                endif
! *** slide across bottom
                xdif = xyz(jn,1) - xslide2
                zbot2 = zmin
                if(abs(xdif).lt.(xmaj+0.0625D0)) then
                  if(xdif.lt.-0.8D0*xmaj) then  !linear end
                    arg = (xdif+xmaj+0.0625D0)/(0.2D0*xmaj+0.0625D0)
                    zbotleft = zmin
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
                xyz(jn,3) = max(zmin,zbot,zbot2)
                zc = zc + xyz(jn,3)*zncn
              enddo
        enddo

      elseif(icase.eq.82) then  !elliptical solid slider

! *** read input
        do j=1,100
          READ(51,*) timein, Tx
        enddo
        
        xslide = Tx + xmaj*cosdip  !slide cg
        xcb = xslide/cosdip

        do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              zncn = 1./float(ncn2)
              zc = 0.
              DO k=1,ncn2
                jn=nen(nn,k)
                
                if (xyz(jn,1).le.1.297D0) then   ! at each time step, calculates Trailing Edge y-pos
                  zmin2=-(xyz(jn,1)*tandip)   ! TE on planar slope
                elseif (xyz(jn,1).ge.1.807D0) then
                  zmin2=-0.435D0                    ! TE on horizontal flume floor
                else                                ! TE on transition curve
                  zmin2=(0.19D0*((xyz(jn,1)-1.297D0)**3))+ &
                      (0.1024D0*((xyz(jn,1)-1.297D0)**2))-(.2728D0*(xyz(jn,1)-1.297D0))-0.34753D0 
                endif

                zbot = -xyz(jn,1)*tandip
! *** slide down ramp if Lx<= 1.297 or Tx<=1.297-2.*xmaj*cosdip
                if(Tx.le.1.297-2.*xmaj*cosdip) then
                  xdif = xyz(jn,1) - xslide
                  if(abs(xdif).lt.(xmaj*cosdip)) then
                    xdif = xyz(jn,1)*cosdip - xcb
                    x2new = xyz(jn,1)*sindip
                    Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
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
                  endif
! *** slide across bottom if Tx>=1.807
                elseif(Tx.ge.1.807) then
                  xdif = xyz(jn,1) - xslide
                  zbot2 = zmin
                  if(abs(xdif).lt.(xmaj)) then
                    arg = 1.D0 - xmaj2*xdif*xdif
                    if(arg.ge.0.) then
                      zbot = zmin + ymaj*sqrt(arg)
                    else
                      write(*,*) ' sqrt out of range'
                      read(*,*)
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
                    do k1=1,3  !iterate
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
                  xdif = xyz(jn,1) - Ox
                  if(abs(xdif).lt.(xmaj*ctheta)) then
!                    xdif = xyz(jn,1)*ctheta - Ox/ctheta  !xcb
!                    x2new = xyz(jn,1)*stheta
                    a2 = (stheta/xmaj)**2 + (ctheta/ymaj)**2
                    b2 = -2.D0*xdif*(stheta/ctheta)/(xmaj**2)
                    Cz = (xdif/ctheta/xmaj)**2 - 1.D0
                    BzonAz = 0.5*B2/A2
                    CzonAz = Cz/A2
                    arg = B2**2 - 4.D0*a2*Cz
                    if(arg.ge.0.) then
                      zbot = oy - xdif*stheta/ctheta - 0.5D0*b2/a2 + 0.5D0*sqrt(arg)/a2
                    else
                      write(*,*) ' sqrt out of range'
                      read(*,*)
                    endif
                  endif
                endif
                xyz(jn,3) = max(zmin,zmin2,zbot)
                zc = zc + xyz(jn,3)*zncn                
              enddo
        enddo
                                

      elseif(icase.eq.83) then  !elliptical solid slider linear ends

! *** read input
        do j=1,100
          READ(51,*) timein, Tx
        enddo
        
        xslide = Tx + xmaj*cosdip  !slide cg
        xcb = xslide/cosdip

        do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              zncn = 1./float(ncn2)
              zc = 0.
              DO k=1,ncn2
                jn=nen(nn,k)
                
                if (xyz(jn,1).le.1.297D0) then   ! at each time step, calculates Trailing Edge y-pos
                  zmin2=-(xyz(jn,1)*tandip)   ! TE on planar slope
                elseif (xyz(jn,1).ge.1.807D0) then
                  zmin2=-0.435D0                    ! TE on horizontal flume floor
                else                                ! TE on transition curve
                  zmin2=(0.19D0*((xyz(jn,1)-1.297D0)**3))+ &
                      (0.1024D0*((xyz(jn,1)-1.297D0)**2))-(.2728D0*(xyz(jn,1)-1.297D0))-0.34753D0 
                endif

                zbot = -xyz(jn,1)*tandip
                zbot2 = zmin
! *** slide down ramp if Lx<= 1.297 or Tx<=1.297-2.*xmaj*cosdip
                if(Tx.le.1.297-2.*xmaj*cosdip) then
                  xdif = xyz(jn,1) - xslide
                  if(abs(xdif).lt.xdif0) then
                    if(xdif.lt.xdif1) then  !linear end
                      arg = (xdif+xdif0)/(xdif1+xdif0)
                      zbotleft = -(xslide - xdif0)*tandip
                      zbotright = -(xslide + xdif1)*tandip + 0.6D0*ymaj/cosdip
                      zbot = zbotleft*(1.D0 - arg) + zbotright*arg
                    elseif(xdif.lt.xdif2) then
                      xdif = xyz(jn,1)*cosdip - xcb
                      x2new = xyz(jn,1)*sindip
                      Bz = -xmaj2*2.*xdif*sindip + ymaj2*2.*xyz(jn,1)*sindip*cosdip
                      Cz = xmaj2*xdif*xdif + ymaj2*x2new*x2new - 1.
                      BzonAz = 0.5*Bz/Az
                      CzonAz = Cz/Az
                      arg = BzonAz*BzonAz-CzonAz
                      if(arg.ge.0.) then
                        zbot = -BzonAz + sqrt(arg)
                      else
                        write(*,*) 'Serious problem in calculating ellipse- sqrt(<0)'
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
                elseif(Tx.ge.1.807) then
                  xdif = xyz(jn,1) - Tx - xmaj
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
                  xdif = xyz(jn,1) - Ox
                  xdif0b = (xmaj+0.0625D0)*ctheta
                  xdif1b = -0.8D0*xmaj*ctheta + 0.6D0*ymaj*stheta
                  xdif2b = 0.8D0*xmaj*ctheta + 0.6D0*ymaj*stheta
!                         Oy =  -Ox*tandip
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
                xyz(jn,3) = max(zmin,zmin2,zbot,zbot2)
                zc = zc + xyz(jn,3)*zncn                
              enddo
        enddo
      else
        write(*,*) 'WRONG slide module - unknown slide option'
      endif

      END SUBROUTINE

!******************************************************************************
!******************************************************************************
