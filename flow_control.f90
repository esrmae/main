module flow_control

  use shared_data

  implicit none

  real(mytype), dimension(:,:), allocatable ::  sig

Contains

!******************************************************************************

  subroutine initial_control
    implicit none 

!   This subroutine initialises wall osc parameters
    
    if (jet_type == 1) then
      call initial_wallosc
    else if (jet_type == 2) then
      call initial_bodyforce
    end if

    return
  end subroutine initial_control

!******************************************************************************

  subroutine initial_wallosc
    implicit none

    integer :: wavenum,step_periods
    real(mytype) :: angle

!   Finds the wall value of the amplitude
    if (wmax == 0.0e0) wmax=wmax_plus*retau/re

!   Finds the wall value of the frequency in time
    if (tperiod == 0.0e0 .and. tperiod_plus .ne. 0.0e0) then 
      omega_plus=2.0e0*pi/tperiod_plus
    else if (tperiod .ne. 0.0e0) then
      omeg = 2.0e0*pi/tperiod
    end if 
    if (omeg == 0.0e0) omeg=omega_plus*(retau**2/re)
    tperiod=2.0e0*pi/omeg

!   Finds the wall value of the frequency in space
    if (kx == 0.0e0) kx=kx_plus*retau
    wavenum=nint(kx*length/(2.0e0*pi))
    kx=real(wavenum,kind=mytype)*2.0e0*pi/length

!   Finds the wall value of the frequency in space for z
    if (kz == 0.0e0) kz=kz_plus*retau
    wavenum=nint(kz*width/(2.0e0*pi))
    kz=real(wavenum,kind=mytype)*2.0e0*pi/width

!   Finds the angle in radians
    angle = theta_s/180.0*pi
!   Finds the component of the angle in the w direction
    fact_w = cos(angle)
!   Finds the component of the angle in the v direction
    fact_v = sin(angle)
!   round off
    if( dabs(theta_s) < very_small ) then
      fact_w = 1.0e0; fact_v = 0.0e0
    end if
    if( dabs(theta_s-90.0e0) < very_small ) then
      fact_w = 0.0e0; fact_v = 1.0e0
    end if
!   Starts time from 0
    if (istart == 1) time=0.0e0
!   Sets dt_fixed to exact multiple of period
    if (ifdt_fixed == 1 .and. omeg.ne.0.0e0) then
      step_periods=int(tperiod/dt_fixed)
      dt_fixed=tperiod/real(step_periods,kind=mytype)
    end if

!   Print data to logfile
    if (nrank == 0) then
      write(21,171)
      write(21,*) '' 
      write(21,172) wmax,omeg,tperiod,kx,kz,fact_w,fact_v,dt_fixed,wavenum,step_periods
      write(21,*) '' 
      write(21,173)
    end if

  171    Format('************************* FLOW CONTROL PARAMETERS **************************')
  172    Format('  wmax=',E10.3,'  omeg=',E10.3,'  tperiod=',E10.3/, &
                '  kx=',E10.3, '  kz=',E10.3,'  fact_w=',E10.3,'  fact_v=',E10.3/, &
                '  dt_fixed=',E10.3,'  wavenum=',I5,'  steps/period=',I5)   
  173    Format('****************************************************************************')

    if (omeg .ne. 0.0e0 .and. kx .ne. 0.0e0 .and. lpstat_z .eq. 1) itravx=1
    if (omeg .ne. 0.0e0 .and. kz .ne. 0.0e0 .and. lpstat_x .eq. 1) itravz=1

    return
  end subroutine initial_wallosc

!******************************************************************************

  subroutine initial_bodyforce
    implicit none

    integer :: wavenum,step_periods
    real(mytype) :: angle

!   Finds the amplitude of force in outer scale
    if (St == 0.0e0) St=St_plus*(retau/re)**2

!   Finds the oscillating frequency in time
    if (tperiod == 0.0e0 .and. tperiod_plus .ne. 0.0e0) then 
      omega_plus=2.0e0*pi/tperiod_plus
    else if (tperiod .ne. 0.0e0) then
      omeg = 2.0e0*pi/tperiod
    end if 
    if (omeg == 0.0e0) omeg=omega_plus*(retau**2/re)
    tperiod=2.0e0*pi/omeg

!   Finds the wave number in streamwise direction
    if (kx == 0.0e0 .and. kx_plus .ne. 0.0e0) then  !Streamwise travelling wave
      kx=kx_plus*retau
      wavenum=nint(kx*length/(2.0e0*pi))
      kx=real(wavenum,kind=mytype)*2.0e0*pi/length
    end if
    if (kz == 0.0e0 .and. kz_plus .ne. 0.0e0) then !Spanwise travelling wave
      kz=kz_plus*retau
      wavenum=nint(kz*width/(2.0e0*pi))
      kz=real(wavenum,kind=mytype)*2.0e0*pi/width
    end if
!   Find the penetration depth in outer scale
    If(hjet == 0.0e0) Then
      hjet=hjet_plus/retau
    End If

!   Starts time from 0
    if (istart == 1) time=0.0e0
!   Sets dt_fixed to exact multiple of period
    if (ifdt_fixed == 1 .and. omeg.ne.0.0e0) then
      step_periods=int(tperiod/dt_fixed)
      if(lpstat_2d == 1 .And. nphase .ne. 1) dt_fixed=tperiod/real(step_periods,kind=mytype)
    end if
!   Initial force array
    fbody(:,:,:)=0.0e0
    call force_distribution

!   Print data to logfile
    if (nrank == 0) then
      write(21,181)
      write(21,*) '' 
      write(21,182) St,omeg,tperiod,kx,kz,hjet,dt_fixed,wavenum,step_periods
      write(21,*) '' 
      write(21,183)
    end if

  181    Format('************************* FLOW CONTROL PARAMETERS **************************')
  182    Format('  St=',E10.3,'  omeg=',E10.3,'  tperiod=',E10.3/, &
                '  kx=',E10.3,'  kz=',E10.3,'  hjet=',E10.3/, &
                '  dt_fixed=',E10.3,'  wavenum=',I5,'  step_periods=',I5)   
  183    Format('****************************************************************************')

    if (omeg .ne. 0.0e0 .and. kx .ne. 0.0e0 .and. lpstat_z .eq. 1) itravx=1
    if (omeg .ne. 0.0e0 .and. kz .ne. 0.0e0 .and. lpstat_x .eq. 1) itravz=1
!    if (kx .ne. 0.0e0 .and. kz .ne. 0.0e0) lpstat_2d=0

    return
  end subroutine initial_bodyforce

!******************************************************************************

  subroutine force_distribution 
    implicit none
    integer :: i,j,k,k1d
    real(mytype) :: dw
    real(mytype), dimension(:), allocatable :: fsta,fend
!	This subroutine calculates the force distribution
!   force is discretized in streamwise direction
    if(idx==1) then
      xfringe(:)=0.0e0
      allocate(fsta(1:ndx),fend(1:ndx))
      dw=length/real(ndx,mytype)
      do j=1,ndx
        fsta(j)=real(j-1,mytype)*dw
        fend(j)=real(j-1,mytype)*dw+dw/2.0e0
      end do
      do i=1,xsize(1); do j=1,ndx
        xfringe(i)=xfringe(i)+Sfun((x(i)-fsta(j))/rise)-Sfun((x(i)-fend(j)) &
                   /rise+1.0e0)
      end do; end do
      deallocate(fsta,fend)
    end if
!   force is discretized in spanwise direction    
    if(idz==1) then
      zfringe(:)=0.0e0
      allocate(fsta(1:ndz),fend(1:ndz))
      dw=width/real(ndz,mytype)
      do j=1,ndz
        fsta(j)=real(j-1,mytype)*dw
        fend(j)=real(j-1,mytype)*dw+dw/2.0e0
      end do
      do k=1,xsize(3); k1d=k+xstart(3)-1; do j=1,ndz
        zfringe(k)=zfringe(k)+Sfun((z(k1d)-fsta(j))/rise)-Sfun((z(k1d)-fend(j)) &
                   /rise+1.0e0)
      end do; end do
      deallocate(fsta,fend)
    end if

    return
  end subroutine force_distribution

!******************************************************************************

  real(mytype) function Sfun(x)
    real(mytype) :: x
    if(x<=0.0e0) then
      Sfun=0.0e0
    else if(x>=1) then
      Sfun=1.0e0
    else
      Sfun=1.0e0/(1.0e0+exp(1.0/(x-1)+1.0e0/x))
    end if
  end function Sfun

!******************************************************************************

  subroutine compute_flow_control 
    implicit none

!	This subroutine calculates the velcoity at wall/s

    allocate(sig(-1:xsize(1)+2,-1:xsize(3)+2))

!   jet_type=1 Wall Oscillation

	if (jet_type == 1) then
	  call wall_osc_calc
	  call calc_yinlet
    	else if (jet_type == 2) then
      	  call fbody_calc
!   update the force on the boundary
          Call update_halo(fbody(:,:,:),2)
          If(ymin(2).Eq.1) fbody(:,-1:0,:)=0.0e0
          If(ymax(2).Eq.1) fbody(:,xsize(2)+1:xsize(2)+2,:)=0.0e0
	end if
!
    deallocate (sig)

    return
  end subroutine compute_flow_control

!**********************************************************************************

  subroutine wall_osc_calc
    implicit none

    integer :: i,k

!   Waits for small velocity to start wall osc

    ampjet = 1.0e0
    !if (istart == 1) then
    !  if (time-time0 .lt. deltat) ampjet=1.0e0*(time-time0)/deltat
    !end if

!   Finds the value of sig, the force of a sine wave with amplitude 1, at time t
    if (amax .ge. 1.0e0) then
      do k=-1,xsize(3)+2; do i=1,xsize(1)
	sig(i,k)=ampjet*sin(kx*x(i)+kz*zc(k)-omeg*time)
      end do; end do
	
!   Can cut off the waves at a value amax 
    else if (amax .lt. 1.0) then
      do k=-1,xsize(3)+2; do i=1,xsize(1)
        sig(i,k)=ampjet*sin(kx*x(i)+kz*zc(k)-omeg*time)
        if (sig(i,k) .le. (-1.0*amax)) then       
          sig(i,k)=-1.0*amax
        end if
      end do; end do
    end if

    return
  end subroutine wall_osc_calc

!*****************************************************************************

  subroutine calc_yinlet
    implicit none

    integer :: i,k,i2
        
    do k=-1,xsize(3)+2; do i=1,xsize(1)
      yinlet1(i,k,2,1,1)=0.0e0
      yinlet1(i,k,2,2,1)=0.0e0
      yinlet1(i,k,2,1,2)=0.0e0
      yinlet1(i,k,2,2,2)=0.0e0
      yinlet1(i,k,3,1,1)=0.0e0
      yinlet1(i,k,3,2,1)=0.0e0
      yinlet1(i,k,3,1,2)=0.0e0
      yinlet1(i,k,3,2,2)=0.0e0
    end do; end do

!   Finds the component of the force in the v and w directions by multiplying the 
!   amplitude by the phase and the factor in each direction

    if (xmin(2) == 1) then 
      do k=-1,xsize(3)+2; do i=1,xsize(1)
!       Calculates the boundary conditions in the v direction for the lower wall
        yinlet1(i,k,2,1,1)= v1(i,1,k)
        yinlet1(i,k,2,1,2)= wmax*sig(i,k)*fact_v

!       Calculates the boundary conditions in the w direction for the lower wall
        yinlet1(i,k,3,1,1)= w1(i,0,k)
        yinlet1(i,k,3,1,2)= wmax*sig(i,k)*fact_w
      end do;end do
    end if

    If(nrank == 0) Write(*,*) 'vb=',wmax*sig(3,3)*fact_v,'wb=',wmax*sig(3,3)*fact_w
    if (xmax(2) == 1) then
      do k=-1,xsize(3)+2; do i=1,xsize(1)
!       Calculates the boundary conditions in the v direction for the upper wall
        yinlet1(i,k,2,2,1)= v1(i,xsize(2)+1,k)
        yinlet1(i,k,2,2,2)= wmax*sig(i,k)*fact_v

!       Calculates the boundary conditions in the w direction for the upper wall
        yinlet1(i,k,3,2,1)= w1(i,xsize(2)+1,k)
        yinlet1(i,k,3,2,2)= wmax*sig(i,k)*fact_w
      end do;end do
    end if

!   periodic in x
    do k=-1,xsize(3)+2; do i=-1,0
      yinlet1(i,k,:,:,:)=yinlet1(xsize(1)+i,k,:,:,:)
    end do; end do
    do k=-1,xsize(3)+2; do i=1,2
      yinlet1(xsize(1)+i,k,:,:,:)=yinlet1(i,k,:,:,:)
    end do; end do

    return
  end subroutine calc_yinlet

!***************************************************************************************

  subroutine fbody_calc
    implicit none
    integer :: i,j,k,j1d,k1d

    if(icub.eq.1) then
      do j=1,xsize(2)
        j1d=j+xstart(2)-1
        if(yc(j1d).le.1.0e0) then !LOWER WALL
          do k=1,xsize(3); k1d=k+xstart(3)-1; do i=1,xsize(1)
            if(iforce == 1) fbody(i,j,k)=St*exp(-yc(j1d)/hjet)*sin(kx*x(i)+kz*z(k1d)-omeg*time)
          end do; end do
        else !UPPER WALL
          do k=1,xsize(3); k1d=k+xstart(3)-1; do i=1,xsize(1)
            if(iforce == 1) fbody(i,j,k)=-St*exp(-(2.0e0-yc(j1d))/hjet)*sin(kx*x(i)+kz*z(k1d)-omeg*time)
          end do; end do
          if(kz.ne.0.0e0) fbody(:,j,:)=-fbody(:,j,:)
        end if
      end do      
    else
      do j=1,xsize(2)      
        j1d=j+xstart(2)-1
        if(yc(j1d).le.1.0e0) then
          do k=1,xsize(3); k1d=k+xstart(3)-1; do i=1,xsize(1)
          if(iforce == 1) then; fbody(i,j,k)=St*exp(-yc(j1d)/hjet)*sin(kx*x(i)+kz*z(k1d)-omeg*time)
!            fbody(i,j,k)=St*exp(-yc(j1d)/hjet)*sin(kz*z(k1d))*sin(-omeg*time)
          else if(iforce == 2) then
              fbody(i,j,k)=St*exp(-yc(j1d)/hjet)*max(sin(kz*z(k1d)),0.0e0)*sin(-omeg*time)
          end if; end do; end do
        end if
      end do
    end if
    if(idx==1) then
      do i=1,xsize(1) 
        fbody(i,:,:)=fbody(i,:,:)*xfringe(i)
      end do
    end if
    if(idz==1) then
      do k=1,xsize(3) 
        fbody(:,:,k)=fbody(:,:,k)*zfringe(k)
      end do
    end if
!    if(nrank==0.and.ii==1) then
!      open(unit=405,file=trim(folder)//'force_distr.dat')
!      do k=1,xsize(3)
!        j1d=k+xstart(3)-1
!        write(405,'(2F20.10)') z(j1d),zfringe(k)
!      end do
!      close(405)
!       do j=1,xsize(2) 
!         write(*,*) j,yc(j),fbody(10,j,10)
!       end do
!    end if 
  end subroutine fbody_calc
!******************************************************************************

 end module flow_control
