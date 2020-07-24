subroutine initialize_parameters(matrices,params)

!This routine reads parameters in glide.in ...
!... calls prior_zc to estimate closure depths
!... calls isotherms to estimate perturbations of closure isotherm
!... applies this correction to the elevations of the ages


use definitions
use omp_lib

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
character :: line*1024,datafile*100
integer i,k,kk,nx,ny,j,nx_dum,ny_dum
integer ii,jj,skip,nx0,ny0
integer ib
double precision tmp1,tmp2,tmp3,tmp4,tmp5
double precision ran,colat,sys
double precision lo,la,spacin
double precision xtime,pi,summ,xl,yl,mean_elev
double precision xx_dum,yy_dum,dist,distmin
double precision xmin,xmax,ymin,ymax
double precision,dimension(:,:),allocatable::lat,lon,topob 
double precision,dimension(:,:),allocatable::lat_full,lon_full,topob_full
double precision,dimension(:),allocatable::topoa 

! "isys" is for the thermocrhonometric sys 1: aft; 2: zft; 3: ahe; 4: zhe

  call random_seed

!read in glide.in

  open (54,file='glide.in',status='unknown')!iostat=ios)
  open (55,status='scratch')
1 read (54,'(a1024)',end=2) line
  if (line(1:1).ne.'$'.and. line(1:1).ne.' ') write (55,'(a)') line
  goto 1
2 close (54)
  rewind (55)
  read (55,'(a5)') params%run
  read (55,*) params%topofile
  read (55,*) nx,ny
  read (55,*) datafile
  read (55,*) params%lon1,params%lon2,params%lat1,params%lat2
  read (55,*) params%edot_mean,params%sigma2
  read (55,*) params%xL,params%angle,params%aspect
  read (55,*) params%deltat
  read (55,*) params%t_total
  read (55,*) params%zl,params%Ts,params%Tb,params%kappa,params%hp
  read (55,*) params%iterM,params%xmu
  
  close(55)

!hp in micro Watts per m **3
params%hp=0.d0

!params%sigma2=(params%sigma2/params%edot_mean)**2
params%sigma2=(params%sigma2)**2
print*,params%sigma2
!spacin defines the resolution of the dummy points, 
!set to large value for test runs
spacin = .05
spacin=.1
spacin=1.5/params%xL
pi=atan(1.)*4.

!Setting up dummy points, at spacin times params%xL

!xl and yl are width and length of the analysis area 
xl = (params%lon2-params%lon1)*111.111*cos(((params%lat1+params%lat2)/2.)*pi/180.) 
yl = (params%lat2-params%lat1)*111.111

!We define the number of dummy points every n km, n is defined as spacin*xL
nx_dum = floor(xl/(spacin*params%xL))+1
ny_dum = floor(yl/(spacin*params%xL))+1
!print*,spacin*params%xL
!stop

!params%dummy=nx_dum*ny_dum

!dx_dum, xl, is therefore (lon2-lon1)/nx_dum
xl = (params%lon2-params%lon1)/float(nx_dum)
yl = (params%lat2-params%lat1)/float(ny_dum)

! reading in data,fission track
open(33,file=datafile,status="old")

!first go through, see which data fall in space (and time) range
k=0
do i=1,20000
 k=k+1
 read(33,*,end=500) lo,la,tmp1,tmp2,tmp3,tmp4,ib
 !print*,lo,la,tmp1,tmp2,tmp3,tmp4
 if (lo.gt.180.) lo = lo - 360.
 if (tmp2.gt.params%t_total.or. &
     lo.lt.params%lon1.or.lo.gt.params%lon2.or. &
     la.lt.params%lat1.or.la.gt.params%lat2) k=k-1
enddo
500 continue
params%n=k-1
rewind(33)

params%contr=37
!params%dummy=params%dummy+params%contr
params%dummy=params%contr

print*,"total ages",params%n,"control points=",params%dummy

!allocate arrays
allocate (matrices%ta(params%n),&
          matrices%a_error(params%n),&
          matrices%zc(params%n),&
          matrices%zcp(params%n),&
          matrices%x(params%n),&
          matrices%elev(params%n),&
          matrices%elev_true(params%n),&
          matrices%x_true(params%n),&
          matrices%y_true(params%n),&
          matrices%misfits(params%n),&
          matrices%syn_age(params%n),&
          matrices%isys(params%n),matrices%depth1(params%n),&
          matrices%iblock(params%n),&
          matrices%y(params%n),matrices%depth(params%n),matrices%zz(params%n))

k=0
do i=1,20000
       k=k+1
       read(33,*,end=300) lo,la,tmp1,tmp2,tmp3,sys,ib 
       if (lo.gt.180.) lo=lo-360.
       if (tmp2.gt.params%t_total.or. &
       lo.lt.params%lon1.or.lo.gt.params%lon2.or. &
       la.lt.params%lat1.or.la.gt.params%lat2) then     
       k=k-1
     else
     call random_number(ran)
       matrices%x(k)=lo
       matrices%y(k)=la
       matrices%elev(k)=tmp1
       matrices%ta(k)=tmp2 !+((ran-0.5)*2.)
!if (matrices%ta(k).lt.10.) matrices%ta(k)=10.
       call random_number(ran)
       matrices%a_error(k)=tmp3
       matrices%isys(k)=sys
       matrices%iblock(k)=ib
       if (k.eq.params%n) exit
     endif
enddo
300 continue
close(33)
k=k-1

!elev is modified to account for perturbation of the isotherm
matrices%elev_true=matrices%elev

!print*,'creating dummy points'
!k=0
!do i=1,nx_dum
!  do j=1,ny_dum
!    k=k+1
!     xx_dum=params%lon1 + float((i-1))*xl
!     yy_dum=params%lat1 + float((j-1))*yl
!     kk=0
!     do jj=1,params%n
!      dist=((xx_dum-matrices%x(jj))**2+(yy_dum-matrices%y(jj))**2)
!      if (dist.lt.(xl*10.*xl*10.)) then
!        kk=1
!        goto 982
!      endif
!     enddo
!982    if (kk.eq.0) k=k-1
!enddo
!enddo
!print*,"reduced number of dummy points=",k,"from=",nx_dum*ny_dum
open(988,file='data/grid_dummy.txt',status='old')
! check number of grid points first
do j=1,1000000
read(988,*,end=514) lo,la,ib
enddo
514 continue
rewind(988)
!read in the grid
k=0
do i=1,j-1
 k=k+1
 read(988,*,end=501) lo,la,ib
 if (lo.gt.180.) lo = lo - 360.
 if (lo.lt.params%lon1.or.lo.gt.params%lon2.or. &
     la.lt.params%lat1.or.la.gt.params%lat2) then 
   !print*,lo,la,ib
   k=k-1
 endif
enddo
501 continue
rewind(988)
params%dummy=k-1+params%contr
print*,'number of dummy point=',params%dummy
allocate(matrices%x_dum(params%dummy),matrices%y_dum(params%dummy))
allocate( matrices%x_dum_true(params%dummy),&
          matrices%y_dum_true(params%dummy),&
          matrices%idum_block(params%dummy))
k=0
do i=1,10000000
  k=k+1
  read(988,*,end=502) lo,la,ib 
 if (lo.gt.180.) lo = lo - 360.
 if (lo.lt.params%lon1.or.lo.gt.params%lon2.or. &
     la.lt.params%lat1.or.la.gt.params%lat2) then
   k=k-1
 else
  matrices%x_dum(k)=lo
  matrices%y_dum(k)=la 
  matrices%idum_block(k)=ib
 endif
  if (k.eq.params%dummy-params%contr) exit
enddo
502 continue
close(988)
!read in control points
!k=k-1
open(999,file="data/control_id.xy",status="old")
do i=1,params%contr
  k=k+1
  read(999,*) matrices%x_dum(k),matrices%y_dum(k),matrices%idum_block(k)
enddo
close(999)
print*,minval(matrices%y_dum),maxval(matrices%y_dum)
print*,minval(matrices%x_dum),maxval(matrices%x_dum)

! positions are later transformed, keep lat lon values
matrices%x_dum_true=matrices%x_dum
matrices%y_dum_true=matrices%y_dum
matrices%x_true=matrices%x
matrices%y_true=matrices%y

!and convert them into local coordinate system
colat=cos(((params%lat1+params%lat2)/2.)*pi/180.)

!convert points to local grid
do i=1,params%dummy
  matrices%x_dum(i) = (matrices%x_dum(i)-params%lon1)*111.11*colat
  matrices%y_dum(i) = (matrices%y_dum(i)-params%lat1)*111.11 
enddo

!give data errors if needed
where (matrices%a_error.eq.0) matrices%a_error=0.1!*matrices%ta

!open topography, should be a bit bigger than data extents
open(49,file=params%topofile,status="old") 

nx0=nx
ny0=ny

skip=1

nx=(nx0-1)/skip+1
ny=(ny0-1)/skip+1
nx=nx0
ny=ny0

allocate(topoa(nx*ny),lon(nx,ny),lat(nx,ny),topob(nx,ny))

allocate(lon_full(nx0,ny0),lat_full(nx0,ny0),topob_full(nx0,ny0))

do j=ny0,1,-1
  do i=1,nx0
    read(49,*) lon_full(i,j),lat_full(i,j),topob_full(i,j)
  enddo
enddo
close(49)

lon=lon_full
lat=lat_full
topob=topob_full

xmin=dble(minval(lon))
xmax=dble(maxval(lon))
ymin=dble(minval(lat))
ymax=dble(maxval(lat))
print*,xmin,xmax,ymin,ymax


! It first allocates the table size
params%m_max=floor((params%t_total)/params%deltat)+1

allocate(matrices%tsteps(100),matrices%tsteps_sum(100))

matrices%tsteps=params%deltat

!define time steps can be !variable in time
do i=1,100
!  matrices%tsteps(i) = 2.5*(1-exp(-(float(i)/10.)))+1.
  ! matrices%tsteps(i) = 2.*(1-exp(-(float(i)/10.)))+1.
  if (i.gt.1) then
    matrices%tsteps_sum(i)=matrices%tsteps_sum(i-1) + matrices%tsteps(i)
  else
    matrices%tsteps_sum(i)=matrices%tsteps(i)
  endif
enddo

!calculate total number of time steps required to describe M_max
summ=0.
k=1
params%m_max=0
do while (summ.lt.params%t_total)
  summ = summ + matrices%tsteps(k)
  k=k+1
enddo
k=k-1
k=k+1
params%m_max = k

allocate(matrices%edot_pr(params%n*params%m_max))
allocate(matrices%edot_pr2(params%n*params%m_max))
allocate(matrices%edot_pr_dum(params%dummy*params%m_max))
allocate(matrices%nsystems(8))


!build prior model, can be complex if required
xtime=params%t_total
do j=1,params%m_max
  xtime=xtime-params%deltat
  do i=1,params%n
    matrices%edot_pr(j+(i-1)*params%m_max) = params%edot_mean
  enddo
  do i=1,params%dummy
    matrices%edot_pr_dum(j+(i-1)*params%m_max) = params%edot_mean
  enddo
enddo

!number of samples with different systems
matrices%nsystems=0
do i=1,params%n
if (sngl(matrices%isys(i)).eq.1) matrices%nsystems(1)=matrices%nsystems(1)+1
if (sngl(matrices%isys(i)).eq.2) matrices%nsystems(2)=matrices%nsystems(2)+1
if (sngl(matrices%isys(i)).eq.3) matrices%nsystems(3)=matrices%nsystems(3)+1
if (sngl(matrices%isys(i)).eq.4) matrices%nsystems(4)=matrices%nsystems(4)+1
if (matrices%isys(i).eq.5) matrices%nsystems(5)=matrices%nsystems(5)+1
if (matrices%isys(i).eq.6) matrices%nsystems(6)=matrices%nsystems(6)+1
if (matrices%isys(i).eq.7) matrices%nsystems(7)=matrices%nsystems(7)+1
if (matrices%isys(i).lt.0.) matrices%nsystems(8)=matrices%nsystems(8)+1
enddo

allocate(matrices%ages(7,5))

!calculate mean elev and temperature
mean_elev = sum(topob)/float(nx*ny)
params%Ts = params%Ts - mean_elev*0.006
params%zl = params%zl + (mean_elev/1000.)

print*,"Temperature of thermal model at z=0",params%Ts,"at an elevation of",mean_elev

call prior_zc(matrices,params)

!matrices%zz contains the noisy closure depths
matrices%zz=matrices%zc
!print*,matrices%zz
!print*,'oops'
!pause

!zc is used for the relief of the closure surface
matrices%zc=0.d0

!remove mean from samples and topography
topob=topob-mean_elev
matrices%elev=matrices%elev-mean_elev

call isotherms(matrices,params,nx,ny,topob,lon,lat)

deallocate(topoa,lon,lat,topob)

!convert to local coordinate system  
do i=1,params%n
 matrices%x(i) = (matrices%x_true(i)-params%lon1)*111.11*colat
 matrices%y(i) = (matrices%y_true(i)-params%lat1)*111.11
enddo

!open(81,file=params%run//"/stuff/aer_pr.aze",status="unknown")

!do i=1,params%n
!  write(81,'(8f12.4,I4)') matrices%x(i),matrices%y(i),matrices%ta(i),matrices%elev(i),matrices%elev_true(i), &
!                          matrices%zc(i),matrices%zz(i),matrices%a_error(i),matrices%isys(i)
!enddo

!close(81)

! remove the topography on closure isotherm from the elevations
do i=1,params%n
   if (matrices%isys(i).gt.0) then
   matrices%elev(i)=(matrices%elev(i)-matrices%zc(i))/1000.
   else
   matrices%elev(i)=0.d0
   matrices%zz(i)=-matrices%isys(i)
   endif
enddo

!!printing the closure depths
print*,"AFT=",matrices%ages(1,2),"km",matrices%ages(1,1),"degC"
print*,"ZFT=",matrices%ages(2,2),"km",matrices%ages(2,1),"degC"
print*,"AHE=",matrices%ages(3,2),"km",matrices%ages(3,1),"degC"
print*,"ZHE=",matrices%ages(4,2),"km",matrices%ages(4,1),"degC"

matrices%zc=matrices%zz
return

end subroutine initialize_parameters
