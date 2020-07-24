program create_grid

use definitions
use omp_lib

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
character :: line*1024,datafile*100
integer i,k,kk,nx,ny,j,nx_dum,ny_dum
integer ii,jj,skip,nx0,ny0
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

params%sigma2=params%sigma2**2.
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
 read(33,*,end=500) lo,la,tmp1,tmp2,tmp3,tmp4
 if (lo.gt.180.) lo = lo - 360.
 if (tmp2.gt.params%t_total.or. &
     lo.lt.params%lon1.or.lo.gt.params%lon2.or. &
     la.lt.params%lat1.or.la.gt.params%lat2) k=k-1
enddo
500 continue
params%n=k-1
rewind(33)

print*,"total ages",params%n,"control points=",params%dummy

!allocate arrays
allocate (matrices%ta(params%n),&
          matrices%a_error(params%n),&
          matrices%zc(params%n),&
          matrices%zcp(params%n),&
          matrices%x(params%n),&
          matrices%elev(params%n),&
          matrices%x_true(params%n),&
          matrices%y_true(params%n),&
          matrices%misfits(params%n),&
          matrices%syn_age(params%n),&
          matrices%isys(params%n),matrices%depth1(params%n),&
          matrices%y(params%n),matrices%depth(params%n),matrices%zz(params%n))

k=0
do i=1,20000
       k=k+1
       read(33,*,end=300) lo,la,tmp1,tmp2,tmp3,sys  
       if (lo.gt.180.) lo=lo-360.
       if (tmp2.gt.params%t_total.or. &
       lo.lt.params%lon1.or.lo.gt.params%lon2.or. &
       la.lt.params%lat1.or.la.gt.params%lat2) then     
       k=k-1
     else
       matrices%x(k)=lo
       matrices%y(k)=la
       matrices%elev(k)=tmp1
       matrices%ta(k)=tmp2
       matrices%a_error(k)=tmp3
       matrices%isys(k)=sys
       if (k.eq.params%n) exit
     endif
enddo
300 continue
close(33)
k=k-1

print*,'creating dummy points'
k=0
do i=1,nx_dum
  do j=1,ny_dum
    k=k+1
     xx_dum=params%lon1 + float((i-1))*xl
     yy_dum=params%lat1 + float((j-1))*yl
     kk=0
     do jj=1,params%n
      dist=((xx_dum-matrices%x(jj))**2+(yy_dum-matrices%y(jj))**2)
      if (dist.lt.(xl*10.*xl*10.)) then
        kk=1
        goto 982
      endif
     enddo
982    if (kk.eq.0) k=k-1
enddo
enddo
params%dummy=k
allocate(matrices%x_dum(params%dummy),matrices%y_dum(params%dummy))
allocate(matrices%x_dum_true(params%dummy),&
         matrices%y_dum_true(params%dummy))
k=0
 do i=1,nx_dum
  do j=1,ny_dum
    k=k+1
     xx_dum=params%lon1 + float((i-1))*xl
     yy_dum=params%lat1 + float((j-1))*yl
     kk=0
     do jj=1,params%n
      dist=((xx_dum-matrices%x(jj))**2+(yy_dum-matrices%y(jj))**2)
      if (dist.lt.(xl*10.*xl*10.)) then
        kk=1
        goto 983
      endif
     enddo
983    if (kk.eq.0) then
        k=k-1
       else
        matrices%x_dum(k) = xx_dum
        matrices%y_dum(k) = yy_dum
       endif
enddo
enddo
print*,'number of dummy points=',params%dummy

!create grid file
open(872,file='grid_noID.txt',status='unknown')
do i=1,params%dummy
  write(872,*) matrices%x_dum(i),matrices%y_dum(i)
enddo


! positions are later transformed, keep lat lon values
!matrices%x_dum_true=matrices%x_dum
!matrices%y_dum_true=matrices%y_dum
!matrices%x_true=matrices%x
!matrices%y_true=matrices%y

!and convert them into local coordinate system
colat=cos(((params%lat1+params%lat2)/2.)*pi/180.)

!convert points to local grid
!do i=1,params%dummy
!  matrices%x_dum(i) = (matrices%x_dum(i)-params%lon1)*111.11*colat
!  matrices%y_dum(i) = (matrices%y_dum(i)-params%lat1)*111.11 
!enddo

!open topography, should be a bit bigger than data extents

!xmin=dble(minval(lon))
!xmax=dble(maxval(lon))
!ymin=dble(minval(lat))
!ymax=dble(maxval(lat))
!print*,xmin,xmax,ymin,ymax

!convert to local coordinate system  
!do i=1,params%n
! matrices%x(i) = (matrices%x_true(i)-params%lon1)*111.11*colat
! matrices%y(i) = (matrices%y_true(i)-params%lat1)*111.11
!enddo

end
