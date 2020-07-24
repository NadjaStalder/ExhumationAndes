subroutine solve_inverse(matrices,params,iter,residu)

! In this subroutine the maximum likelihood estimates of the exhumation rates are obtained.

! edot = edot_pr + CG'(GCG'+Cee) (log(Zc) - log(Aedot_pr))

! First ACA'+Cee is computed adn stored in Y1
! Y1 is inverted using LUD and backsubstitution

! G'(GCG'+Cee) is stored in Y2
! zz stores (log(Zc) - log(Aedot_pr))
 
! The total CG'(GCG'+Cee) (log(Zc) - log(Aedot_pr)) is stored in BB

! Next the dummy points are computed, for the dummy points .....

! edot = edot_pr + C~^G'(ACG'+Cee) (log(Zc) - log(Aedot_pr))
!~ represents a parameter at a control point
!^ represents a parameter at a data point
!Cmm is a row of C^~ dscribing covariance of control point with data points

! H is the inverse operator used to calculate the exhumation rates at the data points
! this is used in other subroutines later

use definitions
use omp_lib

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
integer i,j,k,INFO,ik,jk,iter,iflag
integer id,nproc,chunk,first,last,m
double precision dist,residu,resid,resim
double precision,dimension(:),allocatable::Cmm,zp,xresi,vresi
double precision,dimension(:,:),allocatable:: covinv,test,cdinv
double precision:: rest,summ
double precision:: wa

allocate(matrices%Y1(params%n,params%n))
allocate(matrices%B(params%n),matrices%II(params%n,params%n),matrices%work(1000*params%n)) 
allocate(cdinv(params%n,params%n))

! computes Y1=GCG'+Cee
! firstly it computes the exhumation at the dummy, and then at data points

! note Y2=CG'
!print*,"calculating Y2..."

!matrices%Y2=0.d0
!call dgemm('N','T',params%n*params%m_max,params%n,params%n*params%m_max,1.d0,matrices%cov, & 
!           params%n*params%m_max,matrices%G,params%n,0.d0,matrices%Y2,params%n*params%m_max)

!matrices%Y1=0.d0
!call dgemm('N','N',params%n,params%n,params%n*params%m_max,1.d0,matrices%G,params%n,matrices%Y2,& 
!           params%n*params%m_max,0.d0,matrices%Y1,params%n)

cdinv=1.d-9
do i=1,params%n
 params%sigmaD=max(matrices%a_error(i)/matrices%ta(i),0.4d0)
 !params%sigmaD=max(2.5d0/matrices%ta(i),0.6d0)
 params%sigmaD=.4d0
 matrices%Y1(i,i)=matrices%Y1(i,i)+(params%sigmaD)**2
 cdinv(i,i)=1.d0/(params%sigmaD)**2
enddo

!  Use segtrf+sgetri to invert matrix
!  and find (Y1)**-1

matrices%B=0

!lu decomposing the matrix
!call dgetrf(params%n,params%n,matrices%Y1,params%n,matrices%B,INFO)
!to calculate the inverse
!call dgetri(params%n,matrices%Y1,params%n,matrices%B,matrices%work,1000*params%n,INFO)

!! it computes Y2=G'Cd**(-1)
matrices%Y2=0.d0
call dgemm('T','N',params%n*params%m_max,params%n,params%n,1.d0,matrices%G,params%n, &
           cdinv,params%n,0.d0,matrices%Y2,params%n*params%m_max)

!apriori age (i.e. A*exp(e))
matrices%zz=0.d0
call dgemv('N',params%n,params%n*params%m_max,1.d0,matrices%A,params%n,matrices%edot_pr,1,0.d0,matrices%zz,1)
!matrices%zz=log(matrices%zc+matrices%elev)-log(matrices%zz)
iflag=0
do i=1,params%n
if ((matrices%zc(i)+matrices%elev(i)).lt.0.) then
print*,'negative values',matrices%elev_true(i),matrices%ta(i),matrices%a_error(i)
!matrices%elev(i)=matrices%elev(i)+0.2
iflag=1
endif
enddo
if (iflag.eq.1) stop
matrices%zz = log(matrices%zz) -log(matrices%zc + matrices%elev)
!do i=1,params%n
!print*,i,matrices%zz(i),matrices%ta(i)
!enddo

!BB=Y2*(log(zc)-log(A*eprior))
allocate(matrices%BB(params%n*params%m_max))

matrices%BB=0.d0
call dgemv('N',params%n*params%m_max,params%n,1.d0,matrices%Y2,params%n*params%m_max,matrices%zz,1,0.d0,matrices%BB,1)
!matrices%BB=matmul(matrices%Y2,matrices%zz)

!compute edot at the control points
matrices%edot=0.d0
if (iter.gt.0) then

!$omp parallel private(i,j,k,ik,jk,id,first,last,chunk,dist,Cmm)

  id=omp_get_thread_num()
  nproc=omp_get_num_threads()
  chunk=floor((params%m_max*params%dummy)/float(nproc))

  first=(id*chunk)+1
  last=(id+1)*chunk
  if (id.eq.nproc-1) last=last+mod((params%m_max*params%dummy),nproc)
!print*,'memory allocation',id,nproc
  allocate(Cmm(params%n*params%m_max))
 
do i=first,last

  Cmm=0.d0
  ik = floor(float(i-1)/float(params%m_max)) +1
  do j=1,params%n*params%m_max
    jk = floor(float(j-1)/float(params%m_max))+1
    if (mod(i,params%m_max).eq.mod(j,params%m_max)) then
      dist=sqrt((matrices%x_dum(ik)-matrices%x(jk))**2+(matrices%y_dum(ik)-matrices%y(jk))**2)
      if (matrices%idum_block(ik).ne.matrices%iblock(jk)) dist=1.e6
      Cmm(j)=params%sigma2*exp(-(dist/params%xL))    
    if (Cmm(j).lt.1.d-8) Cmm(j)=0.d0
    endif
  enddo

  do k=1,params%n*params%m_max
    matrices%edot(i) = matrices%edot(i)+Cmm(k)*matrices%BB(k)
  enddo

! matrices%edot(i)=matrices%edot(i)+log(matrices%edot_pr_dum(i))
 matrices%edot(i)=log(matrices%edot_pr_dum(i))-params%xmu*(matrices%edot(i)+(log(matrices%edot_pr_dum(i))-log(params%edot_mean)))
 matrices%edot(i)=exp(matrices%edot(i))
 matrices%edot_pr_dum(i)=matrices%edot(i)
enddo

  deallocate(Cmm)
!$omp end parallel
endif

!where(matrices%edot.lt.0.) matrices%edot=0.
print*,"edot min/max values at dummy= ",minval(matrices%edot),"/",maxval(matrices%edot)
!print*,"edot extremes_pr",minval(matrices%edot_pr_dum),maxval(matrices%edot_pr_dum)
!where(matrices%edot.lt.0.d0) matrices%edot=0.d0

!computing the exhumation rates for the data points
!Y2=CG'
matrices%Y2=0.d0
call dgemm('N','T',params%n*params%m_max,params%n,params%n*params%m_max,1.d0,matrices%cov,&
           params%n*params%m_max,matrices%G,params%n,0.d0,matrices%Y2,params%n*params%m_max)
!H=CG'Cd**-1
matrices%H=0.d0
call dgemm('N','N',params%n*params%m_max,params%n,params%n,1.d0,matrices%Y2,&
           params%n*params%m_max,cdinv,params%n,0.d0,matrices%H,params%n*params%m_max)

!using apriori model to calcuilate depth of sample
!this used to be age*edot_m
matrices%zz=0.d0
call dgemv('N',params%n,params%n*params%m_max,1.d0,matrices%A,params%n,matrices%edot_pr,1,0.d0,matrices%zz,1)

allocate(zp(params%n))
matrices%zz = log(matrices%zz) -log(matrices%zc + matrices%elev)
zp=matrices%zz
!print*,'computing edot....'
! it computes edot=H*(log(zc)-log(ta*exp(epsilon)))
matrices%edot_dat=0.d0
call dgemv('N',params%n*params%m_max,params%n,1.d0,matrices%H,params%n*params%m_max,matrices%zz,1,1.d0,matrices%edot_dat,1)

matrices%edot_dat=log(matrices%edot_pr)-params%xmu*(matrices%edot_dat+(log(matrices%edot_pr)-log(matrices%edot_pr2)))

resid=0.d0
resim=0.d0
open(76,file=params%run//'/resi_log.txt',status='unknown')
do i=1,params%n
write(76,*) sngl(matrices%ta(i)),sngl(matrices%zc(i) + matrices%elev(i)),sngl(zp(i)),sngl(matrices%elev_true(i))
wa=1.d0
!if (matrices%isys(i).eq.1) wa=2.d0
 params%sigmaD=max(matrices%a_error(i)/matrices%ta(i),0.5d0)
resid=resid+zp(i)**2/(params%sigmaD*wa)**2
enddo
close(76)

deallocate(matrices%B,matrices%work,cdinv)
allocate(covinv(params%n*params%m_max,params%n*params%m_max),vresi(params%n*params%m_max),&
         test(params%n*params%m_max,params%n*params%m_max),&
         xresi(params%n*params%m_max),matrices%B(params%n*params%m_max),matrices%work(1000*params%n*params%m_max))
covinv=matrices%cov
do i=1,params%n*params%m_max
xresi(i)=abs(matrices%edot_dat(i)-log(matrices%edot_pr2(i)))
enddo
vresi=0.d0
call dgetrf(params%n*params%m_max,params%n*params%m_max,covinv,params%n*params%m_max,matrices%B,INFO)
call dgetri(params%n*params%m_max,covinv,params%n*params%m_max,matrices%B,matrices%work,1000*params%n*params%m_max,INFO)
call dgemv('N',params%n*params%m_max,params%n,1.d0,covinv,params%n*params%m_max,xresi,1,1.d0,vresi,1)
call dgemm('N','N',params%n*params%m_max,params%n*params%m_max,params%n*params%m_max,&
     1.d0,matrices%cov,params%n*params%m_max,covinv,params%n*params%m_max,&
     0.d0,test,params%n*params%m_max)
do i=1,params%n*params%m_max
  resim=resim+xresi(i)**2/params%sigma2
enddo
deallocate(covinv,vresi,xresi,test)

print*,"residuals",sngl(resim+resid),"resi data",sngl(resid),"resimod",sngl(resim)
residu=(resim+resid)/(params%n)
open(947,file=params%run//'/residualsLOG.txt',status='unknown')
close(947)
open(947,file=params%run//'/residualsLOG.txt',status='old',position='append')
write(947,*) iter,sngl(resim+resid),sngl(resid),sngl(resim)
close(947)
deallocate(zp)

! go back to linear exhumation rates
matrices%edot_dat=exp(matrices%edot_dat)
print*,"edot min/max values at data= ",minval(matrices%edot_dat),"/",maxval(matrices%edot_dat)
matrices%edot_pr=matrices%edot_dat

deallocate(matrices%Y1,matrices%B,matrices%II,matrices%BB,matrices%work)

! update G
do i=1,params%n
  m=1
  k=1
  summ=0.
    !calculate the number of time steps required and the remainder term
    do while (summ.lt.matrices%ta(i))
      summ = summ + matrices%tsteps(m)
      m=m+1
    enddo
    m=m-1
    if (m.eq.0) then
       summ = 0.
       rest = 0.
    else
       summ = summ - matrices%tsteps(m)
       rest = matrices%ta(i) - summ
    endif
    k=0
    do j=(params%m_max-m)+(i-1)*(params%m_max)+1,(i-1)*(params%m_max-1)+params%m_max+(i-1)
      k=k+1
      matrices%A(i,j)=matrices%tsteps(m-k+1)
      if (k.eq.1) matrices%A(i,j)=rest
      matrices%G(i,j)=matrices%edot_pr(j)*matrices%A(i,j)/(matrices%zc(i)+matrices%elev(i))
    enddo
enddo

return
end subroutine solve_inverse

