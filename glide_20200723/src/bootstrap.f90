program bootstrap

use definitions 
implicit none

!include 'mpif.h'

! definition of each type are found in module_definitions.f90

type (matr) matrices
type (parm) params

integer i,ierr,taskid,numtasks

real diff,misfit

!call MPI_INIT(ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
!call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

open(3,file="results.txt",status="unknown")

do i=1,45
  !print*,"calling glide",i
  taskid=i
   call glide(matrices,params,taskid,diff,misfit)
  !print*,"glide called"
  write(3,*) i,diff,misfit
enddo


!call MPI_FINALIZE(ierr)

end

