program reformatfault

double precision x_b,x_e,y_b,y_e
open(8,file='block_boundaries.txt',status='old')
open(7,file='block_boundariesGMT.txt',status='unknown')

do i=1,10060
  read(8,*) x_b,y_b,x_e,y_e
  write (7,*)x_b,y_b
  write (7,*)x_e,y_e
  write (7,'(a)')'>'
enddo

print*,'done writing'
end program
