program foo
    character(len=1024) :: filename
    character :: run*5

run='RUN01'
!open(55,file=run//"/stuff/edot2"/j,status="unknown")
   
    write (filename, "(A5,A13,I2)") run,"/stuff/edot2_", 1

open(55,file=trim(filename),status="unknown")
    print *, trim(filename)
close(55)
end program
