      function secnds(t0)
      external clock
      character*8 it
      call clock(it)
      read(it(1:2),10) ihh
      read(it(4:5),10) imm
      read(it(7:8),10) iss
 10   format(i2)
      secnds=3600.*ihh+60.*imm+iss-t0
      return
      end
