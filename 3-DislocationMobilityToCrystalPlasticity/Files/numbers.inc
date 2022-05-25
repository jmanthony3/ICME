c
c=====================================================================72
c
c------- useful numerical constants
c
      real*8  pzero, pone, ptwo, pthree 
      real*8  pthird, phalf, ptwothrd 
      real*8  sqr23, sqr32, sqr2, sqr3

      real*8  Ident2nd(9)
      real*8  Ident(6), factoTwo(6), factoHalf(6)
  
      data pzero    /0.d0/
      data pone     /1.d0/
      data ptwo     /2.d0/
      data pthree   /3.d0/

      data pthird   /0.333333333333333d0/
      data phalf    /0.5d0/
      data ptwothrd /0.666666666666667d0/

      data sqr23    /0.816496580927726d0/
      data sqr32    /1.224744871391589d0/
      data sqr2     /1.414213562373095d0/
      data sqr3     /1.732050807568877d0/

      data Ident2nd  /1.d0, 0.d0, 0.d0, 
     &                0.d0, 1.d0, 0.d0,
     &                0.d0, 0.d0, 1.d0/

      data Ident     /1.d0, 1.d0, 1.d0, 0.d0,  0.d0,  0.d0/
      data factoTwo  /1.d0, 1.d0, 1.d0, 2.d0,  2.d0,  2.d0/
      data factoHalf /1.d0, 1.d0, 1.d0, 0.5d0, 0.5d0, 0.5d0/
c
c=====================================================================72
c
