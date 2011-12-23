      subroutine bsplvb ( t, lent,jhigh, index, x, left, biatx )
c      implicit none
c     -------------

calculates the value of all possibly nonzero b-splines at  x  of order
c
c               jout  =  dmax( jhigh , (j+1)*(index-1) )
c
c  with knot sequence  t .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  left + jout  , assumed to be nonde-
c        creasing.
c    a s s u m p t i o n  :  t(left)  <  t(left + 1)
c    d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
c
c  jhigh,
c  index.....integers which determine the order  jout = max(jhigh,
c        (j+1)*(index-1))  of the b-splines whose values at  x  are to
c        be returned.  index  is used to avoid recalculations when seve-
c        ral columns of the triangular array of b-spline values are nee-
c        ded (e.g., in  bvalue  or in  bsplvd ). precisely,
c                     if  index = 1 ,
c        the calculation starts from scratch and the entire triangular
c        array of b-spline values of orders 1,2,...,jhigh  is generated
c        order by order , i.e., column by column .
c                     if  index = 2 ,
c        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
c        nerated, the assumption being that  biatx , j , deltal , deltar
c        are, on entry, as they were on exit at the previous call.
c           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
c        the next column of b-spline values is generated.
c
c  w a r n i n g . . .  the restriction   jout <= jmax (= 20)  is
c        imposed arbitrarily by the dimension statement for  deltal and
c        deltar  below, but is  n o w h e r e  c h e c k e d  for .
c
c  x.....the point at which the b-splines are to be evaluated.
c  left.....an integer chosen (usually) so that
c                  t(left) <= x <= t(left+1)  .
c
c******  o u t p u t  ******
c  biatx.....array of length  jout , with  biatx(i)  containing the val-
c        ue at  x  of the polynomial of order  jout  which agrees with
c        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
c        t(left+1)) .
c
c******  m e t h o d  ******
c  the recurrence relation
c
c                       x - t(i)               t(i+j+1) - x
c     b(i,j+1)(x)  =  ----------- b(i,j)(x) + --------------- b(i+1,j)(x)
c                     t(i+j)-t(i)             t(i+j+1)-t(i+1)
c
c  is used (repeatedly) to generate the
c  (j+1)-vector  b(left-j,j+1)(x),...,b(left,j+1)(x)
c  from the j-vector  b(left-j+1,j)(x),...,b(left,j)(x),
c  storing the new values in  biatx  over the old.  the facts that
c            b(i,1) = 1         if  t(i) <= x < t(i+1)
c  and that
c            b(i,j)(x) = 0  unless  t(i) <= x < t(i+j)
c  are used. the particular organization of the calculations follows
c  algorithm (8)  in chapter x of the text.
c

C Arguments
      integer lent, jhigh, index, left
      double precision t(lent),x, biatx(jhigh)
c     dimension     t(left+jout), biatx(jout)
c     -----------------------------------
c current fortran standard makes it impossible to specify the length of
c  t  and of  biatx  precisely without the introduction of otherwise
c  superfluous additional arguments.

C Local Variables
      integer jmax
      parameter(jmax = 20)
      integer i,j,jp1
      double precision deltal(jmax), deltar(jmax),saved,term

      save j,deltal,deltar
      data j/1/
c
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1d0
      if (j .ge. jhigh)                 go to 99
c
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0d0
         do 26 i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
   26       saved = deltal(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
c
   99                                   return
      end
