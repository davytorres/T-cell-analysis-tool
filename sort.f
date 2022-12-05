      subroutine sort(x,xs,n,nmax)
      implicit none
      double precision minx
      double precision x(nmax),xs(nmax)
      logical sorted(nmax)
      integer i,j,jmin,n,nmax
      do i = 1,n
         sorted(i) = .false.
      end do

      do i = 1,n
         minx = 1.d+20
         do j = 1,n
            if (x(j) .lt. minx .and.
     &         .not. sorted(j)) then
               minx = x(j)
               jmin = j
            end if
         end do
         xs(i) = x(jmin)
         sorted(jmin) = .true.
      end do

      return
      end