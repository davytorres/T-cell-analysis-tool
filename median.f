      subroutine median(x,xm,n,nmax)
      double precision x(nmax),xs(nmax),xm
      integer n,n2,nmax

      call sort(x,xs,n,nmax)
      if (mod(n,2) .eq. 0) then
         n2 = n/2
         xm = 0.5d0*(xs(n2) + xs(n2 + 1))
      else
         xm = xs((n+1)/2)
      end if

      return
      end