      subroutine pearsonr(r,slope,x,y,n,sdx,sdy,yint,nmax)
      implicit none
      integer n,i,nmax
      double precision x(nmax),y(nmax)
      double precision sdx,sdy
      double precision sumx,sumy,sumx2,sumy2,sumxy
      double precision r,rn,rtop,rbot,slope
      double precision sigmax,sigmay
      double precision xbar,ybar,yint

      sumx = 0.d0
      sumy = 0.d0
      sumx2 = 0.d0
      sumy2 = 0.d0
      sumxy = 0.d0
      do i = 1,n
         sumx = sumx + x(i)
         sumy = sumy + y(i)
         sumx2 = sumx2 + x(i)**2
         sumy2 = sumy2 + y(i)**2
         sumxy = sumxy + x(i)*y(i)
      end do

      rn = 1.d0/dble(n)

      ybar = sumy*rn
      xbar = sumx*rn

c      print*,'sumx sumy = ',sumx,sumy,rn                                                                                            
      rtop = sumxy - sumx*sumy*rn
c      print*,'rtop = ',rtop                                                                                                         
      sigmax = sumx2 - rn*(sumx**2)
      sigmay = sumy2 - rn*(sumy**2)
      rbot = dsqrt(sigmax*sigmay)
c      print*,'rbot = ',rbot,sigmax,sigmay                                                                                           

      r = rtop/rbot
      sdx = dsqrt(sigmax*rn)
      sdy = dsqrt(sigmay*rn)

      slope = r*dsqrt(sigmay/sigmax)


      yint = ybar - slope*xbar

      return
      end