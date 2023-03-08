      subroutine volume_patrolled_tcell(number_tcells,icount,
     &     jfirstaa,jlastaa,grid, xobv, yobv, zobv, tobv,
     &     max_tcells_allfiles, max_frames_allfiles, ncube, 
     &     volumetcell,dxg, dxg2, vg, time_min, time_max, 
     &     weights_time, radius_tcell_squared)

c     Inputs
c     number_tcells: Number of T cells that will be tracked
c     icount: Number of sets for each T cell
c     jfirstaa: First frame in set
c     jlastaa: Last frame in set
c     grid: Three-dimensional grid that the T cell moves within
c     tobv: timestamp of of when a T cell was observed at a 
c           certain frame
c     xobv: x poistion of a tcell at a certain frame
c     yobv: y poistion of a tcell at a certain frame
c     zobv: z poistion of a tcell at a certain frame
c     ncube: memory allocated for grid
c     time_min: minimum amount of allowed time
c     time_max: maximum amount of allowed time
c     radius_tcell_squared: radius of tcell sqaured
c     weights_time: weighted time based on time spent in frame
c     vg: volume of a grid cell
c     dxg: one edge of the grid cell
c     dxg2: half of dxg
      implicit none
        

      integer number_tcells
      integer ncube
      integer max_tcells_allfiles
      integer max_frames_allfiles
      integer icount(max_tcells_allfiles)
      integer jfirstaa(max_tcells_allfiles,1000)
      integer jlastaa(max_tcells_allfiles,1000)
      double precision time_min, time_max
      double precision xobv(max_tcells_allfiles, max_frames_allfiles)
      double precision yobv(max_tcells_allfiles, max_frames_allfiles)
      double precision zobv(max_tcells_allfiles, max_frames_allfiles)
      double precision tobv(max_tcells_allfiles, max_frames_allfiles)
      logical grid(ncube,ncube,ncube)
      double precision dxg, dxg2
      double precision weights_time(max_tcells_allfiles, 1000)
      double precision radius_tcell_squared
      double precision vg
      
    

      
c     Output
c     volumetcell: volume patrolled per time for each T cell
      double precision volumetcell(max_tcells_allfiles)
c

c     What is missing from the subroutine interface?
c
c     time_vol,time_min,time_max, tobv  ??????

      
c     Local variables
      integer i, jki, j, ijk
      integer ja, ka
      integer ig, jg, kg
      double precision sst,xxxi,yyyi,zzzi
      double precision time_lapse
      double precision time_vol(max_tcells_allfiles)
      integer iig, jjg, kkg
      double precision xc, yc, zc
      double precision distance2
      double precision vtcell
     

      
c     Calculate the volume traversed by each T cell
      do i = 1,number_tcells
                  
         time_vol(i) = 0.d0
         volumetcell(i) = 0.d0

c                  print*,'i icount = ',i,icount(i)
                  do jki = 1,icount(i)

                     ja = jfirstaa(i,jki)
                     ka = jlastaa(i,jki)
                     grid = .false.

                     time_lapse = tobv(i,ka)-tobv(i,ja)
                     if (time_lapse .gt. time_min .and.
     &                   time_lapse .lt. time_max) then
                     
c                       do j = 1,number_frames-1
                        do j = ja,ka-1

                           if (tobv(i,j) .gt. -1.d-6 .and.
     &                         tobv(i,j+1) .gt. -1.d-6) then
c                             (iig,jjg,kkg) are the logical coordinates of the T cell

                              time_vol(i) = time_vol(i) +
     &                        tobv(i,j+1) - tobv(i,j)

                              do ijk = 0,100

c                                Create a line of 101 points between xobv(i,j) and xobv(i,j+1)
                                 sst = dble(ijk)/100.d0
                                 xxxi = xobv(i,j) +
     &                                  sst*(xobv(i,j+1)-xobv(i,j))
                                 yyyi = yobv(i,j) +
     &                                  sst*(yobv(i,j+1)-yobv(i,j))
                                 zzzi = zobv(i,j) +
     &                                  sst*(zobv(i,j+1)-zobv(i,j))                           
                           
                        
                                 iig = (xxxi + 10.d0)/dxg
                                 jjg = (yyyi + 10.d0)/dxg
                                 kkg = (zzzi + 10.d0)/dxg

                                 if (iig + 2 .gt. ncube .or.
     &                               jjg + 2 .gt. ncube .or.
     &                               kkg + 2 .gt. ncube) then
                                    print*,'Increase value of ncube'
                                    stop
                                 end if

                                 if (iig - 2 .lt. 1 .or.
     &                               jjg - 2 .lt. 1 .or.
     &                               kkg - 2 .lt. 1) then
                                    print*,'Increase value of ncube'
                                    stop
                                 end if

                                 do ig = iig-2,iig+2
                                    do jg = jjg-2,jjg+2
                                       do kg = kkg-2,kkg+2
                                          xc = dxg*dble(ig)-10.d0+dxg2
                                          yc = dxg*dble(jg)-10.d0+dxg2
                                          zc = dxg*dble(kg)-10.d0+dxg2
                                          distance2 = (xxxi - xc)**2 +
     &                                                (yyyi - yc)**2 +
     &                                                (zzzi - zc)**2
                                          if (distance2 .lt. 
     &                                       radius_tcell_squared) then
                                             grid(ig,jg,kg) = .true.
                                          end if
                                       end do
                                    end do
                                 end do
                           
                              end do

                           else
                              print*,'ERROR in continuous tracks',i
                              stop
                           end if

                        end do
                     else
                        print*,'Error'
                        stop
                     end if

                     vtcell = 0.d0
                     do kkg = 1,ncube
                        do jjg = 1,ncube
                           do iig = 1,ncube
                              if (grid(iig,jjg,kkg)) then
                                 vtcell = vtcell + vg
                              end if
                           end do
                        end do
                     end do

                     volumetcell(i) = volumetcell(i) +
     &               (vtcell/time_lapse)*weights_time(i,jki)
                     
                  end do
c                  if (icount(i) .gt. 0) then
c                     volumetcell(i) = volumetcell(i)/icount(i)
c     end if

      end do
      end
