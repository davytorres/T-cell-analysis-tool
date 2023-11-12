       subroutine calculate_voulme_traversed(time_vol, volumetcell,
      &                                      max_tcells_allfiles, 
      &                                      max_frames_allfiles icount,
      &                                      jfirstaa, jlastaa, grid,
      &                                      ncube, time_lapse, tobv,
      &                                      time_min, time_max, dxg,
      &                                      dxg2, radius_tcell_squared,
      &                                      numcellsz)

        implicit none

        integer max_tcells_allfiles, max_frames_allfiles
        integer icount(max_tcells_allfiles)
        integer jfirstaa(max_tcells_allfiles,1000)
        integer jlastaa(max_tcells_allfiles,1000)
        integer ncube
        integer numcellsz
        double precision time_vol(max_tcells_allfiles) 
        double precision volumetcell(max_tcells_allfiles)
        double precision tobv(max_tcells_allfiles, max_frames_allfiles)
        double precision time_lapse
        double precision time_min, time_max
        double precision dxg, dxg2
        double precision radius_tcell_squared
        logical grid(ncube,ncube,ncube)

        integer i, jki, ijk
        integer ig, jg, kg
        integer ja, ka
        double precision sst, xxxi, yyyi, zzzi
        double precision iih, jjg, kkg
        double precision xc, yc, zc, distance2
        double precision vtcell
        

c              Calculate the volume traversed by each T cell
               do i = 1,number_tcells
              
c                 Time that cell is tracked & volume that is patrolled
c                 by T cell
                  time_vol(i) = 0.d0
                  volumetcell(i) = 0.d0

c                 print*,'i icount = ',i,icount(i)
                  do jki = 1,icount(i)

c                    First and last frames of set jki
                     ja = jfirstaa(i,jki)
                     ka = jlastaa(i,jki)
                     grid = .false.

c                    Time between first frame and last frame
                     time_lapse = tobv(i,ka)-tobv(i,ja)
                     if (time_lapse .gt. time_min .and.
     &                   time_lapse .lt. time_max) then
                     
c                       do j = 1,number_frames-1
                        do j = ja,ka-1

c                          Make sure frames j and j+1 are valid with no
c                          negative time stamp
                           if (tobv(i,j) .gt. -1.d-6 .and.
     &                         tobv(i,j+1) .gt. -1.d-6) then
c                             (iig,jjg,kkg) are the logical coordinates of the T cell

c                             Sum up time differences between frame j
c                             and j+1 
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
c                  end if


                  if (icount(i) .gt. 0) numcellsz = numcellsz + 1
