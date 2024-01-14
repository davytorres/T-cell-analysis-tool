       subroutine calculate_mean_squared_displacement(ii, i, i90ca, 
     &                                         x90,y90,
     &                                        x90a, 
     &                                                 y90a,
     &                                                     y90_2a, 
     &                                      icount,
     &                                       max_tcells_allfiles,
     &              max_frames_allfiles,
     &                                               jfirstaa,
     &                                               jlastaa,
     &                                              valid, tobv, xobv,
     &                                            yobv, zobv, ym,
     &                                            xtcell,
     &                                                   xtcell2, 
     &                                               ytcell,
     &                                               xytcell,
     &                                               rtcell, slopefit,
     &                                               sdfit, i90c,
     &                                               y90list,y90lista)

       implicit none

       integer ii, i
       integer max_tcells_allfiles
       integer max_frames_allfiles
       integer i90c(1000)
       integer i90ca(1000)
       integer icount(max_tcells_allfiles)
       integer jfirstaa(max_tcells_allfiles,1000)
       integer jlastaa(max_tcells_allfiles,1000)
       integer number_tcells, numcellsz
       double precision x90(1000), y90(1000)
       double precision x90a(1000), y90a(1000), y90_2a(1000)
       double precision tobv(max_tcells_allfiles,max_frames_allfiles)
       double precision xobv(max_tcells_allfiles,max_frames_allfiles)
       double precision yobv(max_tcells_allfiles,max_frames_allfiles)
       double precision zobv(max_tcells_allfiles,max_frames_allfiles)
       double precision time_msd, speed_c, speed_d
       double precision ym
       double precision xtcell, xtcell2, ytcell, xytcell, rtcell
       double precision slopefit, sdfit
       double precision y90list(50,10000), y90lista(10000)
       logical valid


      
       integer j,k, jk
       integer ijka, jki
       integer ja, ka 
       integer jja
       integer jfirst, jlast
       integer i90
       integer rinum
       integer iii, jjj
       double precision distance_msd, distancez, distancef, distanceff
       double precision eltime
       double precision sumx, sumxx, sumy, sumyy, sumxy
       double precision axcell, aycell
       double precision botmsd, slopeindv
       double precision y90alist(50,10000)


                  do ijka = 1,1000
                     i90ca(ijka) = 0
                     x90a(ijka) = 0.d0
                     y90a(ijka) = 0.d0
                     y90_2a(ijka) = 0.d0
                  end do
                  y90alist = 0.d0

                  do jki = 1,icount(i)

                     ja = jfirstaa(i,jki)
                     ka = jlastaa(i,jki)

                     do j = ja,ka-1
                        do k = j+1,ka

                           valid = .true.
                           do jk = j,k
                              if (tobv(i,jk) .le. -1.d0) then
                                 valid = .false.
                              end if
                           end do

                           if (valid) then
                              jfirst = j
                              jlast = k

                              distance_msd =
     &                        (xobv(i,jlast) - xobv(i,jfirst))**2 +
     &                        (yobv(i,jlast) - yobv(i,jfirst))**2 +
     &                        (zobv(i,jlast) - zobv(i,jfirst))**2
                           
                              distancez = 0.d0
                              do jja = jfirst+1,jlast

                                 distancef =
     &                           (xobv(i,jja) - xobv(i,jja-1))**2 +
     &                           (yobv(i,jja) - yobv(i,jja-1))**2 +
     &                           (zobv(i,jja) - zobv(i,jja-1))**2

                                 distancez = distancez + sqrt(distancef)



                                 eltime =
     &                           (tobv(i,jja) - tobv(i,jfirst))/60.d0


                                 distanceff =
     &                           (xobv(i,jja) - xobv(i,jfirst))**2 +
     &                           (yobv(i,jja) - yobv(i,jfirst))**2 +
     &                           (zobv(i,jja) - zobv(i,jfirst))**2

c                                 distanceff = dsqrt(distanceff)
                     
                                 i90 = nint(60.d0*eltime/90.d0)

                                 if (i90 .gt. 0) then
                                    x90a(i90) = x90a(i90) + eltime
                                    y90a(i90) = y90a(i90) +
     &                              distanceff
                                    y90_2a(i90) = y90_2a(i90) +
     &                              distanceff**2
                                    i90ca(i90) = i90ca(i90) + 1


                                    y90alist(i90,i90ca(i90)) =
     &                              distanceff
                           
c                                  x90(i90) = x90(i90) + eltime
c                                  y90(i90) = y90(i90) + distanceff
c                                  y90_2(i90) = y90_2(i90) + distanceff**2
cc                                 z90(i90) = z90(i90) + speed_c
cc                                 z90_2(i90) = z90_2(i90) + speed_c**2                     
c                                  i90c(i90) = i90c(i90) + 1
                                end if

                              end do
                              time_msd = tobv(i,jlast) - tobv(i,jfirst)

                              time_msd = time_msd/60.d0
                              speed_d = distance_msd
                              speed_c = distancez/time_msd

                              


                           end if

                        end do
                     end do



                     
                  end do   

                  sumx = 0.d0
                  sumxx = 0.d0
                  sumy = 0.d0
                  sumxy = 0.d0
                  rinum = 0
                  
c                  do iii = 1,60
                  do iii = 1,7
                     if (i90ca(iii) .gt. 0) then

                        do jjj = 1,i90ca(ii)
                           y90lista(jjj) = y90alist(iii,jjj)
c                        print*,'y90lista(',jjj,') = ',y90lista(jjj),iii
                        end do
                        call median(y90lista,ym,i90ca(ii),10000)
c                        if (ym .gt. 0.d0) then
                        x90a(iii) = x90a(iii)/dble(i90ca(iii))                        
c                        print*,'ym = ',iii,ym
                        y90a(iii) = y90a(iii)/dble(i90ca(iii))
c                       y90a(iii) = ym

                        y90_2a(iii) = y90_2a(iii)/dble(i90ca(iii))
                        axcell = log10(x90a(iii))
                        aycell = log10(y90a(iii))
c                        write(530,*) axcell,aycell
                        write(530,*) x90a(iii),y90a(iii)                        
                        
                        sumx = sumx + axcell
                        sumy = sumy + aycell
                        xtcell = xtcell + axcell
                        ytcell = ytcell + aycell
                        sumxx = sumxx + axcell**2
                        xtcell2 = xtcell2 + axcell**2
                        sumxy = sumxy + axcell*aycell
                        xytcell = xytcell + axcell*aycell
                        rinum = rinum + 1.d0
                        rtcell = rtcell + 1.d0

                        x90(iii) = x90(iii) + x90a(iii)
                        y90(iii) = y90(iii) + y90a(iii)
                        i90c(iii) = i90c(iii) + 1.d0

                        y90list(iii,i90c(iii)) = y90a(iii)
c                        end if
                     
                     end if
                  end do

                  write(530,*)

                  rinum = 1.d0/rinum
                  botmsd = sumxx - sumx*sumx*rinum
                  if (abs(botmsd) .gt. 1.d-6) then
                     slopeindv = (sumxy - sumx*sumy*rinum)/
     &                            botmsd
                     write(529,*) slopeindv
                     if (dabs(slopeindv - slopefit) .lt. sdfit .and.
     &                   rinum .lt. .15 .and.
     &                   y90a(7) .gt. 500.d0 .and.
     &                   y90a(7) .lt. 1000.d0 ) then
                        do iii = 1,7
                           write(983,*)
     &                  log10(x90a(iii)),log10(y90a(iii)),slopeindv
                        end do
                        write(983,*) 
                     end if
                  end if

       return
       end
                  
