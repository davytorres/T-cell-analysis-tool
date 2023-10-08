       
        subroutine identify_tracks(ii, number_tcells, icount, 
     &                             max_tcells_allfiles, max_frames_allfiles, 
     &                             maxnumframes, tobv, time_min, time_max, time_opt, 
     &                             jfirstaa, jlastaa, total_time, weights_time)
 
      implicit none
      integer ii, number_tcells, max_tcells_allfiles, max_frames_allfiles
      integer jfirstaa(max_tcells_allfiles,1000), jlastaa(max_tcells_allfiles,1000)
      integer icount(max_tcells_allfiles)
      integer maxnumframes(max_tcells_allfiles)
      double precision tobv(max_tcells_allfiles,max_frames_allfiles)
      double precision total_time(max_tcells_allfiles)
      double precision weights_time(max_tcells_allfiles,1000)

      integer i, j, k, jki, ja, ka
      integer kstart, jfirst, jlast, jlast_keep
      logical continuetrack
      double precision diff_time_min, diff_time
      double precision time_lapse, time_opt
      double precision time_lapsej(1000)

      do i = 1,number_tcells
           icount(i) = 0
           kstart = 1
           continuetrack = .true.
           do while (continuetrack)
                jfirst = 0
                k = kstart
                do while (jfirst .eq. 0 .and. k .lt. maxnumframes(i))
                     if (tobv(i,k) .gt. -1.d-6) then
                          jfirst = k
                     end if
                     k = k + 1
                end do

                jlast = 0
                if (jfirst .gt. 0) then
                     k = jfirst+1
                     do while (jlast .eq. 0 .and. k .le. maxnumframes(i))
                          if (tobv(i,k) .le. -1.d-6) then
                               jlast = k-1
                          else
                               k = k + 1
                          end if
                     end do
                end if

                if (jfirst .gt. 0) then
                     if (k .gt. maxnumframes(i)) then
                          if (tobv(i,maxnumframes(i)) - tobv(i,jfirst) .gt. time_min) then
                               diff_time_min = 1.d+30
                               do jki = jfirst+1, maxnumframes(i)
                                    time_lapse = tobv(i,jki) - tobv(i,jfirst)
                                    diff_time = abs(time_opt - time_lapse)
                                    if (diff_time .lt. diff_time_min) then
                                         diff_time_min = diff_time
                                         jlast_keep = jki
                                    end if
                               end do

                               if (tobv(i,jlast_keep) - tobv(i,jfirst) .lt. time_max) then
                                    icount(i) = icount(i) + 1
                                    if (icount(i) .gt. 1000) then
                                         print*,'Increase 2nd arg jfirstaa'
                                         stop
                                    end if
                                    jfirstaa(i,icount(i)) = jfirst
                                    jlastaa(i,icount(i)) = jlast_keep
                               else
                                    print*,'Error in continuous tracks A'
                                    print*,tobv(i,jlast_keep) - tobv(i,jfirst)
                                    print*,time_max,i,ii
                                    stop
                               end if
                               if (jlast_keep .ge. maxnumframes(i)) then
                                    continuetrack = .false.
                               else
                                    kstart = jlastaa(i,icount(i)) + 1
                               end if
                          else
                               continuetrack = .false.
                          end if
                     else
                          if (jlast .eq. jfirst) then
                               kstart = jfirst + 2
                          else
                               if (tobv(i,jlast) - tobv(i,jfirst) .gt. time_min) then
                                    diff_time_min = 1.d+20
                                    do jki = jfirst+1, jlast
                                         time_lapse = tobv(i,jki) - tobv(i,jfirst)
                                         diff_time = abs(time_opt - time_lapse)
                                         if (diff_time .lt. diff_time_min) then
                                              diff_time_min = diff_time
                                              jlast_keep = jki
                                         end if
                                    end do

                                    if (tobv(i,jlast_keep) - tobv(i,jfirst) .lt. time_max) then
                                         icount(i) = icount(i) + 1
                                         if (icount(i) .gt. 1000) then
                                              print*,'Increase 2nd arg of jfirstaa'
                                              stop
                                         end if
                                         jfirstaa(i,icount(i)) = jfirst
                                         jlastaa(i,icount(i)) = jlast_keep
                                         kstart = jlastaa(i,icount(i)) + 1
                                    else
                                         print*,'Error in continuous tracks B'
                                         print*,tobv(i,jlast_keep) - tobv(i,jfirst)
                                         print*,time_max
                                         stop
                                    end if
                               else
                                    kstart = jlast + 1
                               end if
                          end if
                     end if
                else
                     continuetrack = .false.
                end if
           end do
      end do

      do i = 1,number_tcells
           total_time(i) = 0.d0
           do jki = 1,icount(i)
                ja = jfirstaa(i,jki)
                ka = jlastaa(i,jki)
                time_lapsej(jki) = tobv(i,ka) - tobv(i,ja)
                if (time_lapsej(jki) .lt. 150.d0) then
                     time_lapsej(jki) = 0.d0
                end if
                total_time(i) = total_time(i) + time_lapsej(jki)
           end do
           if (total_time(i) .gt. 0.d0) then
                do jki = 1,icount(i)
                     weights_time(i,jki) = time_lapsej(jki)/total_time(i)
                end do
           else
                do jki = 1,icount(i)
                     weights_time(i,jki) = 0.d0
                end do
           end if
           wsum = 0.d0
           do jki = 1,icount(i)
                wsum = wsum + weights_time(i,jki)
           end do
      end do

      return
      end 

