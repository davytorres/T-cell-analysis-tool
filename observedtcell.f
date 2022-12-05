
      subroutine observedtcell(kij,ii,max_tcells_allfiles,
     &                         max_frames_allfiles,
     &                         xobv,yobv,zobv,tobv,
     &                         number_tcells,number_frames,
     &                         xmin,xmax,ymin,ymax,zmin,zmax,dt,
     &                         vel,nvel,velavg,
     &                         turnang,velangle,numangles,mxnf,mxnfit,
     &                         velpersist,timepersist,npersist,
     &                         total_tcell_count,
     &                         totalbin,persistbin,probability_after,
     &                         total_num_frames,
     &                         sumi,sumi_tot,
     &                         rtc_count_all_movies,
     &                         volume_all_movies,idsum,maxnumframes)

      implicit none

      integer maxnumframesnew(max_tcells_allfiles)      
      integer maxnumframes_extend(max_tcells_allfiles)      


      integer icountangle,icountcell
      
      integer maxlinesread,mxnf
      integer sumi,sumi_tot,iii,jj
      integer total_tcell_count,kij
      parameter (maxlinesread=20000)
      integer ii,i,j,ij,mxnfit
      integer numangles(mxnf)
      integer iangle
      integer max_tcells_allfiles,max_frames_allfiles
      integer jstart,ibin,ibina
      integer number_tcells,number_frames
      integer idmax_movie
      integer idsum,numframestot
      integer maxbins,ibinc
      integer totalbin(0:300),persistbin(0:300)
      double precision probability_after(0:300)
      double precision prob_a_given_b(0:300)
      double precision ratio_prob(0:300)
      double precision prob_a,prob_b,prob_a_b
      double precision total_num_frames
      double precision timesteptot

      double precision dist,dist2,pi
      double precision xmin,xmax,ymin,ymax,zmin,zmax,dt
      double precision vec1x,vec1y,vec1z
      double precision vec1xi,vec1yi,vec1zi
      double precision vec2x,vec2y,vec2z
      double precision dot,vel2,angle,dista
      double precision velsum,timesum,velinitial
      double precision timedif,distsum
      double precision timesum2,distsum2,timea
      double precision xminz,xmaxz,yminz,ymaxz,zminz,zmaxz
      double precision rangexmax,rangeymax,rangezmax
      double precision rangexm,rangeym,rangezm
      double precision bound,disp2,disp
      double precision dispvprev,dispvafter
      double precision xdisp,ydisp,zdisp,tdisp
      double precision dispaa
      integer maxnumframes_save
      integer maxnumframes(max_tcells_allfiles)
c23456789012345678901234567890123456789012345678901234567890123456789012
      double precision xobv(max_tcells_allfiles,max_frames_allfiles)
      double precision yobv(max_tcells_allfiles,max_frames_allfiles)
      double precision zobv(max_tcells_allfiles,max_frames_allfiles)
      double precision tobv(max_tcells_allfiles,max_frames_allfiles)
      logical tused(max_tcells_allfiles,max_frames_allfiles)      

      double precision xobvc(max_tcells_allfiles,max_frames_allfiles)
      double precision yobvc(max_tcells_allfiles,max_frames_allfiles)
      double precision zobvc(max_tcells_allfiles,max_frames_allfiles)
c      double precision tobvc(max_tcells_allfiles,max_frames_allfiles)      

      double precision xobv_red(max_tcells_allfiles,max_frames_allfiles)
      double precision yobv_red(max_tcells_allfiles,max_frames_allfiles)
      double precision zobv_red(max_tcells_allfiles,max_frames_allfiles)
      double precision tobv_red(max_tcells_allfiles,max_frames_allfiles)            

      integer jfirst,jmin,jjnew,jused,jmax
      logical foundfirst
      double precision tlast,tlast_new
      double precision time_min,timed,time_diff
      double precision time_mina,time_maxa,time_diffa
      
      double precision vel(max_tcells_allfiles,max_frames_allfiles)

      double precision velxa(max_tcells_allfiles,max_frames_allfiles)
      double precision velya(max_tcells_allfiles,max_frames_allfiles)
      double precision velza(max_tcells_allfiles,max_frames_allfiles)
      double precision xpos(max_tcells_allfiles,max_frames_allfiles)
      double precision ypos(max_tcells_allfiles,max_frames_allfiles)
      double precision zpos(max_tcells_allfiles,max_frames_allfiles)
      double precision tpos(max_tcells_allfiles,max_frames_allfiles)

      double precision velavg(max_tcells_allfiles)
      double precision 
     &          velpersist(max_tcells_allfiles,max_frames_allfiles)
      double precision 
     &          timepersist(max_tcells_allfiles,max_frames_allfiles)
      double precision turnang(max_tcells_allfiles*max_frames_allfiles)
      double precision velangle(max_tcells_allfiles*max_frames_allfiles)      

      double precision eps

      double precision rnum,rtc_count
      double precision xsum,ysum,zsum
      double precision xsum2,ysum2,zsum2
      double precision xavg,yavg,zavg
      double precision xstd,ystd,zstd
      double precision xmin_box,xmax_box
      double precision ymin_box,ymax_box
      double precision zmin_box,zmax_box
      double precision volume_box,density_movie
      double precision rtc_count_all_movies
      double precision volume_all_movies
      double precision distx,disty,distz
      double precision dcell2,dcell

      logical prt
      integer icell
      integer nvel(max_tcells_allfiles)
      integer npersist(max_tcells_allfiles)

      open(unit=119,file='mean_half')
      open(unit=131,file='corrangle')
c      open(unit=115,file='directionality_flu')
c      open(unit=116,file='directionality_nflu')      
      open(unit=567,file='bayesian')
      open(unit=889,file='timestep')
      open(unit=419,file='positions_new')
      
      
      eps = 1.d-8
      
      do i = 1,max_tcells_allfiles
         do j = 1,max_frames_allfiles
            xobv(i,j) = -20.d0
            yobv(i,j) = -20.d0
            zobv(i,j) = -20.d0
            tobv(i,j) = -1.d0
            tused(i,j) = .false.
            xobvc(i,j) = -20.d0
            yobvc(i,j) = -20.d0
            zobvc(i,j) = -20.d0

            vel(i,j) = -20.d0
            
            xobv_red(i,j) = -20.d0
            yobv_red(i,j) = -20.d0
            zobv_red(i,j) = -20.d0            
            tobv_red(i,j) = -1.d0            

        end do
      end do

      print*,'readtcell',kij
      call readtcell(ii,max_tcells_allfiles,
     &               max_frames_allfiles,
     &               xobv,yobv,zobv,tobv,idmax_movie,
     &               number_tcells,number_frames,idsum,
     &               maxnumframes)
     
      idsum = idsum + idmax_movie

      xmin = 1.d+20;
      ymin = 1.d+20;
      zmin = 1.d+20;

      xmax = -1.d+20;
      ymax = -1.d+20;
      zmax = -1.d+20;

      rnum = 0.d0

      xsum = 0.d0
      ysum = 0.d0
      zsum = 0.d0

      xsum2 = 0.d0
      ysum2 = 0.d0
      zsum2 = 0.d0

      
      rangexm = -1.d+20
      rangeym = -1.d+20
      rangezm = -1.d+20


      icell = 1
      prt = .false.
      if (prt) then
         if (ii .eq. 1) then
            do i = 1,maxnumframes(icell)
               print*,'tobv(',icell,i,') = ',tobv(icell,i),
     &         xobv(icell,i),tobv(icell,i)-tobv(icell,max(i-1,1))
            end do
         end if
      end if

      do i = 1,number_tcells
c        mannumframesnew(i) keeps track of number of tracks and spacers in reduced track
         maxnumframesnew(i) = 0
         jj = 0
         jjnew = 1000001
         tlast_new = 0.d0
c         print*,'*********************** i = ',i
         maxnumframesnew(i) = 0
c23456789012345678901234567890123456789012345678901234567890123456789012
c        Reuse loop
         do while (jjnew .ge. 2) 
            tlast = tlast_new
c            print*,'jj = ',jj

c           Search for the first unused frame            
            jfirst = 0
            j = 1
            do while (j .le. maxnumframes(i) .and. jfirst .eq. 0)
c               print*,'tobv(',i,j,') = ',tobv(i,j),tused(i,j)
               if (tobv(i,j) .gt. -1.d0 .and. .not. tused(i,j)) then
                  jfirst = j
               end if
               j = j + 1
            end do
c            print*,'jfirst = ',jfirst,jj,i,icell

c           jjnew keeps tracks of how many new points are assigned in this search            
            jjnew = 0

            if (jfirst .gt. 0) then

c              Insert spacer to disconnect series of frames            
               if (tlast .gt. 0.d0) then
c                 Don't insert spacer first time - should not insert two consectutive spacers               
                  jj = jj + 1
                  maxnumframesnew(i) = jj
                  xobv_red(i,jj) = -20.d0
                  yobv_red(i,jj) = -20.d0
                  zobv_red(i,jj) = -20.d0
                  tobv_red(i,jj) = -1.d0
                  if (i .eq. icell .and. prt)
     &    print*,'tobv_red d(',i,jj,') = ',tobv_red(i,jj),xobv_red(i,jj)               
               end if

c              jj keeps track of lastest frame in new reduced track               
               jj = jj + 1
c              jused keeps track of latest frame assigned from original track               
               jused = jfirst
               jjnew = jjnew + 1
               maxnumframesnew(i) = jj
               xobv_red(i,jj) = xobv(i,jfirst)
               yobv_red(i,jj) = yobv(i,jfirst)
               zobv_red(i,jj) = zobv(i,jfirst)
               tobv_red(i,jj) = tobv(i,jfirst)+tlast

               if (i .eq. icell .and. prt)
     &              print*,'tobv_red(',i,jj,') = ',tobv_red(i,jj),
     &              tobv_red(i,jj)-tlast
               tlast_new = tobv_red(i,jj)
               tused(i,jfirst) = .true.

c              do while jfirst = 0               
               do while (jfirst .gt. 0)

                  j = jfirst+1
                  time_mina = 1.d+20
                  time_maxa = -1.d+20
c             print*,'j maxnum = ',j,maxnumframes(i),tobv(i,j),tused(i,j)
                  do while (j .le. maxnumframes(i))
                     if (.not. tused(i,j) .and.
     &                   tobv(i,j) .gt. -1.d0) then                    

                        time_diff = tobv(i,j)-tobv(i,jfirst)-90.d0
                        time_diffa = abs(time_diff)
c                        if (i .eq. icell .and. prt) then
c                           print*,tobv(i,j),tobv(i,jfirst)
c                           print*,'time_diff = ',time_diff
c                        end if
                        if (time_diffa .lt. time_mina) then
                           jmin = j
                           time_mina = time_diffa
                        end if
                        if (time_diff .gt. time_maxa) then
                           jmax = j
                           time_maxa = time_diff
                        end if
                     end if
                     j = j + 1
                  end do
c                  print*,'jmin time_mina = ',jmin,time_mina
c                  if (i .eq. icell .and. prt)
c     &     print*,'jmin = ',jmin,time_mina

c                 Assigning reduced point if                  
                  if (time_mina .le. 30.d0) then
                     jj = jj + 1
                     jused = jmin                     
                     jjnew = jjnew + 1

                     maxnumframesnew(i) = jj
                     xobv_red(i,jj) = xobv(i,jmin)
                     yobv_red(i,jj) = yobv(i,jmin)
                     zobv_red(i,jj) = zobv(i,jmin)
                     tobv_red(i,jj) = tobv(i,jmin)+tlast

                     if (i .eq. icell .and. prt) 
     &               print*,'tobv_red a(',i,jj,') = ',tobv_red(i,jj),
     &               tobv_red(i,jj)-tlast
c                    tlast is the previous max time that offsets new searches                     
c                    tlast_new keeps running track of lastest time assigned
                     tlast_new = tobv_red(i,jj)                     
                     tused(i,jmin) = .true.
                     jfirst = jmin

                  else
                     jfirst = 0
c                     print*,'jused = ',jused,maxnumframes(i),time_maxa
                     
c                    Search for new starting point if                     
                     if (jused .lt. maxnumframes(i) .and.
     &                   time_maxa .gt. 30.d0) then
c                       If time_maxa < 30, there is guaranteed a valid point

c                       Insert spacer                        
                        jj = jj + 1
                        maxnumframesnew(i) = jj
                        xobv_red(i,jj) = -20.d0
                        yobv_red(i,jj) = -20.d0
                        zobv_red(i,jj) = -20.d0
                        tobv_red(i,jj) = -1.d0
                        if (i .eq. icell .and. prt) 
     &    print*,'tobv_red b(',i,jj,') = ',tobv_red(i,jj),xobv_red(i,jj)                        

                        j = jused + 1
                        do while (j .le. maxnumframes(i) .and.
     &                            jfirst .eq. 0)
                           if (.not. tused(i,j) .and.
     &                         tobv(i,j) .gt. -1.d0) then
                              if (tobv(i,j)-tobv(i,jused)-90.d0 .gt.
     &                            30.d0) then
                                 jfirst = j
                              end if
                           end if
                           j = j + 1
                        end do
c                        if (i .eq. icell .and. prt)
c     &                  print*,'jfirst = ',jfirst

                        if (jfirst .gt. 0) then
                           jj = jj + 1
                           jused = jfirst
                           jjnew = jjnew + 1
                           maxnumframesnew(i) = jj
                           xobv_red(i,jj) = xobv(i,jfirst)
                           yobv_red(i,jj) = yobv(i,jfirst)
                           zobv_red(i,jj) = zobv(i,jfirst)
                           tobv_red(i,jj) = tobv(i,jfirst)+tlast

                           if (i .eq. icell .and. prt)
     &    print*,'tobv_red c(',i,jj,') = ',tobv_red(i,jj),xobv_red(i,jj)

                           tused(i,jfirst) = .true.                           
                           tlast_new = tobv_red(i,jj)
                        else
                           print*,'ERROR STOP'
                           stop
                        end if
                     end if
c                    End search for new starting point if
                     
                  end if
c                 End Assigning reduced point if
                  
               end do
c              End do while jfirst = 0               
               

            end if

            if (i .eq. icell .and. prt)
     &      print*,'New ***********************',jjnew,jj
            if (tlast .lt. 1.d-8) then
               maxnumframes_save = maxnumframesnew(i)
            end if
         end do
c        Reuse jjnew loop
c         maxnumframes(i) = maxnumframes_save
         maxnumframes(i) = maxnumframesnew(i)
         
         maxnumframes_extend(i) = maxnumframesnew(i)
c            if (i .eq. icell) stop
      end do



c      goto 8891
      print*,'Reduction completed'
      do i = 1,number_tcells
c     do j = 1,number_frames
         do j = 1,maxnumframes(i)
            xobv(i,j) = xobv_red(i,j)
            yobv(i,j) = yobv_red(i,j)
            zobv(i,j) = zobv_red(i,j)
            tobv(i,j) = tobv_red(i,j)
         end do
      end do


      if (prt) then
      do i = 1,number_tcells
         do j = 1,maxnumframes(i)
            if (i .eq. icell .and. ii .eq. 1) then
               if (j .gt. 1) then
                  if (tobv(i,j-1) .lt. -.9d0) then
                     print*,'tobv(',i,j,') = ',tobv(i,j),
     &                       0.d0,xobv(i,j)
                  else
                     print*,'tobv(',i,j,') = ',tobv(i,j),
     &                    tobv(i,j)-tobv(i,j-1),xobv(i,j)
                  end if
               else
                  print*,'tobv(',i,j,') = ',tobv(i,j),
     &                    tobv(i,j),xobv(i,j)
               end if
            end if
         end do
      end do

      do i = 1,number_tcells
         do j = 1,maxnumframes_extend(i)
            if (i .eq. icell .and. ii .eq. 1) then
               if (j .gt. 1) then
                  if (tobv(i,j-1) .lt. -.9d0) then
                     print*,'tobve(',i,j,') = ',tobv(i,j),
     &                       0.d0,xobv(i,j)
                  else
                     print*,'tobve(',i,j,') = ',tobv(i,j),
     &                    tobv(i,j)-tobv(i,j-1),xobv(i,j)
                  end if
               else
                  print*,'tobve(',i,j,') = ',tobv(i,j),
     &                    tobv(i,j),xobv(i,j)
               end if
            end if
         end do
      end do
      end if

      
            
 8891 continue
c      print*,'xobv = ',xobv(18,100),xobv(18,101)
c      print*,'tobv = ',tobv(18,100),tobv(18,101)
      
      do i = 1,number_tcells
         xminz = 1.d+20
         yminz = 1.d+20
         zminz = 1.d+20

         xmaxz = -1.d+20
         ymaxz = -1.d+20
         zmaxz = -1.d+20
         do j = 1,maxnumframes(i)
c         do j = 1,number_frames
            if (tobv(i,j) .gt. -1.d0) then
c               xmin = min(xmin,xobv(i,j))                                                                             
c               ymin = min(ymin,yobv(i,j))                                                                             
c               zmin = min(zmin,zobv(i,j))                                                                             
c               xmax = max(xmax,xobv(i,j))                                                                             
c               ymax = max(ymax,yobv(i,j))                                                                             
c               zmax = max(zmax,zobv(i,j))                                                                             
               xminz = min(xminz,xobv(i,j))
               yminz = min(yminz,yobv(i,j))
               zminz = min(zminz,zobv(i,j))
               xmaxz = max(xmaxz,xobv(i,j))
               ymaxz = max(ymaxz,yobv(i,j))
               zmaxz = max(zmaxz,zobv(i,j))
            end if
         end do

c         if (i .eq. icell) then
c         print*,'xmin xmax = ',xminz,xmaxz
c         print*,'ymin ymax = ',yminz,ymaxz
c         print*,'zmin zmax = ',zminz,zmaxz
c         end if
         
c       do j = 1,number_frames
         do j = 1,maxnumframes(i)

            xobvc(i,j) = xobv(i,j)
            yobvc(i,j) = yobv(i,j)
            zobvc(i,j) = zobv(i,j)

            if (tobv(i,j) .gt. -1.d0) then
c               if (i .eq. icell) then
c                  print*,'i j xminz = ',i,j,xminz
c               end if
               
               xobv(i,j) = xobv(i,j) - xminz
               yobv(i,j) = yobv(i,j) - yminz
               zobv(i,j) = zobv(i,j) - zminz
               xmin = min(xmin,xobv(i,j))
               ymin = min(ymin,yobv(i,j))
               zmin = min(zmin,zobv(i,j))
               xmax = max(xmax,xobv(i,j))
               ymax = max(ymax,yobv(i,j))
               zmax = max(zmax,zobv(i,j))

               xsum = xsum + xobv(i,j)
               ysum = ysum + yobv(i,j)
               zsum = zsum + zobv(i,j)

               xsum2 = xsum2 + xobv(i,j)**2
               ysum2 = ysum2 + yobv(i,j)**2
               zsum2 = zsum2 + zobv(i,j)**2

               rnum = rnum + 1.d0

            end if
         end do

         rangexmax = xmaxz - xminz
         rangeymax = ymaxz - yminz
         rangezmax = zmaxz - zminz

         rangexm = max(rangexm,rangexmax)
         rangeym = max(rangeym,rangeymax)
         rangezm = max(rangezm,rangezmax)
      end do
c      print*,'xobv = ',xobv(18,100),xobv(18,101)
c      stop

      xavg = xsum/rnum
      yavg = ysum/rnum
      zavg = zsum/rnum

      xstd = sqrt( (xsum2 - xsum**2/rnum)/(rnum - 1.d0)  )
      ystd = sqrt( (ysum2 - ysum**2/rnum)/(rnum - 1.d0)  )
      zstd = sqrt( (zsum2 - zsum**2/rnum)/(rnum - 1.d0)  )

      xmin_box = xavg - xstd
      xmax_box = xavg + xstd
      ymin_box = yavg - ystd
      ymax_box = yavg + ystd
      zmin_box = zavg - xstd
      zmax_box = zavg + xstd


      rtc_count = 0.d0
      do i = 1,number_tcells
c     do j = 1,number_frames
         do j = 1,maxnumframes(i)
            if (tobv(i,j) .gt. -1.d0) then
               if ((xobv(i,j) .ge. xmin_box .and.
     &              xobv(i,j) .le. xmax_box) .and.
     &             (yobv(i,j) .ge. ymin_box .and.
     &              yobv(i,j) .le. ymax_box) .and.
     &             (zobv(i,j) .ge. zmin_box .and.
     &              zobv(i,j) .le. zmax_box)) then
                   rtc_count = rtc_count + 1.d0
               end if
            end if
         end do
      end do


      volume_box = (xmax_box - xmin_box)*
     &             (ymax_box - ymin_box)*
     &             (zmax_box - zmin_box)

      density_movie = rtc_count/volume_box

      rtc_count_all_movies = rtc_count_all_movies + rtc_count
      volume_all_movies = volume_all_movies + volume_box

      do i = 1,number_tcells
         if (kij .eq. 1) then
            total_tcell_count = total_tcell_count + 1
         end if
c         do j = 1,number_frames
c            if (tobv(i,j) .gt. -1.d-6) then
c               if (kij .eq. 1) then
c                  write(119,*) total_tcell_count,
c     &            xobv(i,j),yobv(i,j),zobv(i,j),tobv(i,j)
c               end if
c            end if
c         end do
      end do

      ij = 0

      pi = acos(-1.d0)

      if (kij .eq. 1) then

      iangle = 0

      numframestot = 0
      timesteptot = 0.d0
      do i = 1,number_tcells

         npersist(i) = 0
         j = 1
         icountangle = 0
         icountcell = 0

         jstart = 0

         if (maxnumframes_extend(i) .eq. 0) then
            print*,'ERROR No frames cell ',i
            stop
         end if
         
c         print*,'maxnumframes(',i,') = ',maxnumframes(i)
c     do while (j+2 .le. number_frames)
         do while (j+2 .le. maxnumframes(i))

            if (tobv(i,j) .gt. -1.d-6 .and.
     &          tobv(i,j+1) .gt. -1.d-6 .and.
     &          tobv(i,j+2) .gt. -1.d-6) then
               
               vec1x = xobv(i,j+1) - xobv(i,j)
               vec1y = yobv(i,j+1) - yobv(i,j)
               vec1z = zobv(i,j+1) - zobv(i,j)
               dist2 = vec1x**2 + vec1y**2 + vec1z**2 
               dist = sqrt(dist2)+1.d-16
               dista = dist
               timedif = tobv(i,j+1) - tobv(i,j)
               dispaa = dist/timedif
               timea = timedif
               velinitial = dist/timedif
               vec1x = vec1x/dist
               vec1y = vec1y/dist
               vec1z = vec1z/dist

               if (jstart .eq. 0) then
                  velsum = dist
                  timesum = timedif
                  jstart = 1
                  vec1xi = vec1x
                  vec1yi = vec1y
                  vec1zi = vec1z
               end if

               vec2x = xobv(i,j+2) - xobv(i,j+1)
               vec2y = yobv(i,j+2) - yobv(i,j+1)
               vec2z = zobv(i,j+2) - zobv(i,j+1)

               dist2 = vec2x**2 + vec2y**2 + vec2z**2 
               dist = sqrt(dist2)+1.d-16
               timedif = tobv(i,j+2)-tobv(i,j+1)
               timesum2 = timea + timedif
               distsum2 = dista + dist
               vel2 = dist
               vec2x = vec2x/dist
               vec2y = vec2y/dist
               vec2z = vec2z/dist
                  
               dot =
     &         min(max(vec1x*vec2x+vec1y*vec2y+vec1z*vec2z,-1.d0),1.d0)
               angle = acos(dot)*180.d0/pi

c               if (ii .eq. 9 .and. i .eq. 85) then
c                  print*,'xobv(',i,j,') = ',xobv(i,j),tobv(i,j)
c                  print*,'dist dista timea timedif = ',
c     &            dist,dista,timea,timedif
c               end if

               
c               if (dist .gt. 1.d-12 .and.
c     &             dista .gt. 1.d-12 .and.
c     &             timea .gt. 1.d-8 .and.
c     &              timedif .gt. 1.d-8) then
c               if (dist .gt. 1.d-12 .and.
c     &             dista .gt. 1.d-12 .and.
               if (timea .gt. 1.d-8 .and.
     &              timedif .gt. 1.d-8) then                  
                  if (60.d0*distsum2/timesum2 .gt. 1.d0) then
                     iangle = iangle + 1
                     turnang(iangle) = angle
                     if (angle .le. 90.d0) then
                        icountangle = icountangle + 1
                     end if
                     icountcell = icountcell + 1
                     velangle(iangle) = distsum2/timesum2
                     numangles(ii) = iangle
                  end if
               end if

               if (angle .lt. 90.d0 .and.
     &             dist .gt. 1.d-12 .and.
     &             dista .gt. 1.d-12) then
                  velsum = velsum + vel2
                  timesum = timesum + timedif
                  jstart = jstart + 1
                  vec1x = vec2x
                  vec1y = vec2y
                  vec1z = vec2z
                  dista = dist
               else
                  if (jstart .ge. 2) then
                     npersist(i) = npersist(i) + 1
                     velpersist(i,npersist(i)) = velsum/timesum
                     timepersist(i,npersist(i)) = timesum
                  end if
                  jstart = 0
               end if
            else
               if (jstart .ge. 2) then
                  npersist(i) = npersist(i) + 1
                  velpersist(i,npersist(i)) = velsum/timesum
                  timepersist(i,npersist(i)) = timesum
               end if
               jstart = 0
            end if

            j = j + 1

c        Frame loop            
         end do

         if (icountcell .gt. 5) then
            write(119,*) dble(icountangle)/dble(icountcell),icountcell,
     &      maxnumframes_extend(i)
         else
            if (maxnumframes_extend(i) .gt. 25) then
               print*,'maxframes cell = ',ii,i,maxnumframes_extend(i),
     &                icountcell,icountangle
            end if
         end if
         
         if (jstart .ge. 2) then
            npersist(i) = npersist(i) + 1
            velpersist(i,npersist(i)) = velsum/timesum
            timepersist(i,npersist(i)) = timesum
         end if

         maxbins = 15

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                 \                

c        Bayesian description    

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        do j = 2,number_frames-1
         do j = 2,maxnumframes(i)-1

            if (tobv(i,j-1) .gt. -1.d-6 .and.
     &          tobv(i,j) .gt. -1.d-6 .and.
     &          tobv(i,j+1) .gt. -1.d-6 .and.
     &         tobv(i,j)-tobv(i,j-1) .gt. 1.d-6 .and.
     &         tobv(i,j+1)-tobv(i,j) .gt. 1.d-6) then

               xdisp = xobv(i,j) - xobv(i,j-1)
               ydisp = yobv(i,j) - yobv(i,j-1)
               zdisp = zobv(i,j) - zobv(i,j-1)
               tdisp = tobv(i,j) - tobv(i,j-1) 
               disp2 = xdisp**2 + ydisp**2 + zdisp**2
               disp = sqrt(disp2)
c               if (i .eq. 18 .and. j .eq. 101) then
c                  print*,xobv(i,j-1),yobv(i,j-1),zobv(i,j-1),tobv(i,j-1)
c                  print*,xobv(i,j),yobv(i,j),zobv(i,j),tobv(i,j)
c                  print*,'disp = ',disp,tdisp
c               end if

c              print*,'disp = ',disp
c              Displacement velocity in um/min 

               dispvprev = 60.d0*disp/tdisp

               xdisp = xobv(i,j+1) - xobv(i,j)
               ydisp = yobv(i,j+1) - yobv(i,j)
               zdisp = zobv(i,j+1) - zobv(i,j)
               tdisp = tobv(i,j+1) - tobv(i,j) 
               disp2 = xdisp**2 + ydisp**2 + zdisp**2
               disp = sqrt(disp2)
               dispvafter = 60.d0*disp/tdisp

c              Convert speeds to integers                                                                                           
               ibin = dispvprev
               ibina = dispvafter
               if (ibin .gt. 100) then
                  print*,'ibin = ',ibin,dispvprev,disp,tdisp,ii,i,j
                  stop
               end if

               if (ibin .le. maxbins) then

c                 m in manuscript
                  total_num_frames = total_num_frames + 1.0

c                  do ibinc = 0,maxbins
c                     if (ibina .ge. ibinc-1 .and.
c     &                   ibina .le. ibinc+1) then
c                        probability_after(ibinc) =
c     &                  probability_after(ibinc) + 1.0
c                     end if
c     end do

c                 m_a in manuscript                   
                  probability_after(ibina) =
     &            probability_after(ibina) + 1.0

c                 m_b in manuscript                  
                  totalbin(ibin) =
     &            totalbin(ibin) + 1

c                 m_ab in manuscript                  
                  if (ibina .eq. ibin) then
                     persistbin(ibin) =
     &               persistbin(ibin) + 1
                  end if

               end if

c23456789012345678901234567890123456789012345678901234567890123456789012                                           \                
                                                                                                                                    
            end if
         end do


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                 \                
                                                                                                                                    
c        End: Bayesian description                                                                                 \               
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 


         nvel(i) = 0
         distsum = 0.d0
         timesum = 0.d0

         bound = 1.d0

c     do j = 1,number_frames-1
         do j = 1,maxnumframes(i)-1
            if (tobv(i,j) .gt. -1.d-6 .and.
     &          tobv(i,j+1) .gt. -1.d-6)  then
               
                  vec1x = xobv(i,j+1) - xobv(i,j)
                  vec1y = yobv(i,j+1) - yobv(i,j)
                  vec1z = zobv(i,j+1) - zobv(i,j)
                  dist2 = vec1x**2 + vec1y**2 + vec1z**2 
                  dist = sqrt(dist2)+1.d-12
                  timedif = tobv(i,j+1)-tobv(i,j)
                  distsum = distsum + dist
                  timesum = timesum + timedif
                  vec1x = vec1x/dist
                  vec1y = vec1y/dist
                  vec1z = vec1z/dist
                  if (timedif .gt. 1.d-9) then
                     write(889,*) timedif
                     numframestot = numframestot + 1
                     timesteptot = timesteptot + timedif
                     nvel(i)= nvel(i) + 1
                     vel(i,nvel(i)) = dist/timedif
                     velxa(i,nvel(i)) = vec1x
                     velya(i,nvel(i)) = vec1y
                     velza(i,nvel(i)) = vec1z
                     xpos(i,nvel(i)) = xobvc(i,j)
                     ypos(i,nvel(i)) = yobvc(i,j)
                     zpos(i,nvel(i)) = zobvc(i,j)
                     tpos(i,nvel(i)) = tobv(i,j)
                  end if

            end if
         end do
         if (timesum .gt. 1.d-12) then
            velavg(i) = distsum/timesum
         end if
         
      end do
      end if

      write(901,*) ii,nint(90.d0/(timesteptot/dble(numframestot)))
      
      if (kij .eq. 1 .and. ii .eq. mxnfit) then
c         do ibinc = 1,maxbins-1
         do ibinc = 0,maxbins            
c           prob_b bins will all sum to 1                                                                      
c           prob_a and prob_a_b bins will not all sum to 1                                                   

            prob_a = probability_after(ibinc)/
     &               total_num_frames
            prob_b = dble(totalbin(ibinc))/
     &               total_num_frames
            prob_a_b = dble(persistbin(ibinc))/
     &                 total_num_frames
            prob_a_given_b(ibinc) = prob_a_b/prob_b
            ratio_prob(ibinc) = prob_a_given_b(ibinc)/prob_a
            if (prob_b .gt. .02d0) then
            write(567,*) dble(ibinc)+.5d0,ratio_prob(ibinc),
     &      prob_a_given_b(ibinc),prob_a
            end if
         end do
      end if

      if (kij .eq. 2) then
c     Correlation study                                                                                                      
      do i = 1,number_tcells
         do j = 1,nvel(i)
            do iii = 1,number_tcells
               do jj = 1,nvel(iii)
                  sumi_tot = sumi_tot + 1
                  if (i .ne. iii) then
c     &            (abs(tpos(i,j)-tpos(iii,jj)) .gt. 5.d0*60.d0)) then
                        distx = xpos(i,j)-xpos(iii,jj)
                        disty = ypos(i,j)-ypos(iii,jj)
                        distz = zpos(i,j)-zpos(iii,jj)
                        dcell2 = distx**2+disty**2+distz**2
                        dcell = dsqrt(dcell2)
                        if (dcell .lt. 8.d0) then
                           sumi = sumi + 1
                           dot = velxa(i,j)*velxa(iii,jj) +
     &                           velya(i,j)*velya(iii,jj) +
     &                           velza(i,j)*velza(iii,jj)
                           write(131,*) dot
                        end if
                  end if
               end do
            end do
         end do
      end do
      if (ii .eq. mxnf) then
         print*,'Percentage of points = ',dble(sumi)/dble(sumi_tot),sumi
      end if

      end if

c      rewind(iunit)

      return
      end