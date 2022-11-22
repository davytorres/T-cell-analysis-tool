      program tcell
      implicit none

c     Maximum number of dendritic cells a tcell will visit
      integer max_unique_tcells_visit
c     Maximum number of dendritic cells in a distribution
      integer max_num_dendritic_cells
c     Maximum number of tcells in a single file
      integer max_tcells_allfiles
c     Maximum number of frames for a single tcell
      integer max_frames_allfiles
c     Maximum number of files
      integer mxnf
c     Grid dimension
      integer ncube
      
      parameter (max_unique_tcells_visit = 500)
      parameter (max_num_dendritic_cells = 10000)
      parameter (max_tcells_allfiles = 6000)
      parameter (max_frames_allfiles = 800)
      parameter (ncube = 168)
      parameter (mxnf = 40)

c     Grid used to compute volume traversed by T cell
      logical grid(ncube,ncube,ncube)
      integer mxnfit
      integer sumi,sumi_tot,idsum
      integer numbigiterations,numcf
      integer imap(0:mxnf)
      integer iig,jjg,kkg,ig,jg,kg,ijk,iii
      integer ii,i,j,k,kij,jj,ijka
      integer total_tcell_count
      integer number_tcells
      integer numbertcellstotal,npersistc
      integer num_angles90,num_anglesii

      
      integer jja,jlast_keep,numcellsz
      double precision timef_cons
      double precision botmsd,time_elapsed
      double precision distance_msd,distancez,distancef
      double precision speed_c,speed_d
      double precision standard_dev_y,standard_dev_z
      double precision time_msd,difft,ri90
      double precision diff_time,diff_time_min
      double precision vtcell,wsum
      
      double precision avg_slope,rn_slope
      double precision xpsum,ypsum,xpsum2,denomp
      double precision xypsum2,logv,rsp,slopep
      double precision slopefit,sdfit

      double precision eltime,distanceff
      double precision axcell,aycell,eps
      double precision xtcell,xtcell2
      double precision ytcell,xytcell,rtcell
      double precision xtcellavg,ytcellavg,yint

      integer totalbin(0:300),persistbin(0:300)
      double precision probability_after(0:300)
      double precision total_num_frames

      double precision rbitr
      double precision xrange,yrange,zrange
      double precision dxg,dxg2,vg
      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision diameter_tcell,dt
      double precision time_interval(max_tcells_allfiles,mxnf)
      double precision volumesum
c23456789012345678901234567890123456789012345678901234567890123456789012
      integer maxnumframes(max_tcells_allfiles)
      double precision xobv(max_tcells_allfiles,max_frames_allfiles)
      double precision yobv(max_tcells_allfiles,max_frames_allfiles)
      double precision zobv(max_tcells_allfiles,max_frames_allfiles)
      double precision tobv(max_tcells_allfiles,max_frames_allfiles)
      double precision vel(max_tcells_allfiles,max_frames_allfiles)
      double precision velavg(max_tcells_allfiles)
      double precision time_vol(max_tcells_allfiles)
      double precision sst,xxxi,yyyi,zzzi
      
      
      double precision 
     &          velpersist(max_tcells_allfiles,max_frames_allfiles)
      double precision 
     &          timepersist(max_tcells_allfiles,max_frames_allfiles)
      integer nvel(max_tcells_allfiles)
      integer npersist(max_tcells_allfiles)
      double precision volumetcell(max_tcells_allfiles)
      double precision tcellmaxtime(max_tcells_allfiles)
      double precision tcellmintime(max_tcells_allfiles)
      double precision number_frames_tcell(max_tcells_allfiles,mxnf)
      double precision turnang(max_tcells_allfiles*max_frames_allfiles)
      double precision velangle(max_tcells_allfiles*max_frames_allfiles)            

      double precision rcontact,runique
      double precision time_confined_total,time_tracked_total
      double precision sum_dist_traveled,dist_traveled_total
      
      integer number_framesa(mxnf)
      integer number_tcellsa(mxnf)
      integer numangles(mxnf)
      double precision xmina(mxnf),xmaxa(mxnf)
      double precision ymina(mxnf),ymaxa(mxnf)
      double precision zmina(mxnf),zmaxa(mxnf)
      double precision dta(mxnf)

      double precision volumefile(mxnf)

      integer num_constrainedf
      integer number_frames
      double precision distance2

      double precision xc,yc,zc
      double precision radius_tcell_squared

      double precision dist_traveled
      double precision dist_constrained
      double precision time_constrained_min
      double precision time_constrainedf

      integer jfirst,jjj
      double precision r,slope(max_tcells_allfiles)
      double precision xmsd(max_frames_allfiles)
      double precision ymsd(max_frames_allfiles)
      double precision sdx,sdy
      double precision rtc_count_all_movies
      double precision volume_all_movies
      double precision density_all_movies

      double precision sumx,sumxx,sumy,sumxy,rinum
      double precision slopeall

      integer i90,i90max
      integer i90c(1000)
      integer jfirsta(max_tcells_allfiles)
      integer jlasta(max_tcells_allfiles)
      logical continuetrack
      integer kstart,ja,ka
      integer icount(max_tcells_allfiles)
      integer jfirstaa(max_tcells_allfiles,1000)
      integer jlastaa(max_tcells_allfiles,1000)      
      double precision total_time(max_tcells_allfiles)
      double precision time_lapsej(1000)
      double precision weights_time(max_tcells_allfiles,1000)

      
      double precision x90(1000),y90(1000),z90(1000),y90_2(1000)
      double precision z90_2(1000),y90list(50,10000),y90lista(10000)
      double precision ym,y90alist(50,10000)
      integer i90ca(1000)
      double precision x90a(1000),y90a(1000),y90_2a(1000)      
      double precision slopeindv

      logical valid
      integer jk,num_valid,jki
      double precision meandering
      double precision meandering_ratio
      double precision sum_distance
      double precision time_lapse,speed_avg,speed_cell
      double precision time_min,speed_min,time_max
      double precision displacement_speed
      double precision speed_cell_sum,time_opt
      double precision time_slow,time_total
      double precision time_lapse_avg

      integer jlast
      double precision time_constrainedf_avg,time_constrainedf_sum      
      double precision time_confined,time_tracked
      double precision confined_ratio
      double precision confined_ratio_sum
      
      double precision y90asort(50,600)
      double precision y90cc(50,6000)

      open(unit=91,file='flu_unique')      
      open(unit=92,file='flu_angle')
      open(unit=93,file='flu_avg_velocity')      
      open(unit=94,file='flu_total')
      open(unit=95,file='flu_velocity')
      open(unit=97,file='flu_volume')
      open(unit=98,file='flu_persist')
      open(unit=99,file='flu_time')
      open(unit=127,file='flu_slopes')

      open(unit=130,file='elapsed_time')
      open(unit=765,file='flu_cons')
      open(unit=767,file='flu_cons_time')      
      
      open(unit=529,file='rms_slope')
      open(unit=531,file='rms_slopet')
      open(unit=766,file='displacement')
      open(unit=115,file='directionality_flu')
      open(unit=778,file='slopeall_rms')
      open(unit=779,file='cellspeed_vs_time')
      
c     File used for virtual simulations      
      open(unit=319,file='positions')
      open(unit=321,file='mxnf')

      eps = 1.d-14
      idsum = 0
      sumi = 0
      sumi_tot = 0

      rtc_count_all_movies = 0.d0
      volume_all_movies = 0.d0


      total_tcell_count = 0
      numbigiterations = 1
      
c     Diameter of a tcell in microns
      diameter_tcell = 10.d0
      radius_tcell_squared = (diameter_tcell/2.d0)**2

c     Minimum and maximum in x-, y-, and z- over all files
      do ii = 1,mxnf
          xmina(ii) = 0.d0;
          xmaxa(ii) = 0.d0;
          ymina(ii) = 0.d0;
          ymaxa(ii) = 0.d0;
          zmina(ii) = 0.d0;
          zmaxa(ii) = 0.d0;
          dta(ii) = 0.d0;
      end do

      totalbin = 0
      persistbin = 0
      probability_after = 0.d0
      total_num_frames = 0.d0
      
      do kij = 1,numbigiterations

      do i = 0,mxnf
         imap(i) = i
      end do

      i90max = -1
      do ijka = 1,1000
         i90c(ijka) = 0
         x90(ijka) = 0.d0
         y90(ijka) = 0.d0
         z90(ijka) = 0.d0
         y90_2(ijka) = 0.d0
         z90_2(ijka) = 0.d0
      end do
      y90list = 0.d0

      read(321,*) mxnfit,slopefit,sdfit
      write(319,*) mxnfit

      xtcell = 0.d0
      xtcell2 = 0.d0
      ytcell = 0.d0
      xytcell = 0.d0
      rtcell = 0.d0
      
      do jj = 0,mxnfit

         ii = imap(jj)

         if (ii .gt. 0) then
c23456789012345678901234567890123456789012345678901234567890123456789012   
            call observedtcell(kij,ii,max_tcells_allfiles,
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

            if (kij .eq. 1) then

               write(319,*) number_tcells,number_frames
               do i = 1,number_tcells
                  do j = 1,number_frames
                     write(319,*)
     &               xobv(i,j),yobv(i,j),zobv(i,j),tobv(i,j)
                  end do
               end do

c              observedtcell translates are trajectories so the minimum
c              for all x-,y-, and z-coordinates began at (0.,0.,0.)
               xmina(ii) = xmin
               xmaxa(ii) = xmax
               ymina(ii) = ymin
               ymaxa(ii) = ymax
               zmina(ii) = zmin
               zmaxa(ii) = zmax
               dta(ii) = dt

               xrange = xmaxa(ii) - xmina(ii)
               yrange = ymaxa(ii) - ymina(ii)
               zrange = zmaxa(ii) - zmina(ii)
c              We assume all T cells in the file lie within the
c              range (0:400um, 0:400um, 0:400um)

c               print*,'xrange yrange zrange = ',
c     &         xrange,yrange,zrange
c               stop

               if (xrange .gt. 400.d0 .or.
     &             yrange .gt. 400.d0 .or.
     &             zrange .gt. 400.d0) then
                  print*,'file = ',ii
                  print*,'xrange yrange zrange = ',
     &                    xrange,yrange,zrange
                  stop
               else
c                  print*,'xrange yrange zrange = ',
c     &                    xrange,yrange,zrange
               end if
                  
               
               number_framesa(ii) = number_frames
               number_tcellsa(ii) = number_tcells
               
c               print*,'ii number_tcells = ',ii,number_tcells

c     We assume our grid ranges between
c              (-10:410,-10:410,-10:410)
               
               dxg = 420.d0/dble(ncube)
c              dxg should be 2.5 microns
               dxg2 = dxg/2.d0
c              Volume of 2.5um x 2.5um x 2.5um cube               
               vg = dxg**3


               time_min = 60.d0*9.d0
               time_opt = 60.d0*10.5d0
               time_max = 60.d0*12.d0

               time_min = 0.d0
               time_opt = 1.d+6
               time_max = 1.d+6
c               time_max = 60.d0*12.d0
               numcellsz = 0
c               print*,'here',ii
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
c              Identify continuous tracks within the T cell frames               
               do i = 1,number_tcells

                  icount(i) = 0
                  kstart = 1
                  continuetrack = .true.

c                  do j = 1,maxnumframes(i)
c                     print*,j,xobv(i,j),tobv(i,j)
c                  end do

                  do while (continuetrack)
                  
                     jfirst = 0
                     k = kstart
                     do while (jfirst .eq. 0 .and.
     &                   k .lt. maxnumframes(i))
c                        print*,'tobv(',i,k,') = ',tobv(i,k)
                        if (tobv(i,k) .gt. -1.d-6) then
                           jfirst = k
                        end if
                        k = k + 1
                     end do
c                     print*,'jfirst = ',jfirst

                     jlast = 0
                     if (jfirst .gt. 0) then
                        k = jfirst+1
                        do while (jlast .eq. 0 .and.
     &                       k .le. maxnumframes(i))
c                           print*,'tobv l(',i,k,') = ',tobv(i,k)
                           if (tobv(i,k) .le. -1.d-6) then
                              jlast = k-1
                           else
                              k = k + 1
                           end if
                        end do
                     end if

c                     print*,'k = ',k,jfirst,jlast,maxnumframes(i)

                     if (jfirst .gt. 0) then
c                        print*,'jfirst = ',jfirst,k,maxnumframes(i)
                        if (k .gt. maxnumframes(i)) then
c                  print*,tobv(i,maxnumframes(i))-tobv(i,jfirst),time_min
                           if (tobv(i,maxnumframes(i)) - tobv(i,jfirst)
     &                          .gt. time_min) then
c                  print*,tobv(i,maxnumframes(i))-tobv(i,jfirst),time_min
                              diff_time_min = 1.d+30
c                              print*,'jfirst = ',jfirst
                              do jki = jfirst+1, maxnumframes(i)
                                 time_lapse =
     &                                tobv(i,jki) - tobv(i,jfirst)
                                 diff_time =
     &                                abs(time_opt - time_lapse)
c                            print*,'time_lapse = ',time_lapse,diff_time
                                 if (diff_time .lt. diff_time_min)
     &                           then
                                    diff_time_min = diff_time
                                    jlast_keep = jki
                                 end if
                              end do
c                              print*,'jlast_keep = ',jlast_keep

                              if (tobv(i,jlast_keep) - tobv(i,jfirst)
     &                            .lt. time_max) then
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
c                           print*,'jfirstaa(',i,icount(i),') = ',
c     &                     jfirstaa(i,icount(i)),jlastaa(i,icount(i)),
c     &                             tobv(i,jlast_keep)-tobv(i,jfirst)
c                           stop
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
                              if (tobv(i,jlast) - tobv(i,jfirst)
     &                             .gt. time_min) then

                                 diff_time_min = 1.d+20
                                 do jki = jfirst+1, jlast
                                 time_lapse =
     &                                tobv(i,jki) - tobv(i,jfirst)                                    
c                                    time_lapse =
c     &                              tobv(i,jki) - tobv(i,jfirst)
c                          print*,'time_lapse = ',jki,time_lapse,time_opt                                    
                                    diff_time =
     &                              abs(time_opt - time_lapse)
                                    if (diff_time .lt. diff_time_min)
     &                              then
                                       diff_time_min = diff_time
                                       jlast_keep = jki
                                    end if
                                 end do
c                              print*,'jlast_keep = ',jlast_keep
c                              stop
                                 
                                 if (tobv(i,jlast_keep) - tobv(i,jfirst)
     &                               .lt. time_max) then                                 
                                    icount(i) = icount(i) + 1
                                 if (icount(i) .gt. 1000) then
                                 print*,'Increase 2nd arg of jfirstaa'
                                 stop
                                 end if
                                    jfirstaa(i,icount(i)) = jfirst
                                    jlastaa(i,icount(i)) = jlast_keep
c                             print*,'jfirstaa(',i,icount(i),') = ',
c     &                       jfirstaa(i,icount(i)),jlastaa(i,icount(i)),
c     &                       tobv(i,jlast_keep)-tobv(i,jfirst)
                                    kstart = jlastaa(i,icount(i)) + 1
                                 else
                                print*,'Error in continuous tracks B'
                             print*,tobv(i,jlast_keep) - tobv(i,jfirst)
                                print*,time_max
                                    stop
                                 end if
                              else
                                 kstart = jlast+1
                              end if
                           end if
                        end if

                     else
                        continuetrack = .false.
                     end if

                  end do
c                  if (i .eq. 1) stop
                  
               end do


               do i = 1,number_tcells
                  total_time(i) = 0.d0
                  do jki = 1,icount(i)
                     ja = jfirstaa(i,jki)
                     ka = jlastaa(i,jki)
                     time_lapsej(jki) = tobv(i,ka)-tobv(i,ja)
                     if (time_lapsej(jki) .lt. 150.d0) then
                        time_lapsej(jki) = 0.d0
                     end if
                     total_time(i) = total_time(i) + time_lapsej(jki)
                  end do
                  if (total_time(i) .gt. 0.d0) then
                     do jki = 1,icount(i)
                        weights_time(i,jki) =
     &                  time_lapsej(jki)/total_time(i)
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
c                  if (dabs(wsum-1.d0) .gt. 1.d-8) then
c                     print*,'Error wsum',i,wsum,total_time(i)
c                  end if
c                  if (total_time(i) .lt. 1.d-8) then
c                     print*,'Error total_time = ',i
c                     stop
c                  end if
               end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              Calculate the volume traversed by each T cell
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
c                  end if


                  if (icount(i) .gt. 0) numcellsz = numcellsz + 1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                  
c                 Mean squared displacement calculation                  

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
                  
               end do
               print*,'numcellsz = ',numcellsz,number_tcells

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               

*******************************************************************************************        
                                                                            
c              T cell cannot move more than 10 um from position to be considered in constrained mode                                                          
                                                                                       
               dist_constrained = 5.d0
               time_constrained_min  = 5.d0*60.d0
c              time_constrained_min is 5 minutes                                       
                                                                                        
c*******************************************************************************************         

c*******************************************************************************************                        
c*************************** Compute cell-based speed, displacement speed, and meandering ratio
c*******************************************************************************************
               
               speed_min = 1.d0/60.d0

c               time_min = 0.d0
c               time_max = 1.d+20
               do i = 1,number_tcells
                  displacement_speed = 0.d0
                  meandering = 0.d0
                  speed_cell_sum = 0.d0
                  num_valid = 0
                  time_slow = 0.d0
                  time_total = 0.d0
                  time_lapse_avg = 0.d0

                  sum_dist_traveled = 0.d0
                  dist_traveled_total = 0.d0
                  
                  do jki = 1,icount(i)

                     j = jfirstaa(i,jki)
                     k = jlastaa(i,jki)
                        
                     valid = .true.
                     do jk = j,k
                        if (tobv(i,jk) .le. -1.d0) then
                           valid = .false.
                        end if
                     end do

                     if (valid) then

                        time_lapse = tobv(i,k)-tobv(i,j)
                        if (time_lapse .gt. time_min .and.
     &                      time_lapse .lt. time_max) then

                           sum_distance = 0.d0
                           do jk = j+1,k
                              distance2 =
     &                           (xobv(i,jk)-xobv(i,jk-1))**2 
     &                         + (yobv(i,jk)-yobv(i,jk-1))**2 
     &                         + (zobv(i,jk)-zobv(i,jk-1))**2
                              dist_traveled = sqrt(distance2)
                              sum_distance = sum_distance +
     &                                       dist_traveled
c                                 write(99,*) tobv(i,jk)-tobv(i,jk-1)
                           end do

                           speed_cell = sum_distance/time_lapse
c                          Straight line distance                              
                           distance2 =(xobv(i,k)-xobv(i,j))**2 +
     &                                (yobv(i,k)-yobv(i,j))**2 +
     &                                (zobv(i,k)-zobv(i,j))**2
                           dist_traveled = sqrt(distance2)
                           speed_avg = dist_traveled/time_lapse
c                           write(130,*) time_lapse/60.d0,
c     &                                  60.d0*dist_traveled/time_lapse
                           if (speed_avg .lt. speed_min) then
                              time_slow = time_slow + time_lapse
                           end if
                           time_total = time_total + time_lapse

                           num_valid = num_valid + 1

                           time_lapse_avg = time_lapse_avg + time_lapse
                           
                           meandering_ratio =
     &                     dist_traveled/(sum_distance+1.d-15)

                           sum_dist_traveled = sum_dist_traveled +
     &                     sum_distance
                           dist_traveled_total = dist_traveled_total +
     &                     dist_traveled
                           
                           displacement_speed = displacement_speed +
     &                     speed_avg*weights_time(i,jki)

                           meandering = meandering +
     &                     meandering_ratio*weights_time(i,jki)

                           speed_cell_sum = speed_cell_sum +
     &                     speed_cell*weights_time(i,jki)
                           
                        end if
                           
                     else
                        print*,'Error in frame construction'
                        stop
                     end if
                        
                  end do

                  if (num_valid .ge. 1 .and.
     &                total_time(i) .gt. 1.d-3) then
                     write(130,*) (time_lapse_avg/60.d0)/dble(num_valid)
                     write(766,*) 60.d0*displacement_speed  
c     &               60.d0*displacement_speed/dble(num_valid)
c                     write(115,*) meandering/dble(num_valid)
                     write(115,*) dist_traveled_total/sum_dist_traveled                     
c                     write(93,*) speed_cell_sum*60.d0/dble(num_valid)
                     write(93,*) speed_cell_sum*60.d0                     
                  end if

               end do

*******************************************************************************************
********************* Computing confined ratio  **********************************************
*******************************************************************************************
               
               do i = 1,number_tcells

                  confined_ratio_sum = 0.d0
                  time_confined_total = 0.d0
                  time_tracked_total = 0.d0
                  
                  do jki = 1,icount(i)

                     jfirst = jfirstaa(i,jki)
                     jlast = jlastaa(i,jki)

                     time_confined = 0.d0
                     time_tracked = 0.d0

                     time_lapse = tobv(i,jlast)-tobv(i,jfirst)
                     if (time_lapse .gt. time_min .and.
     &                   time_lapse .lt. time_max) then
                     
                        j = jfirst
                        do while (j .lt. jlast)

                           time_lapse = 0.d0
                           dist_traveled = 0.d0

                           k = j
                           do while (dist_traveled .lt. 5.0 .and.
     &                               k .lt. jlast)
                              k = k + 1
c23456789012345678901234567890123456789012345678901234567890123456789012                           
                              distance2 = (xobv(i,k)-xobv(i,j))**2 +
     &                                    (yobv(i,k)-yobv(i,j))**2 +
     &                                    (zobv(i,k)-zobv(i,j))**2                           
                              dist_traveled = sqrt(distance2)
c                        print*,'dist_traveled = ',
c     &                  j,dist_traveled,tobv(i,k),k
                           end do
                           difft = tobv(i,k)-tobv(i,j)
                           time_tracked = time_tracked + difft
c                     print*,'difft time_tracked = ',difft,time_tracked
                           if (difft .gt. 150.d0) then
c                              print*,'time_confined = ',difft
                              time_confined = time_confined + difft
                           end if
                           j = k
                      
                        end do
                        
                     end if
                     time_confined_total = time_confined_total +
     &               time_confined
                     time_tracked_total = time_tracked_total +
     &               time_tracked
                     confined_ratio = time_confined/time_tracked
                     confined_ratio_sum = confined_ratio_sum +
     &               confined_ratio
                  
                  end do

                  if (icount(i) .gt. 0 .and.
     &                total_time(i) .gt. 1.d-3) then
c                    write(765,*) confined_ratio_sum/dble(icount(i))
                     if (time_tracked_total .gt. 1.d-12) then
                     write(765,*) time_confined_total/time_tracked_total
                     end if
                  end if
                  
               end do
               

*******************************************************************************************
********************* End: Computing ratio  **********************************************
*******************************************************************************************

c23456789012345678901234567890123456789012345678901234567890123456789012                  

*******************************************************************************************
********************* Computing confined time **********************************************
*******************************************************************************************

               do i = 1,number_tcells

                  time_constrainedf_sum = 0.d0
                  num_constrainedf = 0
                  
                  do jki = 1,icount(i)

                     jfirst = jfirstaa(i,jki)
                     jlast = jlastaa(i,jki)

                     time_constrainedf = 0.d0
                     time_constrainedf_avg = 0.d0
                     numcf = 0

                     time_tracked = tobv(i,jlast) - tobv(i,jfirst)
c                    time_constrained_min is 5 minutes

c                     print*,'time_tracked = ',
c     &               time_tracked,time_constrained_min
                     if (time_tracked .gt. time_min .and.
     &                   time_tracked .lt. time_max) then

c                        do j = 1,number_frames
                        do j = jfirst,jlast-1
                           dist_traveled = 0.0
                           k = j
                           timef_cons = 0.d0
c                          dist_constrained is 5 microns                           
                           do while ((dist_traveled .lt.
     &                                dist_constrained) .and.
     &                                k .lt. jlast)
                              k = k + 1
                              distance2 =(xobv(i,k)-xobv(i,j))**2 +
     &                                   (yobv(i,k)-yobv(i,j))**2 +
     &                                   (zobv(i,k)-zobv(i,j))**2
c                             Distance moved by T-cell between frames                                           
                                                                                                                     
                              dist_traveled = sqrt(distance2)
c                              print*,'j k = ',j,k
c                              print*,'dist_traveled dist_cons = ',
c     &                        dist_traveled,dist_constrained,tobv(i,k)
                              if (dist_traveled .ge.
     &                            dist_constrained) then
                                 timef_cons = tobv(i,k)-tobv(i,j)
                                 numcf = numcf + 1
c                                 print*,'timef_cons = ',timef_cons,numcf
                                 time_constrainedf =
     &                           time_constrainedf + timef_cons         
c                                 print*,'time_constrainedf = ',
c     &                         timef_cons,time_constrainedf,numcf
                                 dist_traveled = 2.0*dist_constrained
                              end if
                           end do
                        end do
                     end if

                     if (numcf .gt. 0) then
                        num_constrainedf = num_constrainedf + 1
                        time_constrainedf_avg =
     &                  (time_constrainedf/60.d0)/dble(numcf)

                        time_constrainedf_sum = time_constrainedf_sum +
     &                  time_constrainedf_avg*weights_time(i,jki)
                     end if
                     
                  end do
                     

                  if (num_constrainedf .gt. 0 .and.
     &                total_time(i) .gt. 1.d-3) then
c                     write(767,*)
c     &                    time_constrainedf_sum/dble(num_constrainedf)
                     write(767,*)
     &               time_constrainedf_sum
                  end if                     

               end do
*******************************************************************************************
********************* End computing confined time **********************************************
*******************************************************************************************               


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               num_anglesii = 0
               num_angles90 = 0
               do i = 1,numangles(ii)
                  if (velangle(i)*60.d0 .le. 30.d0) then
                     num_anglesii = num_anglesii + 1
                     write(92,*) turnang(i),velangle(i)*60.d0
                     if (turnang(i) .le. 90.d0) then
                        num_angles90 = num_angles90 + 1
                     end if
                  end if
               end do
c               write(119,*) dble(num_angles90)/dble(num_anglesii)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              Calculate total time each T cell is tracked
               volumesum = 0.d0
               do i = 1,number_tcells
                  tcellmaxtime(i) = -1.d0
                  tcellmintime(i) = 1.d+10
                  number_frames_tcell(i,ii) = 0
                  time_interval(i,ii) = 0.d0
                  do j = 1,number_frames
                     if (tobv(i,j) .gt. -1.d-6) then
c23456789012345678901234567890123456789012345678901234567890123456789012   
                        number_frames_tcell(i,ii) = 
     &                  number_frames_tcell(i,ii) + 1
                        tcellmaxtime(i) = max(tcellmaxtime(i),tobv(i,j))
                        tcellmintime(i) = min(tcellmintime(i),tobv(i,j))
                     end if
                  end do
                  if (tcellmaxtime(i) .lt. tcellmintime(i)) then
                     print*,'tcellmaxtime tcellmintime = ',
     &               tcellmaxtime(i),tcellmintime(i)
                     do j = 1,number_frames
                         print*,'tobv(',i,j,') = ',tobv(i,j)
                     end do
                     stop
                  end if
                  time_interval(i,ii) = tcellmaxtime(i)-tcellmintime(i)

                  if (time_interval(i,ii) .lt. 1.d-12) then
                     time_interval(i,ii) = -1.d0
                     velavg(i) = -1.e+10
                  end if
c                  if (time_vol(i) .gt. 0.d0) then
c                     volumetcell(i) = volumetcell(i)/time_vol(i)
c                  end if
                  volumesum = volumesum + volumetcell(i)
               end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               volumefile(ii) = volumesum/dble(numcellsz)

               do i = 1,number_tcells
                  if (time_interval(i,ii) .gt. 0.d0) then
                     do j = 1,nvel(i)
                        write(95,*) vel(i,j)*60.d0
                     end do
                     if (volumetcell(i) .gt. 0.d0 .and.
     &                   total_time(i) .gt. 1.d-3) then
                        write(97,*) volumetcell(i)
                     end if
                  end if
               end do
                  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               

c           if condition for kij = 1
            end if


c        if condition for ii .ne. value
         end if


c     Loop over files ii
      end do


      sumx = 0.d0
      sumxx = 0.d0
      sumy = 0.d0
      sumxy = 0.d0
      rinum = 0

      do ii = 1,5
c      do ii = 1,60
         if (i90c(ii) .gt. 0) then
            x90(ii) = x90(ii)/dble(i90c(ii))
            y90(ii) = y90(ii)/dble(i90c(ii))

            do jj = 1,i90c(ii)
               y90lista(jj) = y90list(ii,jj)
            end do
            call median(y90lista,ym,i90c(ii),10000)

c            y90(ii) = ym
            
c            if (x90(ii) .gt. 2.5d0) then
c               ri90 = 1.d0/dble(i90c(ii))
c               standard_dev_y = 
c     &         sqrt((y90_2(ii) - ri90*y90(ii)**2)/dble(i90c(ii)-1))
c               standard_dev_z =
c     &         sqrt((z90_2(ii) - ri90*z90(ii)**2)/dble(i90c(ii)-1))            
cc             x90(ii) = x90(ii)/dble(i90c(ii))
c
c               z90(ii) = z90(ii)/dble(i90c(ii))            
               sumx = sumx + log10(x90(ii))
               sumy = sumy + log10(y90(ii))
               sumxx = sumxx + log10(x90(ii))**2
               sumxy = sumxy + log10(x90(ii))*log10(y90(ii))
               rinum = rinum + 1.d0
               write(881,*) x90(ii),y90(ii)
               write(882,*) x90(ii),ym
c               write(778,*) x90(ii),i90c(ii),standard_dev_y               
c               write(779,*) x90(ii),z90(ii),standard_dev_z
c            end if
         end if
      end do
      
      rinum = 1.d0/rinum
      slopeall = (sumxy - sumx*sumy*rinum)/(sumxx - sumx*sumx*rinum)
      print*,'slopeall = ',slopeall
      write(780,*) slopeall

      density_all_movies = rtc_count_all_movies/volume_all_movies
      print*,'density_all_movies = ',density_all_movies

c     Loop over big iterations
      end do
      
      rinum = 1.d0/rtcell
      botmsd = xtcell2 - xtcell*xtcell*rinum
      if (abs(botmsd) .gt. 1.d-6) then
         slopeindv = (xytcell - xtcell*ytcell*rinum)/
     &        botmsd
         ytcellavg = ytcell*rinum
         xtcellavg = xtcell*rinum
         yint = ytcellavg - slopeindv*xtcellavg
         write(531,*) slopeindv
         write(532,*) 0.d0,yint
         write(532,*) 1.1d0,slopeindv*1.1d0 + yint
      end if


      

      numbertcellstotal = 0

      do jj = 1,mxnf
         ii = imap(jj)
         numbertcellstotal = numbertcellstotal +
     &   number_tcellsa(ii)
      end do

      print*,'numbertcellstotal flu  = ',numbertcellstotal
      print*,'Finished'

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
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
         
            
