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

               call identify_tracks(ii, number_tcells, icount,
     &                             max_tcells_allfiles, 
     &                             max_frames_allfiles,
     &                             maxnumframes, tobv, time_min, time_max,
     &                             time_opt, jfirstaa, jlastaa, total_time,
     &                             weights_time)
 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

c              Loop through every ith cell and initialize variables to be 0
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

c                 Loop through the sets of frames for cell i based on
c                 how many times it disappeared from view                   
                  do jki = 1,icount(i)

c                    First and last frames of set jki
                     j = jfirstaa(i,jki)
                     k = jlastaa(i,jki)
                        
                     valid = .true.

c                    Loop through the frames of set jki
c                    If a frame has a time stamp less than 0, it is ignored
                     do jk = j,k
                        if (tobv(i,jk) .le. -1.d0) then
                           valid = .false.
                        end if
                     end do

                     if (valid) then

c                       If time difference between the first and last
c                       frames of set jki are in the range of allowed
c                       time, sum up the distance traveled in set jki  
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

c                          Calclate frame-based speed of cell i in set jki
c                          Then calculate average speed of cell i in jki
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
                           
c                          Calculate total elapsed time in set jki
                           time_total = time_total + time_lapse

c                          Keep count of valid frames
                           num_valid = num_valid + 1

                           time_lapse_avg = time_lapse_avg + time_lapse
                           
c                          Measure the tendency for cell i to deviate
c                          from a straight line path in set jki
                           meandering_ratio =
     &                     dist_traveled/(sum_distance+1.d-15)

c                          Sum up distace traveled between frames and
c                          the distance between the first and last frame
                           sum_dist_traveled = sum_dist_traveled +
     &                     sum_distance
     
c                          Sum up distances between sets    
                           dist_traveled_total = dist_traveled_total +
     &                     dist_traveled
                           
c                          Sum up displacement speed with a weighted
c                          average speed based on how long cell i was
c                          presenent in set jki
                           displacement_speed = displacement_speed +
     &                     speed_avg*weights_time(i,jki)

c                          Sum up meandering ratio in sets where sets
c                          with more elapsed time have a higher weight
                           meandering = meandering +
     &                     meandering_ratio*weights_time(i,jki)

c                          Sum up speed of cell i by summing weighted
c                          speeds
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
         
            
