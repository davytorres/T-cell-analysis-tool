      subroutine readtcell(ii,max_tcells_allfiles,
     &                     max_frames_allfiles,
     &                     xobv,yobv,zobv,tobv,idmax_movie,
     &                     number_tcells,number_frames,idsum,
     &                     maxnumframes)

      implicit none

      character(1) :: c1
      character(2) :: c2
      character(100) :: filenamecsv
      logical idused
      integer max_frames_allfiles,max_tcells_allfiles
      integer iunit,nentries,maxlinesread,idss,ids,id,ifr
      integer number_tcells,number_frames,idsum,iia
      integer maxnumtcells
      integer idmax_movie,i,j,kk      
      integer iiframe,ii
      parameter (maxlinesread=20000)
      integer maxnumframes(max_tcells_allfiles)
      integer numframes(max_tcells_allfiles)      
      integer idcount(max_tcells_allfiles)
      integer map(maxlinesread)
      integer inversemap(max_tcells_allfiles)
      integer iframe(maxlinesread)
      integer tcelldatai(maxlinesread)
      double precision tcelldata(maxlinesread,5)      
      double precision tcelldata5(maxlinesread)

c23456789012345678901234567890123456789012345678901234567890123456789012
      double precision xobv(max_tcells_allfiles,max_frames_allfiles)
      double precision yobv(max_tcells_allfiles,max_frames_allfiles)
      double precision zobv(max_tcells_allfiles,max_frames_allfiles)
      double precision tobv(max_tcells_allfiles,max_frames_allfiles)

      open(unit=419,file='positions_new')
      
      numframes = 0


      iia = ii + 1

      if (iia .lt. 10) then
         write(c1,10) iia
 10      format(I1)
         filenamecsv = 'file'//c1//'.txt'
      else
         write(c2,20) iia
 20            format(I2)
         filenamecsv = 'file'//c2//'.txt'
      end if

      iunit = 10+iia
      open(unit = iunit,file=filenamecsv)
c      read(iunit,*) nentries
      print*,'iunit = ',iunit
      nentries = 0
      do
      read(iunit,*,END=15)
         nentries = nentries + 1
      end do
 15   rewind(iunit)

      print*,'filenamecsv = ',
     &       filenamecsv,nentries

      do i = 1,max_tcells_allfiles
         idcount(i) = 0
      end do
      number_frames = -1

      number_tcells = 0
      map = 0
      inversemap = 0

      print*,'nentries = ',nentries
      if (nentries .gt. maxlinesread) then
         print*,'nentries = ',nentries
         print*,'maxlinesread = ',maxlinesread
         stop
      end if

      do i = 1,nentries
         read(iunit,*) tcelldatai(i),
     &                 tcelldata(i,2),
     &                 tcelldata(i,3),tcelldata(i,4),
     &                 idss,tcelldata5(i)
         tcelldatai(i) = tcelldatai(i)+1
         numframes(tcelldatai(i)) = 
     &   numframes(tcelldatai(i)) + 1
         iiframe = numframes(tcelldatai(i))
         if (iiframe .gt. max_frames_allfiles) then
            print*,'iiframe = ',iiframe
            print*,'max_frames_allfiles = ',max_frames_allfiles
            stop
         end if

c        This is done so id's are consecutive and
c        number_tcells count is accurate
         ids = tcelldatai(i)
         idused = .false.
         if (number_tcells .gt. 0) then
c           Go through list of names assigned
            do kk = 1,number_tcells
               if (inversemap(kk) .eq. ids) then
                  idused = .true.
               end if
            end do
         end if
         if (.not. idused) then
            number_tcells = number_tcells + 1
            if (ids .gt. maxlinesread .or.
     &          number_tcells .gt. max_tcells_allfiles) then
               print*,'number_tcells = ',number_tcells
               print*,'ids = ',ids
               print*,'Increase dimension of map'
               stop
            end if
            map(ids) = number_tcells
            inversemap(number_tcells) = ids
         end if

         iframe(i) = iiframe
         number_frames = max(number_frames,iiframe)
      end do
      
      do i = 1,max_tcells_allfiles
         do j = 1,max_frames_allfiles
            xobv(i,j) = -20.d0
            yobv(i,j) = -20.d0
            zobv(i,j) = -20.d0
            tobv(i,j) = -1.d0
        end do
      end do

      maxnumtcells = 0
      maxnumframes = 0
      idmax_movie = -100000
      do i = 1,nentries
         id = map(tcelldatai(i))
         ifr = iframe(i)
         if (id .gt. max_tcells_allfiles) then
            print*,'Increase max_tcells_allfiles'
            stop
         end if
         if (ifr .gt. max_frames_allfiles) then
            print*,'Increase max_frames_allfiles'
            stop
         end if

         maxnumtcells = max(maxnumtcells,id)
         maxnumframes(id) = max(maxnumframes(id),ifr)
         xobv(id,ifr) = tcelldata(i,2)
         yobv(id,ifr) = tcelldata(i,3)
         zobv(id,ifr) = tcelldata(i,4)
         tobv(id,ifr) = tcelldata5(i)
c         if (id .eq. 2) then
c            print*,'tobv(',id,ifr,') = ',tobv(id,ifr)
c         end if

         if (tobv(id,ifr) .gt. -1.d0) then
            write(419,*)
     &      id+idsum,xobv(id,ifr),yobv(id,ifr),zobv(id,ifr),tobv(id,ifr)
            idmax_movie = max(idmax_movie,id)
         end if

      end do

      end
