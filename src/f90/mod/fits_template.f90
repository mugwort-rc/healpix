!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
! template for routine SP/DP overloading for module fitstools
!
!    subroutine input_map_KLOAD
!    subroutine read_bintab_KLOAD
!    subroutine read_conbintab_KLOAD
!    subroutine write_bintab_KLOAD
!    subroutine write_asctab_KLOAD
!    subroutine dump_alms_KLOAD
!    subroutine write_alms_KLOAD
!    subroutine read_alms_KLOAD
!    subroutine read_bintod_KLOAD
!    subroutine write_bintabh_KLOAD
!
!
! K M A P   : map kind                 either SP or DP
!
! K L O A D : suffixe of routine name, to be replaced by either s (SP) or d (DP)
  !=======================================================================
  subroutine input_map_KLOAD(filename, map, npixtot, nmaps, fmissval, header, units, extno)
    !=======================================================================
    !     reads fits file
    !     filename = fits file (input)
    !     map      = data read from the file (ouput) = real*4 array of size (npixtot,nmaps)
    !     npixtot  = length of the map (input)
    !     nmaps     = number of maps
    !     fmissval  = OPTIONAL argument (input) with the value to be given to missing
    !             pixels, its default value is 0
    !     header    = OPTIONAL (output) contains extension header
    !     units     = OPTIONAL (output) contains the maps units
    !     extno     = OPTIONAL (input)  contains the unit number to read from (0 based)
    !
    !     modified Feb 03 for units argument to run with Sun compiler
    !=======================================================================

    INTEGER(I4B),     INTENT(IN)           :: npixtot, nmaps
    REAL(KMAP),       INTENT(OUT), dimension(0:,1:) :: map
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    REAL(KMAP),       INTENT(IN), OPTIONAL :: fmissval
    CHARACTER(LEN=*), INTENT(OUT), dimension(1:), OPTIONAL :: header
    CHARACTER(LEN=*), INTENT(OUT), dimension(1:), OPTIONAL :: units
    INTEGER(I4B)                   , optional,  INTENT(IN) :: extno

    INTEGER(I4B) :: i,imap
    REAL(KMAP)     :: fmissing, fmiss_effct
    INTEGER(I4B) :: imissing
    integer(I4B) :: obs_npix, nlheader

    LOGICAL(LGT) :: anynull
    integer(I4B), dimension(:), allocatable :: pixel
!     real(KMAP),     dimension(:), allocatable :: signal
    real(SP),     dimension(:), allocatable :: signal
    integer(I4B) :: status
    integer(I4B) :: type_fits, nmaps_fits, maxpix, minpix
    CHARACTER(LEN=80)  :: units1
    CHARACTER(LEN=80),dimension(1:10)  :: unitsm
!    CHARACTER(LEN=80),dimension(:), allocatable  :: unitsm
    integer(I4B) :: extno_i
    !-----------------------------------------------------------------------

    units1    = ' '
    unitsm(:) = ' '
    fmiss_effct = 0.
    imissing = 0
    if (PRESENT(fmissval)) fmiss_effct = fmissval
    if (PRESENT(header)) then
       nlheader = size(header)
    else
       nlheader = 0
    endif
    extno_i = 0
    if (present(extno)) extno_i = extno

    obs_npix = getsize_fits(filename, nmaps = nmaps_fits, type=type_fits, extno=extno_i)

    !       if (nmaps_fits > nmaps) then 
    !          print*,trim(filename)//' only contains ',nmaps_fits,' maps'
    !          print*,' You try to read ',nmaps
    !       endif
    if (type_fits == 0 .or. type_fits == 2) then ! full sky map (in image or binary table)
       if (present(header)) then
          call read_bintab(filename, map(0:,1:), &
               & npixtot, nmaps, fmissing, anynull, header=header(1:), &
               & units=unitsm(1:), extno=extno_i)
       else
          call read_bintab(filename, map(0:,1:), &
               & npixtot, nmaps, fmissing, anynull, units=unitsm(1:), &
               &  extno=extno_i)
       endif
       if (present(units)) then
          units(1:size(units)) = unitsm(1:size(units))
       endif
       do imap = 1, nmaps
          anynull = .true.
          if (anynull) then
             imissing = 0
             do i=0,npixtot-1
                if ( ABS(map(i,imap)/fmissing -1.) < 1.e-5 ) then
                   map(i,imap) = fmiss_effct
                   imissing = imissing + 1
                endif
             enddo
          endif
       enddo
       if (imissing > 0) write(*,'(a,1pe11.4)') 'blank value : ' ,fmissing
    else if (type_fits == 3 .and. (nmaps == 1 .or. nmaps == 3)) then
       do imap = 1, nmaps
          obs_npix = getsize_fits(filename, extno = imap-1)       
          ! one partial map (in binary table with explicit pixel indexing)
          allocate(pixel(0:obs_npix-1), stat = status)
          allocate(signal(0:obs_npix-1), stat = status)
          call read_fits_cut4(filename, obs_npix, pixel, signal, header=header, units=units1, extno=imap-1)
          if (present(units)) units(imap) = trim(units1)
          maxpix = maxval(pixel)
          minpix = maxval(pixel)
          if (maxpix > (npixtot-1) .or. minpix < 0) then
             print*,'map constructed from file '//trim(filename)//', with pixels in ',minpix,maxpix
             print*,' wont fit in array with ',npixtot,' elements'
             call fatal_error
          endif
          map(:,imap)        = fmiss_effct
          map(pixel(:),imap) = signal(:)
          imissing = npixtot - obs_npix
          deallocate(signal)
          deallocate(pixel)
       enddo
    else 
       print*,'Unable to read the ',nmaps,' required  map(s) from file '//trim(filename)
       print*,'file type = ',type_fits
       call fatal_error
    endif
    !-----------------------------------------------------------------------
    if (imissing > 0) then
       write(*,'(i7,a,f7.3,a,1pe11.4)') &
            &           imissing,' missing pixels (', &
            &           (100.*imissing)/npixtot,' %),'// &
            &           ' have been set to : ',fmiss_effct
    endif

!    deallocate(unitsm)

    RETURN
  END subroutine input_map_KLOAD
  !=======================================================================
  subroutine read_bintab_KLOAD(filename, map, npixtot, nmaps, nullval, anynull, header, units, extno)
  !=======================================================================
    !     Read a FITS file
    !     This routine is used for reading MAPS by anafast.
    !     modified Feb 03 for units argument to run with Sun compiler
    !     Jan 2005, EH, improved for faster writting
    !=======================================================================
    character(len=*),                          intent(IN)  :: filename
    integer(I4B),                              intent(IN)  :: npixtot, nmaps
    real(KMAP),      dimension(0:,1:),         intent(OUT) :: map
    real(KMAP),                                intent(OUT) :: nullval
    logical(LGT),                              intent(OUT) :: anynull
    character(LEN=*), dimension(1:), optional, intent(OUT) :: header
    character(LEN=*), dimension(1:), optional, intent(OUT) :: units
    integer(I4B)                   , optional, intent(IN) :: extno

    integer(I4B) :: nl_header, len_header, nl_units, len_units
    integer(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    integer(I4B) :: group, firstpix, i, npix_old
    real(KMAP)   :: blank, testval
    real(DP)     :: bscale,bzero
    character(len=80) :: comment
    logical(LGT) :: extend
    integer(I4B) :: nmove, hdutype
    integer(I4B) :: frow, imap
    integer(I4B) :: datacode, width
    LOGICAL(LGT) ::  anynull_i

    integer(I4B),     parameter            :: maxdim=20 !number of columns in the extension
    integer(i4b), dimension(1:maxdim) :: npix, repeat
    integer(i8b), dimension(1:maxdim) :: i0, i1
    integer(i4b)                      :: nrow2read, nelem

    integer(I4B)                           :: nrows, tfields, varidat
    character(len=20), dimension(1:maxdim) :: ttype, tform, tunit
    character(len=20)                      :: extname
    !-----------------------------------------------------------------------
    status=0

    unit = 150
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    bscale = 1.0d0
    bzero = 0.0d0
    blank = -2.e25
    nullval = bscale*blank + bzero

    nl_header = 0
    if (present(header)) then
       nl_header = size(header)
       len_header = 80
    endif

    nl_units = 0
    if (present(units)) then
       nl_units = size(units)
       len_units = min(80,len(units(1))) ! due to SUN compiler bug
    endif

    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    if (naxis > 0) then ! there is an image
       !        determine the size of the image (look naxis1 and naxis2)
       call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

       !        check that it found only NAXIS1
       if (nfound == 2 .and. naxes(2) > 1) then
          print *,'multi-dimensional image'
          print *,'expected 1-D data.'
          call fatal_error
       end if

       if (nfound < 1) then
          call printerror(status)
          print *,'can not find NAXIS1.'
          call fatal_error
       endif

       npix(1)=naxes(1)
       if (npix(1) /= npixtot) then
          print *,'WARNING: found ',npix(1),' pixels in '//trim(filename)
          print *,'         expected ',npixtot
          npix(1) = min(npix(1), npixtot)
          print *,'         only ',npix(1),' will be read'
       endif

       call ftgkyd(unit,'BSCALE',bscale,comment,status)
       if (status == 202) then ! BSCALE not found
          bscale = 1.0d0
          status = 0
       endif
       call ftgkyd(unit,'BZERO', bzero, comment,status)
       if (status == 202) then ! BZERO not found
          bzero = 0.0d0
          status = 0
       endif
       call f90ftgky_(unit, 'BLANK', blank, comment, status)
!        if (KMAP == SP) then
!           call ftgkye(unit, 'BLANK', sdummy, comment, status) ; blank = sdummy
!        endif
!        if (KMAP == DP) then
!           call ftgkyd(unit, 'BLANK', ddummy, comment, status) ; blank = ddummy
!        endif
       if (status == 202) then ! BLANK not found 
          ! (according to fitsio BLANK is integer)
          blank = -2.e25
          status = 0
       endif
       nullval = bscale*blank + bzero

       !        -----------------------------------------

       group=1
       firstpix = 1
       call f90ftgpv_(unit, group, firstpix, npix(1), nullval, map(0:npix(1)-1,1), anynull, status)
       ! if there are any NaN pixels, (real data)
       ! or BLANK pixels (integer data) they will take nullval value
       ! and anynull will switch to .true.
       ! otherwise, switch it by hand if necessary
       testval = 1.e-6 * ABS(nullval)
       do i=0, npix(1)-1
          if (ABS(map(i,1)-nullval) < testval) then
             anynull = .true.
             goto 111
          endif
       enddo
111    continue

    else if (extend) then ! there is an extension
       nmove = +1
       if (present(extno)) nmove = +1 + extno
       call ftmrhd(unit, nmove, hdutype, status)

       call assert (hdutype==2, 'this is not a binary table')

       !        reads all the keywords
       call ftghbn(unit, maxdim, &
            &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &        status)

       if (tfields < nmaps) then
          print *,'found ',tfields,' maps in file '//trim(filename)
          print *,'expected ',nmaps
          call fatal_error
       endif

       !        finds the bad data value
!        if (KMAP == SP) then
!           call ftgkye(unit, 'BAD_DATA', sdummy, comment, status) ; nullval = sdummy
!        endif
!        if (KMAP == DP) then
!           call ftgkyd(unit, 'BAD_DATA', ddummy, comment, status) ; nullval = ddummy
!        endif
       call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
       if (status == 202) then ! bad_data not found
          if (KMAP == SP) nullval = s_bad_value ! default value
          if (KMAP == DP) nullval = d_bad_value ! default value
          status = 0
       endif

       if (nl_header > 0) then
          do i=1,nl_header
             header(i)(1:len_header) = ""
          enddo
          call get_clean_header(unit, header, filename, status)
          status = 0
       endif

       if (nl_units > 0) then
          do i=1,nl_units
             units(i)(1:len_units) = 'unknown' ! default
          enddo
          do imap = 1, min(nmaps, nl_units)
             units(imap)(1:len_units) = adjustl(tunit(imap))
          enddo
       endif

       npix_old = npixtot
       do imap = 1, nmaps
          !parse TFORM keyword to find out the length of the column vector
          call ftbnfm(tform(imap), datacode, repeat(imap), width, status)
          npix(imap) = nrows * repeat(imap)
          if (npix(imap) /= npixtot .and. npix_old /= npix(imap)) then
             print *,'WARNING: found ',npix(imap),' pixels in '//trim(filename)//', column ',imap
             print *,'         expected ',npixtot,' or ',npix_old
             npix_old = npix(imap)
             npix(imap) = min(npix(imap), npixtot)
             print *,'         only  ',npix(imap),' will be read'
          endif
       enddo

       call ftgrsz(unit, nrow2read, status)
       nrow2read = max(nrow2read, 1)
       firstpix  = 1  ! starting position in FITS within row, 1 based
       i0(:) = 0_i8b  ! starting element in array, 0 based
       do frow = 1, nrows, nrow2read
          do imap = 1, nmaps
             i1(imap) = min(i0(imap) + nrow2read * repeat(imap), int(npix(imap),i8b)) - 1_i8b
             nelem = i1(imap) - i0(imap) + 1
             call f90ftgcv_(unit, imap, frow, firstpix, nelem, &
                  & nullval, map(i0(imap):i1(imap),imap), anynull_i, status)
             anynull = anynull .or. anynull_i
             i0(imap) = i1(imap) + 1_i8b
          enddo
       enddo
       ! sanity check
       do imap = 1, nmaps
          if (i0(imap) /= npix(imap)) then
             call fatal_error('something wrong during piece wise reading')
          endif
       enddo

    else ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension')
    endif
    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_bintab_KLOAD

  !=======================================================================
  subroutine read_conbintab_KLOAD(filename, alms, nalms, units, extno)
    !=======================================================================
    !     Read a FITS file containing alms values
    !
    !     slightly modified to deal with vector column 
    !     in binary table       EH/IAP/Jan-98
    !
    !     Used by synfast when reading a binary file with alms for cons.real.
    !                        FKH/Apr-99
    !
    !   extno : optional, number of extension to be read, default=0: first extension
    !     Jan 2005, EH, improved for faster reading
    !=======================================================================
    CHARACTER(LEN=*),                   INTENT(IN) :: filename
    INTEGER(I4B),                       INTENT(IN) :: nalms !, nlheader
    REAL(KMAP), DIMENSION(0:nalms-1,1:6), INTENT(OUT) :: alms
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units
    INTEGER(I4B)                   , optional,  INTENT(IN) :: extno

    REAL(KMAP)                                        :: nullval
    LOGICAL(LGT)                                    ::  anynull

    INTEGER(I4B), DIMENSION(:), allocatable :: lm
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: npix
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: frow, imap
    INTEGER(I4B) :: datacode, repeat, width
    integer(I4B) :: i, l, m
    integer(i4b) :: nrow2read, nelem
    integer(i8b) :: i0, i1

    INTEGER(I4B), PARAMETER :: maxdim=20 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20), dimension(1:maxdim) :: ttype, tform, tunit
    CHARACTER(LEN=20) :: extname
    character(len=*), parameter :: code="read_conbintab"

    !-----------------------------------------------------------------------
    status=0

    unit = 150
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    alms=0.  ! set output to 0.
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    if (.not.extend) then
       print*,'No extension!'
       call fatal_error
    endif

    ! go to assigned extension (if it exists)
    nmove = +1
    if (present(extno)) nmove = +1 + extno
    call ftmrhd(unit, nmove, hdutype, status)
    if (status > 0) then
       ! if required extension not found: 
       ! print a warning, fill with dummy values, return to calling routine
       print*,code//' WARNING: the extension ',extno,' was not found in ',trim(filename)
       alms(0:nalms-1,1)=0.    ! l = 0
       alms(0:nalms-1,2)=1.    ! m = 1
       alms(0:nalms-1,3:6)=0.
       status = 0
       call ftclos(unit, status)        !     close the file
       return
    endif

    if (hdutype /= 2) then ! not a binary table
       print*, 'this is not a binary table'
       call fatal_error
    endif

    !        reads all the keywords
    call ftghbn(unit, maxdim, &
         &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
         &        status)

    if ((tfields/=3).and.(tfields/=5)) then
       print *,'found ',tfields,' columns in the file'
       print *,'expected 3 or 5'
       call fatal_error
    endif
    !        finds the bad data value
!     if (KMAP == SP) call ftgkye(unit,'BAD_DATA',nullval,comment,status)
!     if (KMAP == DP) call ftgkyd(unit,'BAD_DATA',nullval,comment,status)
    call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (KMAP == SP) nullval = s_bad_value ! default value
       if (KMAP == DP) nullval = d_bad_value ! default value
       status = 0
    endif

    !          if (nlheader > 0) then
    !             header = ""
    !             status = 0
    !             call get_clean_header(unit, header, filename, status)
    !          endif

    if (present(units)) then
       units(1) = tunit(2)
    endif

    !parse TFORM keyword to find out the length of the column vector
    call ftbnfm(tform(1), datacode, repeat, width, status)
    npix = nrows * repeat
    if (npix /= nalms) then
       print *,'found ',npix,' alms'
       print *,'expected ',nalms
       call fatal_error
    endif

    call ftgrsz(unit, nrow2read, status)
    nrow2read = max(nrow2read, 1)
    nelem = nrow2read * repeat
    i0 = 0_i8b
    allocate(lm(0:nelem-1))
    do frow = 1, nrows, nrow2read
       i1 = min(i0 + nrow2read * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       ! first column -> index
       call ftgcvj(unit, 1, frow, 1, nelem, i_bad_value, &
         &        lm(0), anynull, status)
       call assert(.not. anynull, 'There are undefined values in the table!')
       ! other columns -> a(l,m)
       do imap = 2, tfields
          call f90ftgcv_(unit, imap, frow, 1, nelem, nullval, &
               &        alms(i0:i1,imap+1), anynull, status)
          call assert (.not. anynull, 'There are undefined values in the table!')
       enddo
       ! recoding of the mapping, EH, 2002-08-08
       ! lm = l*l + l + m + 1
       do i = i0, i1
          l = int(sqrt(   real(lm(i-i0)-1, kind = DP)  ) )
          m = lm(i-i0) - l*(l+1) - 1
          ! check that round-off did not screw up mapping
          if (abs(m) > l) then
             print*,'Inconsistent l^2+l+m+1 -> l,m mapping'
             print*,l, m, l*(l+1)+m+1, lm(i-i0)
             call fatal_error
          endif
          alms(i,1) = real( l, kind=KMAP)
          alms(i,2) = real( m, kind=KMAP)
       enddo
       i0 = i1 + 1_i8b
    enddo
    deallocate(lm)
    ! sanity check
    if (i0 /= npix) then
       print*,'something wrong during piece-wise reading'
       call fatal_error
    endif

!     !reads the columns
!     ! first column : index -> (l,m)
!     allocate(lm(0:nalms-1))
!     column = 1
!     frow = 1
!     firstpix = 1
!     npix = nrows * repeat
!     if (npix /= nalms) then
!        print *,'found ',npix,' alms'
!        print *,'expected ',nalms
!        call fatal_error
!     endif
!     call ftgcvj(unit, column, frow, firstpix, npix, i_bad_value, &
!          &        lm(0), anynull, status)
!     call assert (.not. anynull, 'There are undefined values in the table!')
!     ! recoding of the mapping, EH, 2002-08-08
!     ! lm = l*l + l + m + 1
!     do i = 0, nalms - 1
!        l = int(sqrt(   real(lm(i)-1, kind = DP)  ) )
!        m = lm(i) - l*(l+1) - 1
!        ! check that round-off did not screw up mapping
!        if (abs(m) > l) then
!           print*,'Inconsistent l^2+l+m+1 -> l,m mapping'
!           print*,l, m, l*(l+1)+m+1, lm(i)
!           call fatal_error
!        endif
!        alms(i,1) = real( l, kind=KMAP )
!        alms(i,2) = real( m, kind=KMAP )
!     enddo
!     deallocate(lm)
!
!     do imap = 2, tfields
!        !parse TFORM keyword to find out the length of the column vector
!        call ftbnfm(tform(imap), datacode, repeat, width, status)
!
!        !reads the columns
!        column = imap
!        frow = 1
!        firstpix = 1
!        npix = nrows * repeat
!        if (npix /= nalms) then
!           print *,'found ',npix,' alms'
!           print *,'expected ',nalms
!           call fatal_error
!        endif
!        call f90ftgcv_(unit, column, frow, firstpix, npix, nullval, &
!             &        alms(0:npix-1,imap+1), anynull, status)
!        call assert (.not. anynull, 'There are undefined values in the table!')
!     enddo
    !     close the file
    call ftclos(unit, status)


    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)
    return
  end subroutine read_conbintab_KLOAD

  !=======================================================================
  subroutine write_bintab_KLOAD(map, npix, nmap, header, nlheader, filename, extno)
    !=======================================================================
    !     Create a FITS file containing a binary table extension with 
    !     the temperature map in the first column
    !     written by EH from writeimage and writebintable 
    !     (fitsio cookbook package)
    !
    !     slightly modified to deal with vector column (ie TFORMi = '1024E')
    !     in binary table       EH/IAP/Jan-98
    !
    !     simplified the calling sequence, the header sould be filled in
    !     before calling the routine
    !
    !     July 21, 2004: SP version
    !     Jan 2005, EH, improved for faster writting
    !=======================================================================
    INTEGER(I4B),     INTENT(IN) :: npix, nmap, nlheader
!     REAL(KMAP),       INTENT(IN), DIMENSION(0:npix-1,1:nmap), target :: map
    REAL(KMAP),       INTENT(IN), DIMENSION(0:npix-1,1:nmap) :: map
    CHARACTER(LEN=*), INTENT(IN), DIMENSION(1:nlheader)      :: header
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    INTEGER(I4B)    , INTENT(IN), optional :: extno

    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 20 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    integer(i4b) :: repeat, nrow2write, nelem
    integer(i8b) :: i0, i1
    INTEGER(I4B) :: frow,  felem, colnum, hdutype
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn, extno_i
    character(len=filenamelen) sfilename
    character(len=1)          :: pform 
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'
    status=0
    unit = 100
    blocksize=1

    extno_i = 0
    if (present(extno)) extno_i = extno

    if (extno_i == 0) then
       !*************************************
       !     create the new empty FITS file
       !*************************************
       call ftinit(unit,filename,blocksize,status)

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       !     writes supplementary keywords : none

       !     write the current date
       call ftpdat(unit,status) ! format (yyyy-mm-dd)

       !     ----------------------
       !     image : none
       !     ----------------------

       !     ----------------------
       !     extension
       !     ----------------------
    else

       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1,blocksize,status)
       call ftmahd(unit,1+extno_i,hdutype,status)

    endif

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    tfields  = nmap
    repeat   = 1024
    tform(1:nmap) = '1024'//pform
    if (npix < 1024) then ! for nside <= 8
       repeat = 1
       tform(1:nmap) = '1'//pform
    endif
    nrows    = npix / repeat ! naxis1
    ttype(1:nmap) = 'simulation'   ! will be updated
    tunit(1:nmap) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORM1 at the desired place
    if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
    if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          if (itn > tfields) goto 10
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          call putrec(unit,header(i), status)    ! write new TTYPE1
          status = 0
          call ftpkys(unit,'TFORM'//stn,tform(1),comment,status) ! and write new TFORM1 right after
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
10     continue
       status = 0
    enddo

    !     write the extension buffer by buffer
    call ftgrsz(unit, nrow2write, status)
    nrow2write = max(nrow2write, 1)
    felem  = 1  ! starting position in FITS (element), 1 based
    i0 = 0_i8b  ! starting element in array, 0 based
    do frow = 1, nrows, nrow2write
       i1 = min(i0 + nrow2write * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       do colnum = 1, nmap
          call f90ftpcl_(unit, colnum, frow, felem, nelem, map(i0:i1, colnum), status)
       enddo
       i0 = i1 + 1_i8b
    enddo
    ! sanity check
    if (i0 /= npix) then
       call fatal_error("something wrong during piece wise writting")
    endif

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_bintab_KLOAD
  !=======================================================================
  subroutine write_asctab_KLOAD &
       &  (clout, lmax, ncl, header, nlheader, filename)
    !=======================================================================
  !     Create a FITS file containing an ASCII table extension with 
  !     the measured power spectra
  !     written by EH from writeimage and writeascii 
  !     (fitsio cookbook package)
  !
  !     clout = power spectra with l in (0:lmax)
  !     ncl = number of spectra
  !     header = FITS header to be put on top of the file
  !     nlheader = number of lines of the header
  !     filename = FITS output file name
  !=======================================================================
    !
    INTEGER(I4B),      INTENT(IN)           :: lmax, ncl,nlheader
    REAL(KMAP),         INTENT(IN)           ::  clout(0:lmax,1:ncl)
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)           :: filename

    INTEGER(I4B) :: bitpix,naxis,naxes(1)
    LOGICAL(LGT) :: simple,extend
    CHARACTER(LEN=10) ::  card

    INTEGER(I4B), PARAMETER :: nclmax = 20
    INTEGER(I4B) ::  status,unit,blocksize,tfields,nrows,rowlen
    INTEGER(I4B) ::  nspace,tbcol(nclmax),colnum,frow,felem
    CHARACTER(LEN=16) :: ttype(nclmax),tform(nclmax),tunit(nclmax),extname
    CHARACTER(LEN=80) :: comment, card_tbcol
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn, i
    character(len=6)        :: form
    !=======================================================================
    if (KMAP == SP) form = 'E15.7'
    if (KMAP == DP) form = 'D24.15'

    status=0

    unit = 87

    !     open the FITS file, with write access
    !      readwrite=1
    !      call ftopen(unit,filename,readwrite,blocksize,status)

    blocksize=1
    call ftinit(unit,filename,blocksize,status)

    !     -----------------------------------------------------
    !     initialize parameters about the FITS image
    simple=.true.
    bitpix=32     ! integer*4
    naxis=0       ! no image
    naxes(1)=0
    extend=.true. ! there is an extension

    !     ----------------------
    !     primary header
    !     ----------------------
    !     write the required header keywords
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !     writes supplementary keywords : none

    !     write the current date
    call ftpdat(unit,status) ! (format ccyy-mm-dd)

    !     ----------------------
    !     image : none
    !     ----------------------

    !     ----------------------
    !     extension
    !     ----------------------

    !     append a new empty extension onto the end of the primary array
    call ftcrhd(unit,status)

    !     define parameters for the ASCII table
    nrows   = lmax+1
    tfields = ncl
    tform(1:ncl) = trim(adjustl(form))
    ttype(1:ncl) = 'power spectrum' ! is updated by the value given in the header
    tunit(1:ncl) = '' ! optional, will not appear
    extname      = '' ! optional, will not appear

    !     calculate the starting position of each column, and the total row length
    nspace=1
    call ftgabc(tfields,tform,nspace,rowlen,tbcol,status)

    !     write the required header parameters for the ASCII table
    call ftphtb(unit,rowlen,nrows,tfields,ttype,tbcol,tform,tunit, &
         &            extname,status)

    !     write the header literally, putting TFORM1 at the desired place
    comment = ''
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          ! discard at their original location:
          call ftdkey(unit,card(1:6),status)  !         old TTYPEi
          status = 0
          call ftdkey(unit,'TFORM'//stn,status) !           TFORMi
          status = 0
          call ftgcrd(unit,'TBCOL'//stn,card_tbcol,status)
          status = 0
          call ftdkey(unit,'TBCOL'//stn,status) !           TBCOLi
          ! and rewrite
          status = 0
          call putrec(unit,card_tbcol,status) !             TBCOLi
          status = 0
          call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! TFORMi right after
          status = 0
          call putrec(unit,header(i), status)   !           TTYPEi
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
10     continue
       status = 0
    enddo

    frow=1
    felem=1
    do colnum = 1, ncl
       call f90ftpcl_(unit, colnum, frow, felem, nrows, clout(0:nrows-1,colnum), status)  
    enddo

    !     close the FITS file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_asctab_KLOAD
  !=======================================================================
  !   DUMP_ALM
  !
  !     Create/extend a FITS file containing a binary table extension with 
  !     the a_lm coefficients in each extension.
  !
  !     Format of alm FITS file: in each extension, 3 columns: 
  !            index=l*l+l+m+1, real(a_lm), imag(a_lm)
  !      the a_lm are obtained using the complex Y_lm
  !      only the modes with m>=0 are stored (because the map is assumed real)
  !
  !     The extensions contain, in general, T, E (= G) and B (= C) in that order
  !     First extension has extno = 0
  !
  !     Adapted from write_bintab, FKH/Apr-99
  !   SP/DP overloaded, 2004-12, EH
  !    reduced size of intermediate arrays (alms_out, lm)
  !
  subroutine dump_alms_KLOAD(filename, alms, nlmax, header, nlheader, extno)
    !=======================================================================
    INTEGER(I4B),      INTENT(IN) :: nlmax, nlheader, extno
    COMPLEX(KMAPC),      INTENT(IN), DIMENSION(0:,0:) :: alms
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename
    ! local variables
    INTEGER(I4B),   DIMENSION(:),   allocatable :: lm
    REAL(KMAP),     DIMENSION(:,:), allocatable :: alms_out
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i,l,m,cnt,hdutype, nmmax
    LOGICAL(LGT) ::  simple,extend, found_lmax, found_mmax
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 20 !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat
    INTEGER(I4B) :: frow,  felem
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn
    character(len=filenamelen) sfilename
    character(len=1) :: pform
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'
    
    nmmax = size(alms,2) - 1
    if (nmmax < 0 .or. nmmax > nlmax) call fatal_error('inconsistent Lmax and Mmax in dump_alms')
!!    npix=((nlmax+1)*(nlmax+2))/2
    npix = ((nmmax+1)*(2*nlmax-nmmax+2))/2

    status=0
    unit = 100
    blocksize=1

    if (extno==0) then
       !*********************************************
       !     create the new empty FITS file
       !*********************************************
       call ftinit(unit,filename,blocksize,status)

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       !     writes supplementary keywords : none

       !     write the current date
       call ftpdat(unit,status) ! format (yyyy-mm-dd)

       !     ----------------------
       !     image : none
       !     extension
       !     ----------------------
    else

       !*********************************************
       !     reopen an existing file and go to the end
       !*********************************************
       ! remove the leading '!' (if any) when reopening the same file
       sfilename = adjustl(filename)
       if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)
       call ftopen(unit,sfilename,1,blocksize,status)
       call ftmahd(unit,1+extno,hdutype,status)

    endif

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    nrows    = npix  ! naxis1
    tfields  = 3
    tform(1)='1J'
    tform(2:3) = '1'//pform
    ttype(1) = 'index=l^2+l+m+1'
    ttype(2) = 'alm (real)'
    ttype(3) = 'alm (imaginary)' 
    tunit(1:3) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORM1 at the desired place
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          call putrec(unit,header(i), status)           ! write new TTYPE1
          if (itn==1) then
             comment = 'data format of field: 4-byte INTEGER'
          else
             if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
             if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
          endif
          status = 0
          call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    call ftukyj(unit, 'MIN-LPOL', 0,     'Minimum L multipole order',  status)
    call ftukyj(unit, 'MAX-LPOL', nlmax, 'Maximum L multipole order',  status)
    call ftukyj(unit, 'MAX-MPOL', nmmax, 'Maximum M multipole degree', status)

    allocate(lm      (0:nmmax))
    allocate(alms_out(0:nmmax,1:2))
    !     write the extension one column by one column
    frow   = 1  ! starting position (row)
    felem  = 1  ! starting position (element)
    do l = 0, nlmax
       cnt = 0
       do m = 0, min(l,nmmax)
          lm(cnt) = l**2 + l + m + 1
          alms_out(cnt,1)=REAL( alms(l,m))
          alms_out(cnt,2)=AIMAG(alms(l,m))
          cnt = cnt + 1
       enddo
       call ftpclj(unit, 1, frow, felem, cnt, lm(0),         status)
       call f90ftpcl_(unit, 2, frow, felem, cnt, alms_out(0:cnt-1,1), status)
       call f90ftpcl_(unit, 3, frow, felem, cnt, alms_out(0:cnt-1,2), status)
       frow = frow + cnt
    enddo
    deallocate(lm)
    deallocate(alms_out)

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)


    return
  end subroutine dump_alms_KLOAD
  !=======================================================================
  subroutine write_alms_KLOAD(filename, nalms, alms, ncl, header, nlheader, extno)
    !=======================================================================
    !     Writes alms from to binary FITS file, FKH/Apr-99
    !     ncl is the number of columns, in the output fits file,
    !     either 3 or 5 (with or without errors respectively)
    !
    !    input array (real)                   FITS file
    !     alms(:,1) = l                      )
    !     alms(:,2) = m                      )---> col 1: l*l+l+m+1
    !     alms(:,3) = real(a_lm)              ---> col 2
    !     alms(:,4) = imag(a_lm)              ---> col 3
    !     alms(:,5) = real(delta a_lm)        ---> col 4
    !     alms(:,6) = imag(delta a_lm)        ---> col 5
    !
    !=======================================================================
    
    INTEGER(I4B), INTENT(IN) :: nalms, nlheader, ncl, extno
    REAL(KMAP),        INTENT(IN), DIMENSION(0:nalms-1,1:(ncl+1)) :: alms
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)               :: filename

    INTEGER(I4B), DIMENSION(:), allocatable :: lm
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i,hdutype, lmax, mmax, lmin
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 20 !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat, repeat
    INTEGER(I4B) :: frow,  felem, colnum, stride, istart, iend, k
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn
    integer(I4B) :: l, m
    character(len=1) :: pform
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'

    status=0
    unit = 100

    !     create the new empty FITS file
    blocksize=1

    if (extno==1) then

       call ftinit(unit,filename,blocksize,status)

       !     -----------------------------------------------------
       !     initialize parameters about the FITS image
       simple=.true.
       bitpix=32     ! integer*4
       naxis=0       ! no image
       naxes(1)=0
       extend=.true. ! there is an extension

       !     ----------------------
       !     primary header
       !     ----------------------
       !     write the required header keywords
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       !     writes supplementary keywords : none

       !     write the current date
       call ftpdat(unit,status) ! format ccyy-mm-dd

       !     ----------------------
       !     image : none
       !     ----------------------

       !     ----------------------
       !     extension
       !     ----------------------

    else

       call ftopen(unit,filename,1,blocksize,status)
       call ftmahd(unit,extno,hdutype,status)

    endif

    !     creates an extension
    call ftcrhd(unit, status)

    !     writes required keywords
    nrows    = nalms  ! naxis1
    tfields  = ncl
    repeat   = 1
    tform(1)='1J'
    tform(2:ncl) = '1'//pform
    ttype(1) = 'index=l^2+l+m+1'
    ttype(2) = 'alm (real)'
    ttype(3) = 'alm (imaginary)' 
    if (ncl>3) then
       ttype(4) = 'error (real)'
       ttype(5) = 'error (imaginary)'
    endif
    tunit(1:ncl) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)

    !     write the header literally, putting TFORM1 at the desired place
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          call putrec(unit,header(i), status)           ! write new TTYPE1
          if (itn==1) then
             comment = 'data format of field: 4-byte INTEGER'
          else
             if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
             if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
          endif
          status = 0
          call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo

    lmax = nint(maxval(alms(:,1)))
    lmin = nint(minval(alms(:,1)))
    mmax = nint(maxval(alms(:,2)))
    call ftukyj(unit, 'MIN-LPOL', lmin, 'Minimum L multipole order',  status)
    call ftukyj(unit, 'MAX-LPOL', lmax, 'Maximum L multipole order',  status)
    call ftukyj(unit, 'MAX-MPOL', mmax, 'Maximum M multipole degree', status)

    !     write the extension by blocks of rows ! EH, Dec 2004
    felem  = 1  ! starting position (element)
!!!    stride = 1000 ! 1000 rows at a time
    call ftgrsz(unit, stride, status) ! find optimal stride in rows
    stride = max( stride, 1)
    allocate(lm(0:stride-1))
    do k = 0, (nalms-1)/(stride * repeat)
       istart = k * (stride * repeat)
       iend   = min(nalms, istart + stride * repeat) - 1
       do i = istart, iend ! recode the (l,m) -> lm mapping, EH, 2002-08-08
          l = nint(alms(i,1))
          m = nint(alms(i,2))
          lm(i-istart) = (l+1)*l + m + 1
       enddo

       frow = istart/repeat + 1
       npix = iend - istart + 1
       call ftpclj(unit, 1, frow, felem, npix, lm(0), status)
       do colnum = 2, ncl
          call f90ftpcl_(unit, colnum, frow, felem, npix, alms(istart:iend,colnum+1), status)
       enddo
    enddo
    deallocate(lm)

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine write_alms_KLOAD
  !=======================================================================
  subroutine read_alms_KLOAD(filename, nalms, alms, ncl, header, nlheader, extno) 
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column 
    !     in binary table       EH/IAP/Jan-98
    !
    !     Used by synfast when reading a binary file with alms for cons.real.
    !                        FKH/Apr-99
    !
    !     called by fits2alms
    !     Jan 2005, EH, improved for faster reading
    !=======================================================================
    CHARACTER(LEN=*),               INTENT(IN) :: filename
    INTEGER(I4B),                        INTENT(IN) :: nalms, ncl,nlheader,extno
    REAL(KMAP), DIMENSION(0:nalms-1,1:(ncl+1)), INTENT(OUT) :: alms
    CHARACTER(LEN=80), INTENT(OUT), DIMENSION(1:nlheader) :: header
    REAL(KMAP)                                        :: nullval
    LOGICAL(LGT)                                    ::  anynull

    INTEGER(I4B), DIMENSION(:), allocatable :: lm
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: npix
    CHARACTER(LEN=80) :: comment ! , record
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype ! , nkeys , nspace
    INTEGER(I4B) :: frow, imap
    INTEGER(I4B) :: datacode, repeat, width
    integer(I4B) :: i, l, m
    integer(i4b) :: nrow2read, nelem
    integer(i8b) :: i0, i1

    INTEGER(I4B), PARAMETER :: maxdim=20 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname

    !-----------------------------------------------------------------------
    status=0
    header=''
    unit = 150
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    alms=0.
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    call assert (extend, 'No extension!')
    nmove = +extno
    call ftmrhd(unit, nmove, hdutype, status)
    !cc         write(*,*) hdutype

    call assert(hdutype==2, 'this is not a binary table')

    header = ""
    call get_clean_header( unit, header, filename, status)


    !        reads all the keywords
    call ftghbn(unit, maxdim, &
         &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
         &        status)

    if (tfields<ncl) then
       print *,'found ',tfields,' columns in the file'
       print *,'expected ',ncl
       call fatal_error
    endif
    !        finds the bad data value
!     if (KMAP == SP) call ftgkye(unit,'BAD_DATA',nullval,comment,status)
!     if (KMAP == DP) call ftgkyd(unit,'BAD_DATA',nullval,comment,status)
    call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (KMAP == SP) nullval = s_bad_value ! default value
       if (KMAP == DP) nullval = d_bad_value ! default value
       status = 0
    endif
    !parse TFORM keyword to find out the length of the column vector
    call ftbnfm(tform(1), datacode, repeat, width, status)
    npix = nrows * repeat
    if (npix /= nalms) then
       print *,'found ',npix,' alms'
       print *,'expected ',nalms
       call fatal_error
    endif

    call ftgrsz(unit, nrow2read, status)
    nrow2read = max(nrow2read, 1)
    nelem = nrow2read * repeat
    i0 = 0_i8b
    allocate(lm(0:nelem-1))
    do frow = 1, nrows, nrow2read
       i1 = min(i0 + nrow2read * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       ! first column -> index
       call ftgcvj(unit, 1, frow, 1, nelem, i_bad_value, &
         &        lm(0), anynull, status)
       call assert(.not. anynull, 'There are undefined values in the table!')
       ! other columns -> a(l,m)
       do imap = 2, ncl
          call f90ftgcv_(unit, imap, frow, 1, nelem, nullval, &
               &        alms(i0:i1,imap+1), anynull, status)
          call assert (.not. anynull, 'There are undefined values in the table!')
       enddo
       ! recoding of the mapping, EH, 2002-08-08
       ! lm = l*l + l + m + 1
       do i = i0, i1
          l = int(sqrt(   real(lm(i-i0)-1, kind = DP)  ) )
          m = lm(i-i0) - l*(l+1) - 1
          ! check that round-off did not screw up mapping
          if (abs(m) > l) then
             print*,'Inconsistent l^2+l+m+1 -> l,m mapping'
             print*,l, m, l*(l+1)+m+1, lm(i-i0)
             call fatal_error
          endif
          alms(i,1) = real( l, kind=KMAP)
          alms(i,2) = real( m, kind=KMAP)
       enddo
       i0 = i1 + 1_i8b
    enddo
    deallocate(lm)
    ! sanity check
    if (i0 /= npix) then
       print*,'something wrong during piece-wise reading'
       call fatal_error
    endif

    !     close the file
    call ftclos(unit, status)


    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)
    return
  end subroutine read_alms_KLOAD
  !**************************************************************************
  SUBROUTINE read_bintod_KLOAD(filename, tod, npixtot, ntods, firstpix, nullval, anynull, &
                          header, extno)
  !**************************************************************************
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column (ie TFORMi = '1024E')
    !     in binary table       EH/IAP/Jan-98
    !
    !     This routine is used for reading TODS by anafast.
    !     Modified to start at a given pix numb OD & RT 02/02
    !     Modified to handle huge array (npix_tot > 2^32) OD & EH 07/02
    !      2002-07-08 : bugs correction by E.H. 
    !=======================================================================
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*),               INTENT(IN)  :: filename
    INTEGER(I8B)   ,                INTENT(IN)  :: npixtot,firstpix
    INTEGER(I4B),                   INTENT(IN)  :: ntods
    REAL(KMAP), DIMENSION(0:,1:),     INTENT(OUT) :: tod
    REAL(KMAP),                       INTENT(OUT) :: nullval
    LOGICAL(LGT),                   INTENT(OUT) :: anynull
    character(len=*), dimension(1:),intent(out), optional :: header
    INTEGER(I4B),                   INTENT(IN),  OPTIONAL :: extno
    
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: npix_32 !,firstpix_32
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: column, frow, itod
    INTEGER(I4B) :: datacode, repeat, width
    
    INTEGER(I4B), PARAMETER :: maxdim=20 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat,felem
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim),tunit(maxdim), extname
 
    INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp, i0, i1
    
    !-----------------------------------------------------------------------
    status=0
    
    unit = 150
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .FALSE.
    
    readwrite=0
    CALL ftopen(unit,filename,readwrite,blocksize,status)
    IF (status .GT. 0) CALL printerror(status)
    !     -----------------------------------------
    
    !     determines the presence of image
    CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
    IF (status .GT. 0) CALL printerror(status)
    
    !     determines the presence of an extension
    CALL ftgkyl(unit,'EXTEND', extend, comment, status)
    IF (status .GT. 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code
    
    IF (naxis .GT. 0) THEN ! there is an image
       print*,'WARNING : Image is ignored in '//trim(filename)
    ENDIF
    IF (extend) THEN ! there is an extension

       nmove = +1
       if (present(extno)) nmove = +1 + extno
       CALL ftmrhd(unit, nmove, hdutype, status)
       !cc         write(*,*) hdutype

       call assert(hdutype==2, 'this is not a binary table')

       ! reads all the keywords
       CALL ftghbn(unit, maxdim, &
            &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &        status)

       IF (tfields .LT. ntods) THEN
          PRINT *,'found ',tfields,' tods in the file'
          PRINT *,'expected ',ntods
          call fatal_error
       ENDIF

       if (present(header)) then
          header = ""
          status = 0
          call get_clean_header(unit, header, filename, status)
       endif

       ! finds the bad data value
!        if (KMAP == SP) CALL ftgkye(unit,'BAD_DATA',nullval,comment,status)
!        if (KMAP == DP) CALL ftgkyd(unit,'BAD_DATA',nullval,comment,status)
       call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
       IF (status .EQ. 202) THEN ! bad_data not found
          if (KMAP == SP) nullval = s_bad_value ! default value
          if (KMAP == DP) nullval = d_bad_value ! default value
          status = 0
       ENDIF

       IF (npixtot .LT. nchunk_max) THEN

          DO itod = 1, ntods

             !parse TFORM keyword to find out the length of the column vector (repeat)
             CALL ftbnfm(tform(itod), datacode, repeat, width, status)
             frow = (firstpix)/repeat+1          ! 1 based 
             felem = firstpix-(frow-1)*repeat+1  ! 1 based 
             npix_32 = npixtot 

             !reads the columns
             column = itod
             CALL f90ftgcv_(unit, column, frow, felem, npix_32, nullval, &
                  &        tod(0:npix_32-1,itod), anynull, status)
          END DO

       ELSE

          q = (npixtot-1)/nchunk_max
          DO iq = 0,q
             IF (iq .LT. q) THEN
                npix_tmp = nchunk_max
             ELSE
                npix_tmp = npixtot - iq*nchunk_max
             ENDIF
             firstpix_tmp = firstpix + iq*nchunk_max
             npix_32 = npix_tmp
             i0 = firstpix_tmp-firstpix
             i1 = i0 + npix_tmp - 1_i8b

             DO itod = 1, ntods
                ! parse TFORM keyword to find out the length of the column vector
                CALL ftbnfm(tform(itod), datacode, repeat, width, status)
                frow = (firstpix_tmp)/repeat+1          ! 1 based 
                felem = firstpix_tmp-(frow-1)*repeat+1  ! 1 based 
                CALL f90ftgcv_(unit, itod, frow, felem, npix_32, nullval, &
                     &      tod(i0:i1,itod), anynull, status)
             END DO

          ENDDO

       ENDIF

    ELSE ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension')
    ENDIF

    ! close the file
    CALL ftclos(unit, status)

    ! check for any error, and if so print out error messages
    IF (status .GT. 0) CALL printerror(status)

    RETURN

  END SUBROUTINE read_bintod_KLOAD
  !=======================================================================
  
  !======================================================================================
  SUBROUTINE write_bintabh_KLOAD(tod, npix, ntod, header, nlheader, filename, extno, firstpix, repeat)
    !======================================================================================

    ! =================================================================================
    !     Create a FITS file containing a binary table extension in the first extension
    !
    !     Designed to deal with Huge file, (n_elements > 2^31)
    !
    !     OPTIONNAL NEW PARAMETERS:
    !     firstpix : position in the file of the first element to be written (start at 0) 
    !                default value =0
    !                8 bytes integer
    !                if NE 0 then suppose that the file already exists
    !
    !     repeat   : lenght of vector per unit rows and colomns of the first binary extension
    !                default value = 12000 (\equiv 1 mn of PLANCK/HFI data)
    !                4 byte integer
    ! 
    !     OTHER PARAMETERS
    !     unchanged as compare to the standard write_bintab of the HEALPIX package except 
    !     npix which is a 8 bytes integer
    !
    !     Adapted from write_bintab
    !                                           E.H. & O.D. @ IAP 07/02
    !
    !     Requires a compilation in 64 bits of the CFITSIO 
    !     Note that the flag -D_FILE_OFFSETS_BITS=64 has to be added 
    !         (cf page CFITIO 2.2 User's guide  Chap 4, section 4-13)
    ! 
    ! 2002-07-08 : bugs correction by E.H. 
    !    (uniform use of firstpix_tmp, introduction of firstpix_chunk)
    !==========================================================================================

    USE healpix_types
    IMPLICIT NONE

    INTEGER(I8B)     , INTENT(IN)           :: npix
    INTEGER(I8B)     , INTENT(IN), OPTIONAL :: firstpix
    INTEGER(I4B)     , INTENT(IN), OPTIONAL :: repeat
    INTEGER(I4B)     , INTENT(IN)           :: ntod,nlheader
    REAL(KMAP)       , INTENT(IN), DIMENSION(0:npix-1,1:ntod) :: tod
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader)      :: header
    CHARACTER(LEN=*),  INTENT(IN)           :: filename
    INTEGER(I4B)    , INTENT(IN)     , OPTIONAL :: extno

    INTEGER(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(1),repeat_tmp,repeat_fits
    INTEGER(I4B) :: i,npix_32
    LOGICAL(LGT) :: simple,extend
    CHARACTER(LEN=80) :: comment, ch

    INTEGER(I4B), PARAMETER :: maxdim = 20 !number of columns in the extension
    INTEGER(I4B)      :: nrows,tfields,varidat
    INTEGER(I4B)      :: frow,felem,colnum,readwrite,width,datacode,hdutype
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) :: card
    CHARACTER(LEN=2)  :: stn
    INTEGER(I4B)      :: itn  

    INTEGER(I4B)      :: extno_i
    character(len=filenamelen) :: sfilename
    INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp,firstpix_chunk, i0, i1
    character(len=1) :: pform
    !-----------------------------------------------------------------------

    if (KMAP == SP) pform = 'E'
    if (KMAP == DP) pform = 'D'

    IF (.NOT. PRESENT(repeat) ) THEN 
       repeat_tmp = 1
       if (mod(npix,1024_i8b) == 0) then
          repeat_tmp = 1024
       elseif (npix >= 12000) then
          repeat_tmp = 12000 
       endif
    ELSE 
       repeat_tmp = repeat
    ENDIF
    IF (.NOT. PRESENT(firstpix) ) THEN 
       firstpix_tmp = 0 
    ELSE 
       firstpix_tmp = firstpix
    ENDIF

    extno_i = 0
    if (present(extno)) extno_i = extno

    status=0
    unit = 100
    blocksize=1

    ! remove the leading '!' (if any) when reopening the same file
    sfilename = adjustl(filename)
    if (sfilename(1:1) == '!') sfilename = sfilename(2:filenamelen)

    ! create the new empty FITS file

    IF (firstpix_tmp .EQ. 0) THEN

       if (extno_i == 0) then 
          CALL ftinit(unit,filename,blocksize,status)

          ! -----------------------------------------------------
          ! Initialize parameters about the FITS image
          simple=.TRUE.
          bitpix=32     ! integer*4
          naxis=0       ! no image
          naxes(1)=0
          extend=.TRUE. ! there is an extension

          !     ----------------------
          !     primary header
          !     ----------------------
          !     write the required header keywords
          CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

          !     writes supplementary keywords : none

          !     write the current date
          CALL ftpdat(unit,status) ! format ccyy-mm-dd

       !     ----------------------
       !     image : none
       !     ----------------------

       !     ----------------------
       !     extension
       !     ----------------------
       else

          !*********************************************
          !     reopen an existing file and go to the end
          !*********************************************
          call ftopen(unit,sfilename,1,blocksize,status)
          call ftmahd(unit,1+extno_i,hdutype,status)

       endif

       !     creates an extension
       CALL ftcrhd(unit, status)

       !     writes required keywords
       nrows    = npix / repeat_tmp ! naxis1
       tfields  = ntod
       WRITE(ch,'(i8)') repeat_tmp
       tform(1:ntod) = TRIM(ADJUSTL(ch))//pform

       IF (npix .LT. repeat_tmp) THEN
          nrows = npix
          tform(1:ntod) = '1'//pform
       ENDIF
       ttype(1:ntod) = 'simulation'   ! will be updated
       tunit(1:ntod) = ''      ! optional, will not appear
       extname  = ''      ! optional, will not appear
       varidat  = 0

       CALL ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
            &     extname, varidat, status)

       !     write the header literally, putting TFORM1 at the desired place
       DO i=1,nlheader
          card = header(i)
          IF (card(1:5) == 'TTYPE') THEN ! if TTYPE1 is explicitely given
             stn = card(6:7)
             READ(stn,'(i2)') itn
             ! discard at their original location:
             CALL ftmcrd(unit,'TTYPE'//stn,'COMMENT',status)  ! old TTYPEi and 
             CALL ftmcrd(unit,'TFORM'//stn,'COMMENT',status)  !     TFORMi
             CALL ftprec(unit,header(i), status)           ! write new TTYPE1
             if (KMAP == SP) comment = 'data format of field: 4-byte REAL'
             if (KMAP == DP) comment = 'data format of field: 8-byte REAL'
             CALL ftpkys(unit,'TFORM'//stn,tform(1),comment,status) ! and write new TFORM1 right after
          ELSEIF (header(i).NE.' ') THEN
             CALL ftprec(unit,header(i), status)
          ENDIF
10        CONTINUE
       ENDDO

    ELSE
       ! The file already exists
       readwrite=1
       CALL ftopen(unit,sfilename,readwrite,blocksize,status)
       CALL ftmahd(unit,2+extno_i,hdutype,status) 

       CALL ftgkys(unit,'TFORM1',tform(1),comment,status)
       CALL ftbnfm(tform(1),datacode,repeat_fits,width,status)

       IF (repeat_tmp .NE. repeat_fits) THEN
          WRITE(*,*) 'WARNING routine write_bintabh'
          WRITE(*,*) 'Inexact repeat value. Use the one read in the file'
       ENDIF

    ENDIF


    IF (npix .LT. nchunk_max) THEN ! data is small enough to be written in one chunk

       frow = (firstpix_tmp)/repeat_tmp + 1
       felem = firstpix_tmp-(frow-1)*repeat_tmp+1
       npix_32 = npix 

       DO colnum = 1, ntod
          call f90ftpcl_(unit, colnum, frow, felem, npix_32, tod(0:npix_32-1,colnum), status)
       END DO

    ELSE ! data has to be written in several chunks

       q = (npix-1)/nchunk_max
       DO iq = 0,q
          IF (iq .LT. q) THEN
             npix_tmp = nchunk_max
          ELSE
             npix_tmp = npix - iq*nchunk_max
          ENDIF
          firstpix_chunk = firstpix_tmp + iq*nchunk_max
          frow  = (firstpix_chunk)/repeat_tmp+1
          felem =  firstpix_chunk-(frow-1)*repeat_tmp+1 
          npix_32 = npix_tmp
          i0 = firstpix_chunk - firstpix_tmp
          i1 = i0 + npix_tmp - 1_i8b
          DO colnum = 1, ntod
             call f90ftpcl_(unit, colnum, frow, felem, npix_32, &
                  &          tod(i0:i1, colnum), status)
          END DO
       ENDDO

    ENDIF

    ! ----------------------
    ! close and exit
    ! ----------------------

    ! close the file and free the unit number
    CALL ftclos(unit, status)

    ! check for any error, and if so print out error messages
    IF (status .GT. 0) CALL printerror(status)

    RETURN

  END SUBROUTINE write_bintabh_KLOAD
  ! ==============================================================================

