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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computes the spherical harmonics coefficients (a_lm)
! and the power spectrum from a full sky map in HEALPIX pixelisation
! Written and developed by E. Hivon (efh@ipac.caltech.edu) and K. Gorski
! (krzysztof.m.gorski@jpl.nasa.gov) based on HEALPIX pixelisation developed
! by K. Gorski
!
! Copyright 1997 by Eric Hivon and Krzysztof M. Gorski.
!  All rights reserved.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=======================================================================
  !=======================================================================
  !  EXTERNAL LIBRARY:
  !     this code uses the CFITSIO library that can be found at
  !     http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
  !
  !  RELATED LITTERATURE:
  !     about HEALPIX   : see Gorski et al, 2005, ApJ, 622, 759
  !     about this code : see Hivon & Gorski, in preparation
  !
  !  HISTORY:
  !     April-October 1997, Eric Hivon, TAC (f77)
  !               Jan 1998 : translation in Fortran 90
  !               Aug 1998 : version 1.0.0
  !               Apr 1999 : precomp. plms, higher order iterations and alm dumping, FKH
  !               Aug 2000 : CMBFAST convention for polarisation
  !               Sep 2002 : implement new parser
  !
  !  FEEDBACK:
  !     for any questions : efh@ipac.caltech.edu
  !
  !=======================================================================
  !     version 2.0.0
  !=======================================================================
  ! this file can not be compiled on its own.
  ! It must be inserted into the file synfast.f90 by the command  include
  !
  integer(I4B),             parameter :: KALMC = KIND((1.0_KMAP, 1.0_KMAP)) ! precision of alm arrays
  integer(I4B),             parameter :: KCL   = KMAP  ! precision of c(l) arrays
  !-------------------------------------------------------------------------

  integer(I4B) :: nsmax, nlmax, nmmax

  complex(kind=KALMC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  complex(kind=KALMC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC2
  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_TQU
  real(kind=KMAP),     DIMENSION(:,:),   ALLOCATABLE :: map_TQU2
  real(kind=KMAP),     DIMENSION(:,:), target, ALLOCATABLE :: mask, mask_tmp
  real(kind=KMAP),     dimension(:), pointer :: pmask, pmask_tmp
  real(kind=KCL),      DIMENSION(:,:),   ALLOCATABLE :: clout
  real(kind=KCL),      DIMENSION(:),     ALLOCATABLE :: fact_norm
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: w8ring_TQU
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: dw8
  real(kind=DP),       DIMENSION(:,:),   ALLOCATABLE :: plm
  integer(kind=I4B) :: status

  integer(kind=I4B) :: l, mlpol, polar_fits, polcconv
  integer(kind=I4B) :: n_plm, plm_nside, plm_lmax, plm_pol, plm_mmax
  integer(kind=I4B) :: npixtot, nmaps
  integer(kind=I4B) :: ordering, order_type
  integer(kind=I4B) :: i, j
  integer(kind=I4B) :: iter, iter_order, polar, n_pols
  integer(kind=I4B) :: simul_type
  integer(kind=I4B) :: n_rings, n_rings_read, won, nmw8
  integer(kind=i4b) :: type

  real(kind=KMAP), dimension(1:3) :: mean, rms
  real(kind=KMAP)               :: fmissval, fmiss_mask
  real(kind=DP)                 :: pix_size_arcmin
  real(kind=DP)                 :: theta_cut_deg, cos_theta_cut
  real(kind=DP)                 :: nullval
  real(kind=DP), dimension(1:2) :: zbounds
  real(kind=DP)                 :: fsky
  real(kind=SP)                 :: clock_time, time0, time1
  real(kind=SP)                 :: ptime, ptime0, ptime1
  real(kind=DP)                 :: min_mask, max_mask

  character(len=80), DIMENSION(1:120) :: header
  integer :: ncl, nlheader

  character(len=FILENAMELEN)          :: infile
  character(len=FILENAMELEN)          :: outfile
!  character(len=FILENAMELEN)          :: parafile = ''
  character(len=FILENAMELEN)          :: outfile_alms
  character(len=FILENAMELEN)          :: infile_plm
  character(len=FILENAMELEN)          :: infile_w8
  character(len=FILENAMELEN)          :: w8name
  character(len=FILENAMELEN)          :: def_dir, def_file
  character(len=FILENAMELEN)          :: usr_dir, usr_file
  character(len=FILENAMELEN)          :: final_file
  character(len=FILENAMELEN)          :: healpixdir
  character(len=FILENAMELEN)          :: description
  character(len=FILENAMELEN)          :: maskfile
  character(len=100)                  :: chline, chline1
  character(len=20)                   :: coordsys
!  character(len=5)                    :: sstr
  LOGICAL(kind=LGT) :: bad, ok, polarisation, bcoupling, do_mask
!   character(len=*), PARAMETER :: code = "ANAFAST"
  character(len=*), PARAMETER :: version = "2.0.0"
  character(len=80), dimension(1:3)   :: units_map, units_pow

  type(paramfile_handle) :: handle

  real(kind=DP), dimension(0:3) :: mono_dip
  real(kind=DP)                 :: theta_dip, phi_dip, long_dip, lat_dip, ampl_dip
  integer(kind=I4B)             :: lowlreg

  integer(kind=I4B)        :: npix_mask, order_mask, nside_mask, nmasks, i_mask
  logical(LGT)             :: pessim_mask
  !-----------------------------------------------------------------------
  !                    get input parameters and create arrays
  !-----------------------------------------------------------------------
  call wall_clock_time(time0)
  call cpu_time(ptime0)
  !     --- read parameters interactively if no command-line arguments
  !     --- are given, otherwise interpret first command-line argument
  !     --- as parameter file name and read parameters from it:

  !     --- announces program, default banner ***
  PRINT *, " "
  PRINT *, "                "//code//" "//version
  write(*,'(a)') " *** Analysis and power spectrum calculation for a Temperature/Polarisation map ***"
  if (KMAP == SP) print*,'                Single precision outputs'
  if (KMAP == DP) print*,'                Double precision outputs'
  PRINT *, " "
  handle = parse_init(parafile)

  !     --- choose temp. only or  temp. + pol. ---
  description = concatnl( &
       & " Do you want to analyse", &
       & " 1) a Temperature only map", &
       & " 2) Temperature + polarisation maps")
  simul_type = parse_int(handle, 'simul_type', default=1, vmin=1, vmax=2, descr=description)
  if (simul_type .eq. 1) polarisation = .false.
  if (simul_type .eq. 2) polarisation = .true.
  bcoupling = polarisation ! <<<<
!!  print*, " Bcoupling :  ",bcoupling

  !     --- gets the file name for the map ---
  chline = "''"
  call getEnvironment("HEALPIX",healpixdir)
  if (trim(healpixdir)/='') chline = trim(healpixdir)//'/test/map.fits'
  description = concatnl( &
       & "", &
       & " Enter input file name (Map FITS file): ")
  infile = parse_string(handle, 'infile',default=chline, descr=description, filestatus='old')

  !     --- finds out the pixel number of the map and its ordering and coordinates ---
  npixtot = getsize_fits(infile, nmaps = nmaps, ordering=ordering, nside=nsmax,&    
       &              mlpol=mlpol, type = type, polarisation = polar_fits, &
       &             coordsys=coordsys, polcconv=polcconv)
  if (nsmax<=0) then
     print*,"Keyword NSIDE not found in FITS header!"
     call fatal_error(code)
  endif
  if (type == 3) npixtot = nside2npix(nsmax) ! cut sky input data set
  if (nsmax/=npix2nside(npixtot)) then
     print*,"FITS header keyword NSIDE does not correspond"
     print*,"to the size of the map!"
     call fatal_error(code)
  endif

  if (polarisation .and. (nmaps >=3) .and. polar_fits == -1) then
     print*,"The input fits file MAY NOT contain polarisation data."
     print*,"Proceed at your own risk"
  endif

  if (polarisation .and. (nmaps<3 .or. polar_fits ==0)) then
     print *,"The file does NOT contain polarisation maps"
     print *,"only the temperature field will be analyzed"
     polarisation = .false.
  endif

  if (polarisation .and. (polcconv == 2)) then
     print *,"The input map contains polarized data in the IAU coordinate convention (POLCCONV)"
     print *,code//" can not proceed with these data"
     print *,"See Healpix primer for details."
     call fatal_error(code)
  endif

  if (polarisation .and. (polcconv == 0)) then
     print *,"WARNING: the polarisation coordinate convention (POLCCONV) can not be determined"
     print *,"         COSMO will be assumed. See Healpix primer for details."
  endif

  !     --- check ordering scheme ---
  if ((ordering/=1).and.(ordering/=2)) then
     PRINT*,"The ordering scheme of the map must be RING or NESTED."
     PRINT*,"No ordering specification is given in the FITS-header!"
     call fatal_error(code)
  endif
  
  order_type = ordering

  !     --- gets the L range for the analysis ---
  WRITE(chline1,"(a,i5)") " The map has Nside = ",nsmax
  if (mlpol<=0) WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",3*nsmax-1,")"
  if (mlpol> 0) WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",mlpol,")"
  description = concatnl(&
       & "", &
       & chline1, &
       & " Enter the maximum l range (l_max) for the simulation. ", &
       & chline )
  nlmax = parse_int(handle, 'nlmax', default=2*nsmax, descr=description)

  ! ----------  default parameters ----------------
  nmmax   = nlmax
  npixtot = nside2npix(nsmax)
  pix_size_arcmin = 360.*60./SQRT(npixtot*PI)
  !   -------------------------------------

  !     --- ask for pixel mask ---
  description = concatnl(&
       & "", &
       & " Enter FITS file containing pixel mask and/or weight you want to use: ", &
       & " the map(s) will be multiplied by this mask(s) before analysis ", &
       & " (If '', all valid pixels are used with their current value) ")
  maskfile = parse_string(handle, 'maskfile', default="''", descr=description, filestatus='old')
  do_mask = (trim(maskfile) /= "")

  !     --- gets the cut applied to the map ---
  description = concatnl(&
       & "", &
       & " Enter the symmetric cut around the equator in DEGREES : ", &
       & " (One ignores data within |b| < b_cut)     0 <= b_cut = ")
  theta_cut_deg = parse_double(handle, 'theta_cut_deg', &
       &                       vmin=0.0_dp, default=0.0_dp, descr=description)
  cos_theta_cut = SIN(theta_cut_deg/180.d0*PI) !counted from equator instead of from pole

  zbounds = (/ cos_theta_cut , -cos_theta_cut /)
  if (theta_cut_deg<1e-4) zbounds = (/ -1.0_dp, 1.0_dp /) !keep all sphere if cut not set
  fsky = (zbounds(2)-zbounds(1))/2.0_dp
  if (fsky <= 0.0_dp) fsky = 1.0_dp + fsky
  write(*,"(a,f6.1,a)") "One keeps ",100.*fsky," % of the original map"

  !     --- ask about monopole/dipole removal ---
  description = concatnl(&
       & "", &
       & " To improve the temperature map analysis (specially on a cut sky) ", &
       & " you have the option of first regressing out the best fit monopole and dipole", &
       & "   NOTE : the regression is made on valid (unflagged) pixels, ", &
       & " out of the symmetric cut (if any)", &
       & "  Do you want to : ", &
       & " 0) Do the analysis on the raw map", &
       & " 1) Remove the monopole (l=0) first ", &
       & " 2) Remove the monopole and dipole (l=1) first")
  lowlreg = parse_int(handle, 'regression', vmin=0, vmax=2, default=0, descr=description)

  !          --- precomputed plms  ---
40 continue
  n_plm   = 0
  plm_pol = 0
  description = concatnl(&
       & "", &
       & " Enter name of file with precomputed p_lms :", &
       & "(eg, plm.fits) (if '' is typed, no file is read)" )
  infile_plm = parse_string(handle, 'plmfile', default="''", descr=description, filestatus='old')

  if (trim(infile_plm) /= "") then
     ! check that the plm file matches the current simulation (for nside,lmax,mmax)
     call read_par(infile_plm, plm_nside, plm_lmax, plm_pol, mmax=plm_mmax)

     if ((plm_pol /= 1).and.(.not.polarisation)) plm_pol=1

     if (plm_nside /= nsmax .or. plm_lmax /= nlmax .or. plm_mmax /= nmmax) then
        print *," "//code//"> Plm file does not match map to simulate in Nside, Lmax or Mmax ! "
        if (handle%interactive) goto 40
        call fatal_error(code)
     endif
     n_plm = (nmmax+1)*(2*nlmax+2-nmmax)*nsmax
     print *," "
  endif

50 continue
  !     --- gets the output power spectrum filename ---
  description = concatnl(&
       & "", &
       & " Enter Output Power Spectrum file name (eg, cl_out.fits) :", &
       & "  (or !cl_out.fits to overwrite an existing file)" )
  outfile = parse_string(handle, "outfile", &
       default="!cl_out.fits", descr=description, filestatus="new")

  !     --- gets the output alm-filename ---
  description = concatnl(&
       & "", &
       & " Enter output file name for alms(eg: alms.fits): ", &
       & " (or !alms.fits to overwrite an existing file): ", &
       & " (If '', the alms are not written to a file) ")
  outfile_alms = parse_string(handle, 'outfile_alms', default="''", descr=description, filestatus='new')

  if (trim(outfile)=="" .and. trim(outfile_alms)=="") then
     print*,'ERROR: No output file provided'
     if (handle%interactive) goto 50
     call fatal_error(code)
  endif
  !-----------------------------------------------------------------------
  !                      ask for weights
  !-----------------------------------------------------------------------

  description = concatnl(&
       & "", &
       & " Do you want to use ring weights in the analysis (1=yes; 2,0=no) ?")
  won = parse_int(handle, 'won', vmin=0, vmax=2, default=0, descr=description)

  infile_w8=""

  if (won.eq.1) then

     ! default weight file name
!     write (sstr,"(I5.5)") nsmax
!      w8name="weight_ring_n"//trim(sstr)//".fits"
     w8name="weight_ring_n"//trim(string(nsmax,"(i5.5)"))//".fits"

     def_file = trim(w8name)
     def_dir  = concatnl("","../data","./data","..")
     call getEnvironment("HEALPIX",healpixdir)
     if (trim(healpixdir) .ne. "") then
        def_dir = concatnl(&
             & def_dir, &
             & healpixdir, &
             & trim(healpixdir)//"data", &
             & trim(healpixdir)//"/data", &
             & trim(healpixdir)//char(92)//"data") !backslash
     endif

22   continue
     final_file = ''
     ok = .false.
     ! if interactive, try default name in default directories first
     if (handle%interactive) ok = scan_directories(def_dir, def_file, final_file)
     if (.not. ok) then
        ! not found, ask the user
        description = concatnl("",&
             &        " Could not find weight file", &
             &        " Enter the directory where this file can be found:")
        usr_dir = parse_string(handle,'w8filedir',default='',descr=description)
        if (trim(usr_dir) == '') usr_dir = trim(def_dir)
        description = concatnl("",&
             &        " Enter the name of the weight file:")
        usr_file = parse_string(handle,'w8file',default=def_file,descr=description)
        ! look for new name in user provided or default directories
        ok   = scan_directories(usr_dir, usr_file, final_file)
        ! if still fails, crash or ask again if interactive
        if (.not. ok) then
           print*,' File not found:'//trim(usr_file)
           if (handle%interactive) goto 22
           call fatal_error(code)
        endif
     endif
     infile_w8 = final_file
  endif

  !     --- gets the order of iteration ---
  description = concatnl(&
       & "", &
       & " Do you want : ", &
       & " 0) a standard analysis", &
       & " 1,2,3,4....) an iterative analysis", &
       & " (enter order of iteration, 3rd order is usually optimal)")
  iter_order=parse_int(handle, 'iter_order', vmin=0, default=0, descr=description)
  if (.not. handle%interactive) then
     select case (iter_order)
     case (0)
        PRINT*," Standard analysis"
     case default
        PRINT*," Iterative analysis"
     end select
  end if

  PRINT *," "
  call parse_summarize(handle,code=lcode,prec=KMAP)
  call brag_openmp()

  !-----------------------------------------------------------------------
  !              allocate space for arrays
  !-----------------------------------------------------------------------
  polar = 0
  if (polarisation) polar = 1
  n_pols = 1 + 2*polar ! either 1 or 3

  ALLOCATE(map_TQU(0:npixtot-1,1:n_pols),stat = status)
  call assert_alloc(status,code,"map_TQU")

  ALLOCATE(w8ring_TQU(1:2*nsmax,1:n_pols),stat = status)
  call assert_alloc(status,code,"w8ring_TQU")

  ALLOCATE(dw8(1:2*nsmax,1:n_pols),stat = status)
  call assert_alloc(status,code,"dw8")

  ALLOCATE(fact_norm(0:nlmax),stat = status)
  call assert_alloc(status,code,"fact_norm")

  if (n_plm/=0) then
     ALLOCATE(plm(0:n_plm-1,1:plm_pol),stat = status)
     call assert_alloc(status,code,"plm")
  endif
  !-----------------------------------------------------------------------
  !                      reads the map
  !-----------------------------------------------------------------------
  PRINT *,"      "//code//"> Inputting map "

  fmissval = 0.0
  if (lowlreg > 0) fmissval = real(HPX_DBADVAL,kind=KMAP)
  call input_map(infile, map_TQU(0:,1:n_pols), &
       &  npixtot, n_pols, fmissval=fmissval, units=units_map(1:n_pols))
  units_map(1) = adjustl(units_map(1))
  units_pow(1) = trim(units_map(1))//'^2'


  !-----------------------------------------------------------------------
  !                      reads the mask
  !-----------------------------------------------------------------------
  if (do_mask) then
     fmiss_mask = 0.0_KMAP ! set to mask=0 missing pixels
     pessim_mask = .true.
     PRINT *,"      "//code//"> Inputting mask from "//trim(maskfile)
     npix_mask = getsize_fits(maskfile, ordering=order_mask, nside=nside_mask, nmaps=nmasks)

     if (nside_mask <= 0) then
        print*,"Keyword NSIDE not found in FITS header!"
        call fatal_error(code)
     endif
     if ((order_mask /= 1).and.(order_mask /= 2)) then
        print*,"The ordering scheme of the mask must be RING or NESTED."
        print*,"No ordering specification is given in the FITS header!"
        call fatal_error(code)
     endif

     if (nside_mask == nsmax) then
        ! mask has same resolution as map: read and proceed
        allocate(mask(0:npix_mask-1,1:nmasks),stat=status)
        call assert_alloc(status,code,"mask")
        call input_map(maskfile, mask, npix_mask, nmasks, fmissval = fmiss_mask)
     else
        ! mask has different resolution than map: read
        allocate(mask_tmp(0:npix_mask-1,1:nmasks),stat=status)
        call assert_alloc(status,code,"mask")
        call input_map(maskfile, mask_tmp, npix_mask, nmasks, fmissval = fmiss_mask)
        ! then degrade or prograde
        print*,'Modifying mask to match map resolution:'
        write(*,'(a,i6,a,i6)') 'Nside:  ',nside_mask,' --> ',nsmax
        do i_mask = 1, nmasks
           pmask_tmp => mask_tmp(:,i_mask)
           pmask     => mask(:,i_mask)
           if (order_mask == 2) then ! nested scheme
              call udgrade_nest(pmask_tmp, nside_mask, pmask, nsmax, &
                   &            fmissval= fmiss_mask, pessimistic=pessim_mask)
           else ! ring scheme
              call udgrade_ring(pmask_tmp, nside_mask, pmask, nsmax, &
                   &            fmissval= fmiss_mask, pessimistic=pessim_mask)
           endif
        enddo
        deallocate(mask_tmp)
        nside_mask = nsmax ; npix_mask = npixtot
     endif

     !--------------------------------------
     ! **** what if nmasks /= nmaps ??? ****
     !---------------------------------------
     do i_mask = 1, nmasks
        max_mask = maxval(mask(:,i_mask))
        min_mask = minval(mask(:,i_mask))
        if (max_mask > 0. .and. min_mask >=0.) then
           write(*,'(a,i2,a)') 'Field #', i_mask,' is OK'
           exit
        endif
        write(*,'(a,i2,a,2(1pg11.3))') 'Field #', i_mask, &
             & ' of '//trim(maskfile)//' has range ', min_mask, max_mask
        print*,' not a valid mask, will try next field'
        if (i_mask == nmasks) then
           print*,'Error : No valid mask found in '//trim(maskfile)
           call fatal_error(code)
        endif
     enddo

     ! keep only first valid mask found
     pmask => mask(:,i_mask)
     nmasks = 1

     ! make sure mask and sky map have same ordering
     if (order_mask /= order_type) then 
        PRINT *,"      "//code//"> Reorder mask "      
        if (order_mask == 2) call convert_nest2ring (nside_mask, pmask)
        if (order_mask == 1) call convert_ring2nest (nside_mask, pmask)
        order_mask = order_type
     endif

  else ! create mask place holder
     allocate(pmask(0:1))
  endif 

  !-----------------------------------------------------------------------
  !                      remove dipole
  !-----------------------------------------------------------------------

  if (lowlreg > 0) then
     PRINT *,"      "//code//"> Remove monopole (and dipole) from Temperature map"
     call remove_dipole(nsmax, map_TQU(0:npixtot-1,1), ordering, lowlreg, &
          & mono_dip(0:), zbounds, fmissval=fmissval, mask=pmask)

     if (fmissval /= 0.0) then
        do j=1,n_pols
           do i=0, npixtot-1
!               if (abs(map_TQU(i,j)/fmissval-1.0) < 1.e-6*fmissval) then
              if (abs(map_TQU(i,j)-fmissval) <= abs(1.e-6*fmissval)) then
                 map_TQU(i,j) = 0.0_kmap
              endif
           enddo
        enddo
     endif

     write(unit=*,fmt=8000)  " Monopole = ", mono_dip(0)," "//trim(units_map(1))
     if (lowlreg > 1) then
        call vec2ang( mono_dip(1:3), theta_dip, phi_dip)
        ampl_dip = sqrt(sum(mono_dip(1:3)**2))
        long_dip =      phi_dip  /PI*180.0_dp
        lat_dip  = 90.0-theta_dip/PI*180.0_dp
        write(unit=*,fmt=8000)  " Dipole   = ",ampl_dip," "//trim(units_map(1))
        write(unit=*,fmt=8002) "(long.,lat.) = (",long_dip,lat_dip,") Deg"
     endif
8000 format(a,1pg13.3,a)
8002 format(a,1pg9.3,', ',1pg9.3,a)

  endif
  !-----------------------------------------------------------------------
  !                      Input weights
  !-----------------------------------------------------------------------

  dw8=0.0_dp ! read as DP, even in SP in FITS file
  if (trim(infile_w8)/="") then
     n_rings = 2 * nsmax
     n_rings_read = getsize_fits(infile_w8, nmaps=nmw8)

     if (n_rings_read /= n_rings) then
        PRINT *," "
        print *,"wrong ring weight file:"//trim(infile_w8)
        call fatal_error(code)
     endif
     nmw8 = min(nmw8, n_pols)

     PRINT *,"      "//code//"> Inputting Quadrature ring weights "
     call input_map(infile_w8, dw8, n_rings, nmw8, fmissval=0.0_dp)
  endif
  w8ring_TQU = 1.0_dp + dw8

  !-----------------------------------------------------------------------
  !                 convert to RING if necessary
  !-----------------------------------------------------------------------
  if (order_type == 2) then
     PRINT *,"      "//code//"> Convert Nest -> Ring "
     call convert_nest2ring (nsmax, map_TQU)
     if (do_mask) call convert_nest2ring (nsmax, pmask)
  endif

  !-----------------------------------------------------------------------
  !                   multiply MAP by MASK
  !-----------------------------------------------------------------------
  if (do_mask) then
     do i=1,n_pols
        map_TQU(:,i) = map_TQU(:,i) * pmask
     enddo
  endif

  !-----------------------------------------------------------------------
  !                    precomputed plms
  !-----------------------------------------------------------------------

  if (n_plm/=0) then
     if (plm_pol.eq.1) then
        PRINT*,"      "//code//"> Reading precomputed scalar P_lms "
     else
        PRINT*,"      "//code//"> Reading precomputed tensor P_lms "
     endif
     call read_dbintab(infile_plm,plm,n_plm,plm_pol,nullval,anynull=bad)
     if (bad) call fatal_error("Missing points in P_lm-file!")
  endif

  !-----------------------------------------------------------------------
  !                        map to alm
  !-----------------------------------------------------------------------

  ALLOCATE(alm_TGC(1:n_pols, 0:nlmax, 0:nmmax),stat = status)
  call assert_alloc(status,code,"alm_TGC")
  if (iter_order > 0) then
     ALLOCATE(alm_TGC2(1:n_pols, 0:nlmax, 0:nmmax),stat = status)
     call assert_alloc(status,code,"alm_TGC2")

     ALLOCATE(map_TQU2(0:npixtot-1,1:n_pols),stat = status)
     call assert_alloc(status,code,"map_TQU2")

     alm_TGC2 = 0.0_KALMC      ! final alm
     map_TQU2 = map_TQU        ! input map
     alm_TGC  = 0.0_KALMC      ! intermediate alm
  endif

  PRINT *,"      "//code//"> Analyse map "

  do iter = 0, iter_order
     ! perform a Jacobi iterative scheme
     ! map -> alm
     select case (plm_pol+polar*10)
     case(0) ! Temperature only
        call map2alm(nsmax, nlmax, nmmax, map_TQU(:,1), alm_TGC, zbounds, w8ring_TQU)
     case(1) ! T with precomputed Ylm
        call map2alm(nsmax, nlmax, nmmax, map_TQU(:,1), alm_TGC, zbounds, w8ring_TQU,plm(:,1))
     case(10) ! T+P
        call map2alm(nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring_TQU)
     case(11) ! T+P with precomputed Ylm_T
        call map2alm(nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring_TQU,plm(:,1))
     case(13) ! T+P with precomputed Ylm
        call map2alm(nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring_TQU,plm)
     end select

     if (iter_order == 0) exit
     if (iter == 0 .and. iter_order > 0) then
        PRINT *,"      "//code//"> Iterated map analysis "
        PRINT *,"       input map rms"
        do i=1,n_pols
           mean(i) = sum(map_TQU(:, i))/npixtot
           rms (i) = sqrt( sum((map_TQU(:,i)-mean(i))**2)/(npixtot-1) )
        enddo
        write(*,9111) 0,rms(1:n_pols)
        PRINT *,"       residual map rms"
     endif
9111 format(i3,3(g20.5))
     map_TQU = 0.0_KMAP
     ! alm_2 = alm_2 + alm : increment alms
     alm_TGC2 = alm_TGC2 + alm_TGC

     if (iter == iter_order) exit

     ! alm_2 -> map : reconstructed map
     select case (plm_pol+polar*10)
     case(0)
        call alm2map(nsmax, nlmax, nmmax, alm_TGC2, map_TQU(:,1))
     case(1)
        call alm2map(nsmax, nlmax, nmmax, alm_TGC2, map_TQU(:,1), plm(:,1))
     case(10)
        call alm2map(nsmax, nlmax, nmmax, alm_TGC2, map_TQU)
     case(11)
        call alm2map(nsmax, nlmax, nmmax, alm_TGC2, map_TQU, plm(:,1))
     case(13)
        call alm2map(nsmax, nlmax, nmmax, alm_TGC2, map_TQU, plm)
     end select

     ! map = map_2 - map : residual map
     map_TQU = map_TQU2 - map_TQU

     do i=1,n_pols
        mean(i) = sum(map_TQU(:, i))/npixtot
        rms (i) = sqrt( sum((map_TQU(:,i)-mean(i))**2)/(npixtot-1) )
     enddo
     write(*,9111) iter,rms(1:n_pols)
  enddo

  if (iter_order > 0) then
     alm_TGC = alm_TGC2 ! final alms
     map_TQU = map_TQU2 ! input map
     deallocate(alm_TGC2, map_TQU2)
  endif

  deallocate( map_TQU )
  deallocate(w8ring_TQU)
  deallocate(dw8)

  !-----------------------------------------------------------------------
  !                  power spectrum calculation
  !-----------------------------------------------------------------------

  PRINT *,"      "//code//"> Compute power spectrum "

  ncl = 1
  if (polarisation) ncl = 4
  if (bcoupling) ncl = 6
  allocate(clout(0:nlmax,1:ncl))

  call alm2cl(nlmax, nmmax, alm_TGC, alm_TGC, clout)
  !-----------------------------------------------------------------------
  !                  outputs the power spectrum
  !-----------------------------------------------------------------------

  if (trim(outfile)/="") then
     PRINT *,"      "//code//"> outputting Power Spectrum "

     do l = 0, nlmax
        fact_norm(l) = 1.0_KCL !  Normalisation=1, Output 'pure' C_ls
     enddo
     do i = 1, ncl
        clout(:,i) = clout(:,i) * fact_norm(:)
     enddo

     header = ""
     call add_card(header,"COMMENT","-----------------------------------------------")
     call add_card(header,"COMMENT","     Map Analysis Specific Keywords      ")
     call add_card(header,"COMMENT","-----------------------------------------------")
     call add_card(header,"EXTNAME","'ANALYSED POWER SPECTRUM'")
     call add_card(header,"COMMENT"," POWER SPECTRUM : C(l) ")
     call add_card(header)
     call add_card(header,"CREATOR",code,        "Software creating the FITS file")
     call add_card(header,"VERSION",version,     "Version of the simulation software")
     call add_card(header,"CUT-SKY",theta_cut_deg," [deg] Symmetric Cut Sky")
     call add_card(header,"NITERALM",iter_order," Number of iterations for a_lms extraction")
     call add_card(header,"POLAR",polarisation," Polarisation included (True/False)")
     call add_card(header,"BCROSS",bcoupling," Magnetic cross terms included (True/False)")
     call add_card(header,"MAX-LPOL",nlmax," Maximum L multipole order")
     call add_card(header,"LREGRESS",lowlreg," Regression of low multipoles (0/1/2)")
     call add_card(header)
     call add_card(header,"HISTORY","Input Map = "//TRIM(infile))
     call add_card(header,"NSIDE",nsmax," Resolution Parameter of Input Map")
     call add_card(header,"COORDSYS",coordsys," Input Map Coordinates")
     call add_card(header)
     call add_card(header,"TTYPE1", "TEMPERATURE","Temperature C(l)")
     units_pow(1) = ''''//trim(units_pow(1))//'''' ! add quotes to keep words together
     call add_card(header,"TUNIT1", units_pow(1),"square of map units")
     call add_card(header)
     !
     if (polarisation) then
        call add_card(header,"TTYPE2", "GRADIENT","Gradient (=ELECTRIC) polarisation C(l)")
        call add_card(header,"TUNIT2", units_pow(1),"power spectrum units")
        call add_card(header)
        !
        call add_card(header,"TTYPE3", "CURL","Curl (=MAGNETIC) polarisation C(l)")
        call add_card(header,"TUNIT3", units_pow(1),"power spectrum units")
        call add_card(header)
        !
        call add_card(header,"TTYPE4", "G-T","Gradient-Temperature (=CROSS) terms")
        call add_card(header,"TUNIT4", units_pow(1),"power spectrum units")
        call add_card(header)
        !
        if (bcoupling) then
           call add_card(header,"TTYPE5", "C-T","Curl-Temperature terms")
           call add_card(header,"TUNIT5", units_pow(1),"power spectrum units")
           call add_card(header)
           !
           call add_card(header,"TTYPE6", "C-G","Curl-Gradient terms")
           call add_card(header,"TUNIT6", units_pow(1),"power spectrum units")
           call add_card(header)
        endif
        !
        call add_card(header,"COMMENT","The polarisation power spectra have the same definition as in CMBFAST")
        call add_card(header)
     endif

     nlheader = SIZE(header)
     call write_asctab (clout, nlmax, ncl, header, nlheader, outfile)
  endif
  !-----------------------------------------------------------------------
  !                  outputs the alm
  !-----------------------------------------------------------------------
  if (trim(outfile_alms)/="") then

     PRINT *,"      "//code//"> outputting alms "

     ! write each component (T,G,C) in a different extension of the same file
     do i=1,1+polar*2
        header = ""
        call add_card(header,"COMMENT","-----------------------------------------------")
        call add_card(header,"COMMENT","     Map Analysis Specific Keywords      ")
        call add_card(header,"COMMENT","-----------------------------------------------")
        if (i == 1) then
           call add_card(header,"EXTNAME","""ANALYSED a_lms (TEMPERATURE)""")
        elseif (i == 2) then
           call add_card(header,"EXTNAME","'ANALYSED a_lms (GRAD / ELECTRIC component)'")
        elseif (i == 3) then
           call add_card(header,"EXTNAME","'ANALYSED a_lms (CURL / MAGNETIC component)'")
        endif
        call add_card(header)
        call add_card(header)
        call add_card(header,"CREATOR",code,        "Software creating the FITS file")
        call add_card(header,"VERSION",version,     "Version of the simulation software")
        call add_card(header,"CUT-SKY",theta_cut_deg," [deg] Symmetric Cut Sky")
        call add_card(header,"NITERALM",iter_order,"Number of iterations for a_lms extraction")
        call add_card(header,"POLAR",polarisation," Polarisation included (True/False)")
        call add_card(header,"LREGRESS",lowlreg," Regression of low multipoles (0/1/2)")
        call add_card(header)
        call add_card(header,"HISTORY","Input Map = "//TRIM(infile))
        call add_card(header,"NSIDE",nsmax," Resolution Parameter of Input Map")
        call add_card(header,"COORDSYS",coordsys," Input Map Coordinates")
        if (i == 1) then
           call add_card(header,"COMMENT","Temperature a_lms")
        else
           call add_card(header,"COMMENT","Polarisation a_lms")
        endif
        call add_card(header,"COMMENT"," The real and imaginary part of the a_lm with m>=0")
        call add_card(header)
        call add_card(header,"COMMENT","Input Map = "//TRIM(infile))
        call add_card(header)
        call add_card(header,"TTYPE1", "INDEX"," i = l^2 + l + m + 1")
        call add_card(header,"TUNIT1", "   "," index")
        call add_card(header)
        call add_card(header,"TTYPE2", "REAL"," REAL a_lm")
        units_map(1) = ''''//trim(units_map(1))//'''' ! add quotes to keep words together
        call add_card(header,"TUNIT2", units_map(1)," alm units")
        call add_card(header)
        !
        call add_card(header,"TTYPE3", "IMAG"," IMAGINARY a_lm")
        call add_card(header,"TUNIT3", units_map(1)," alm units")
        call add_card(header)
        !
        call add_card(header)
        call add_card(header)
        call add_card(header)
        !
        nlheader = SIZE(header)
        call dump_alms (outfile_alms,alm_TGC(i,0:nlmax,0:nlmax),nlmax,header,nlheader,i-1)
     enddo
  endif

  !-----------------------------------------------------------------------
  !                      deallocate memory for arrays
  !-----------------------------------------------------------------------

  DEALLOCATE( alm_TGC )

  !-----------------------------------------------------------------------
  !                      report card
  !-----------------------------------------------------------------------
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  write(*,9000) " "
  write(*,9000) " Report Card for "//code//" analysis run"
  write(*,9000) "----------------------------------------"
  write(*,9000) " "
  if (.not. polarisation) then
     write(*,9000) "      Temperature alone"
  else
     write(*,9000) "    Temperature + Polarisation"
  endif
  write(*,9000) " "
  write(*,9000) " Input map            : "//TRIM(infile)
  write(*,9005) " Multipole range      : 0 <= l <= ", nlmax
  write(*,9010) " Number of pixels     : ", npixtot
  write(*,9020) " Pixel size in arcmin : ", pix_size_arcmin
  write(*,9020) " Symmetric cut in deg : ", theta_cut_deg
  if (do_mask) write(*,9000) &
       &        " Mask read from       : "//trim(maskfile)
  write(*,9000) " Output power spectrum: "//TRIM(outfile)
  if (outfile_alms/="") write(*,9000) " Output alm           : "//TRIM(outfile_alms)
  if (won == 1) write(*,9000) " Ring weight file     : "//trim(infile_w8)
  write(*,9030) " Clock and CPU time [s] : ", clock_time, ptime

  !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------

  write(*,9000) " "
  write(*,9000) " "//code//"> normal completion"

9000 format(a)
9005 format(a,i8)
9010 format(a,i16)
9020 format(a,g20.5)
9030 format(a,f11.2,f11.2)

