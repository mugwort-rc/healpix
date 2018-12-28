 //=======================================================================;
  public void int[] query_disc ( int nside, vector0, double radius, boolean nest, boolean inclusive) {
    //=======================================================================;
    //;
    //      query_disc (Nside, Vector0, Radius, Listpix, Nlist[, Nest, Inclusive]);
    //      ----------;
    //      routine for pixel query in the RING or NESTED scheme;
    //      all pixels within an angular distance Radius of the center;
    //;
    //     Nside    = resolution parameter (a power of 2);
    //     Vector0  = central point vector position (x,y,z in double precision);
    //     Radius   = angular radius in RADIAN (in double precision);
    //     Listpix  = list of pixel closer to the center (angular distance) than Radius;
    //     Nlist    = number of pixels in the list;
    //     nest  (OPT), :0 by default, the output list is in RING scheme;
    //                  if set to 1, the output list is in NESTED scheme;
    //     inclusive (OPT) , :0 by default, only the pixels whose center;
    //                       lie in the triangle are listed on output;
    //                  if set to 1, all pixels overlapping the triangle are output;
    //;
    //      * all pixel numbers are in {0, 12*Nside*Nside - 1}
    //     NB : the dimension of the listpix array is fixed in the ing;
    //     routine and should be large enough for the specific configuration;
    //;
    //      lower level public void s ed by getdisc_ring :;
    //       (you don't need to know them);
    //      ring_num (nside, ir);
    //      --------;
    //      in_ring(nside, iz, phi0, dphi, listir, nir, nest=nest);
    //      -------;
    //;
    // v1.0, EH, TAC, ??;
    // v1.1, EH, Caltech, Dec-2001;
    //=======================================================================;
;
    int  irmin, irmax, ilist, iz, ip, nir, npix;
    double  norm_vect0;
    double  x0, y0, z0, radius_eff;
    double  a, b, c, cosang;
    double  dth1, dth2;
    double  phi0, cosphi0, cosdphi, dphi;
    double  rlat0, rlat1, rlat2, zmin, zmax, z;
    int  listir;
    int  status;
    character(len=*), parameter :: code = "QUERY_DISC";
    int  list_size, nlost;
    boolean;
;
    //=======================================================================;
;
    list_size = size(listpix);
    //     ---------- check inputs ----------------;
    npix = 12 * nside * nside;
;
    if (radius < 0.0_dp || radius > PI){
       //write(unit=*,fmt="(a)") code//"> the angular radius is in RADIAN ";
       //write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] ";
       throw new Exception("> program abort ");
    }
;
    do_inclusive = .false.;
    if (present(inclusive)){
       if (inclusive == 1) do_inclusive = .true.;
    }
;
    //     --------- allocate memory -------------;
    ALLOCATE( listir(0: 4*nside-1), STAT = status);
    if (status /= 0){
       //write(unit=*,fmt="(a)") code//"> can not allocate memory for listir :";
       throw new Exception("> program abort ");
    }
;
    dth1 = 1.0_dp / (3.0_dp*real(nside,kind=dp)**2);
    dth2 = 2.0_dp / (3.0_dp*real(nside,kind=dp));
;
    radius_eff = radius;
    if (do_inclusive){
       // increase radius by half pixel size;
       radius_eff = radius + PI / (4.0_dp*nside);
    }
    cosang = Math.cos(radius_eff);
;
    //     ---------- circle center -------------;
    norm_vect0 =  Math.sqrt(DOT_PRODUCT(vector0,vector0));
    x0 = vector0(1) / norm_vect0;
    y0 = vector0(2) / norm_vect0;
    z0 = vector0(3) / norm_vect0;
;
    phi0=0.0_dp;
    if ((x0/=0.0_dp)||(y0/=0.0_dp)) phi0 = Math.atan2 (y0, x0)  // in ]-Pi, Pi];
    cosphi0 = Math.cos(phi0);
    a = x0*x0 + y0*y0;
;
    //     --- coordinate z of highest and lowest points in the disc ---;
    rlat0  = Math.asin(z0)    // latitude in RAD of the center;
    rlat1  = rlat0 + radius_eff;
    rlat2  = rlat0 - radius_eff;
    if (rlat1 >=  halfpi){
       zmax =  1.0_dp;
    } else {
       zmax = Math.sin(rlat1);
    }
    irmin = ring_num(nside, zmax);
    irmin = MAX(1, irmin - 1) // start from a higher point, to be safe;
;
    if (rlat2 <= -halfpi){
       zmin = -1.0_dp;
    } else {
       zmin = Math.sin(rlat2);
    }
    irmax = ring_num(nside, zmin);
    irmax = MIN(4*nside-1, irmax + 1) // go down to a lower point;
;
    ilist = -1;
;
    //     ------------- loop on ring number ---------------------;
    do iz = irmin, irmax;
;
       if (iz <= nside-1){      // north polar cap
          z = 1.0_dp  - real(iz,kind=dp)**2 * dth1;
       } else { if (iz <= 3*nside){    // tropical band + equat.
          z = real(2*nside-iz,kind=dp) * dth2;
       } else {
          z = - 1.0_dp + real(4*nside-iz,kind=dp)**2 * dth1;
       }
;
       //        --------- phi range in the disc for each z ---------;
       b = cosang - z*z0;
       c = 1.0_dp - z*z;
       if ((x0==0.0_dp).and.(y0==0.0_dp)){
          cosdphi=-1.0_dp;
          dphi=PI;
          goto 500;
       }
       cosdphi = b / Math.sqrt(a*c);
       if (Math.abs(cosdphi) <= 1.0_dp){
          dphi = Math.acos (cosdphi) // in [0,Pi];
       } else {
          if (cosphi0 < cosdphi) goto 1000 // out of the disc;
          dphi = PI // all the pixels at this elevation are in the disc;
       }
500    continue;
;
       //        ------- finds pixels in the disc ---------;
        in_ring(nside, iz, phi0, dphi, listir, nir, nest);
;
       //        ----------- merge pixel lists -----------;
       nlost = ilist + nir + 1 - list_size;
       if ( nlost > 0 ){
          //print*,code//"> listpix is too short, it will be truncated at ",nir;
          //print*,"                         pixels lost : ", nlost;
          nir = nir - nlost;
       }
       do ip = 0, nir-1;
          ilist = ilist + 1;
          listpix(ilist) = listir(ip);
       }//enddo
;
1000   continue;
    }//enddo
;
    //     ------ total number of pixel in the disc --------;
    nlist = ilist + 1;
;
;
    //     ------- deallocate memory and exit ------;
    DEALLOCATE(listir);
;
    return;
  end public void  query_disc;
