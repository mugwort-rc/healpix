 //=======================================================================;
  function ring_num (nside, z) result(ring_num_result);
    //=======================================================================;
    //     returns the ring number in {1, 4*nside-1}
    //     from the z coordinate;
    //=======================================================================;
    int  ring_num_result;
    double  z;
    int  nside;

    int  iring;
    //=======================================================================;

    //     ----- equatorial regime ---------;
    iring = NINT( nside*(2.0_dp-1.500_dp*z));

    //     ----- north cap ------;
    if (z > twothird){
       iring = NINT( nside* Math.sqrt(3.0_dp*(1.0_dp-z)));
       if (iring == 0) iring = 1;
    }

    //     ----- south cap -----;
    if (z < -twothird   ){
       iring = NINT( nside* Math.sqrt(3.0_dp*(1.0_dp+z)));
       if (iring == 0) iring = 1;
       iring = 4*nside - iring;
    }

    ring_num_result = iring;

    return;
  end function ring_num;
