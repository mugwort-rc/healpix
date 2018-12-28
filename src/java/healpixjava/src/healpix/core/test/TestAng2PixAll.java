package healpix.core.test;

import java.util.*;
import java.text.*;
import junit.framework.*;
import healpix.tools.*;
import healpix.core.*;

public class TestAng2PixAll extends TestCase {

  public void testAllRing128() throws Exception{
      testAll(128,true);
  }
  public void testAllNest128() throws Exception{
      testAll(128,false);
  }
   public  void testAll(int nside, boolean ring) throws Exception{

    int ind=0;
      int length  = 12*nside*nside;
      DecimalFormat form = new DecimalFormat("#.###");
      int ns2 = nside*nside;
      int ecount = 0;
      for(int i = 0; i < length; i++){
	  AngularPosition pos =null;
	  if (ring) 
	  	pos = Healpix.pix2ang_ring(nside,i);
	  else
	  	pos = Healpix.pix2ang_nest(nside,i);
	 
	  int pix = 0;
	  if (ring) 
	  	pix = Healpix.ang2pix_ring(nside,pos.theta(), pos.phi());
	  else
	  	pix = Healpix.ang2pix_nest(nside,pos.theta(), pos.phi());
	  assertEquals(" theta "+form.format(Math.cos(pos.theta()))+" phi/pi "+form.format(pos.phi()/Healpix.pi),
		i,pix);
      }
      
    }

  }


