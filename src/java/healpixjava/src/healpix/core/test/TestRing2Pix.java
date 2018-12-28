package healpix.core.test;

import java.util.*;
import healpix.tools.*;
import healpix.core.*;
import junit.framework.*;

public class TestRing2Pix extends TestCase
{

  public void testRing2Pix14t2() throws Exception{
      int pix = Ring2Pix(2,14);
      assertEquals("Wrong pix",23,pix);
      
  }
  public int Ring2Pix(int nside, int pix) throws Exception
  {
	int npix = Healpix.ring2nest(nside,pix);
	System.out.println("Nest:"+npix+" = Ring:"+pix + " Nside is"+nside);
	System.out.println("Finished  ::" +new Date());
        return npix;
  }

}

