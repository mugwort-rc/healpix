package healpix.core.test;

import java.util.*;
import java.awt.*;
import healpix.core.*;
import junit.framework.*;

public class TestPixXY extends TestCase
{

  public static void testPixXY8() throws Exception 
  {
	  int nside = 8;

			NestedMap nm = new NestedMap(nside);
			int nopix = nm.nopixface();
			for (int f=0;f<12;f++) {
				System.out.println("Face: "+f);
			   doFace( nm.getFace(f), (f*nopix), nopix);
                        }
  }

	static void doFace(NestedFace face, int start, int end) throws Exception {
			for (int i=start;i<start+end;i++) {
				try {
					Point p = face.pix2xy(i);
					int pix = face.xy2pix(p.x,p.y);
                                        assertEquals("xy incorrect",i,pix);
				} catch (Exception ex) {
					System.out.println("Skipping "+i+" "+ex);
				}
			}
	}

}

