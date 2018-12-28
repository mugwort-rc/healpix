package healpix.core.test;

import java.util.*;
import java.text.DecimalFormat;
import healpix.core.*;
import junit.framework.*;

public class TestPixBox extends TestCase
{

  public  void testPixBox()
  {
	  int nside = 8;
	  int pix = 10;
	  int ind=0;
	  int xr= 3 , yr=3;
		System.out.println("Getting "+xr+","+yr+" pixels around "+pix);

		try {
			NestedMap nm = new NestedMap(nside);
			int[][] neigh = nm.box(pix,xr,yr);
			output(neigh);
		} catch (Exception e) {
			System.err.println("Cant do it");
			e.printStackTrace();
		}
		  System.out.println("Finished  ::" +new Date());
            }

	 void output(int[][] ar)  {
			DecimalFormat form = new DecimalFormat("000");
			for (int y=0;y< ar[0].length;y++){
				for (int x=0;x<ar.length;x++) {
					if (ar[x][y]<0) 
					   System.out.print("  -1");
					else
					   System.out.print(" "+form.format(ar[x][y]));
				}
				System.out.println();
			}
	}

}

