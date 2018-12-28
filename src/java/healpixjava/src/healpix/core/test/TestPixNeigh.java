package healpix.core.test;

import java.util.*;
import healpix.core.*;
import junit.framework.*;

public class TestPixNeigh extends TestCase
{

  public void testPixNeigh(){
      doPixNeigh(8,100);
  }
  
  public void doPixNeigh(int nside,int pix)
  {
	
		System.out.println("Using nside="+nside);

		try {
			NestedMap nm = new NestedMap(nside);
			if(pix < 0) {
			for (int f=0;f<12;f++) {
				doFace(nm.getFace(f),f*nm.nopixface(), nm.nopixface());
			}
			} else {
				int facen = pix/(nside*nside);
				NestedFace face = nm.getFace(facen);
				int[] neigh =face.neighbours(pix);
				System.out.print ("Neighbours of "+pix+"in face:"+face.faceNum()+" pix2xy"+ face.pix2xy(pix) +" are ");
				for (int j = 0; j < neigh.length; j++) {
					 System.out.print (" "+neigh[j]);
				}
				System.out.println();
				
			}
		} catch (Exception e) {
			System.err.println("Cant do it");
			e.printStackTrace();
		}
			  System.out.println("Finished  ::" +new Date());
  }

	public  void doFace(NestedFace face,int start,int nopix) throws Exception{
			for (int i=start;i<start+nopix;i++) {
				try {
				int[] neigh = face.neighbours(i);
				System.out.print ("Neighbours of "+i+"in face:"+face.faceNum()+" pix2xy"+ face.pix2xy(i) +" are ");
				for (int j = 0; j < neigh.length; j++) {
					 System.out.print (" "+neigh[j]);
				}
				System.out.println();
				} catch (Exception ex) {
					System.out.println("Skipping "+i+" "+ex);
					ex.printStackTrace();
				}
			}
	}

}

