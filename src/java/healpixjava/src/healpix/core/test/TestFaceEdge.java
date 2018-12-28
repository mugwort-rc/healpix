package healpix.core.test;

import java.util.*;
import healpix.tools.*;
import healpix.core.*;
import junit.framework.*;

public class TestFaceEdge extends TestCase{

    public void testEdges() throws Exception{
        doEdges(8,11);
    }
  public  void doEdges(int nside,int facen)throws Exception{
      
      NestedMap nm = new NestedMap(nside);
      
      NestedFace nf10 = nm.getFace(facen);
      int cpix4 = nf10.xy2pix((nside-1)/2, 0);

      for(int x = 0; x < nside; x++){
		int pix  = nf10.xy2pix(x, 0);
		int pixa  = nf10.xy2pix(0, x);
		System.out.println("("+x+",0)"+pix+" (0,"+x+")"+pixa);
      }
      System.out.println();
    
  }
  
  static final void output(double[][] box){
    for(int x = 0; x < box.length; x++){
	for(int y = 0; y < box[0].length; y++){
	  System.out.print(box[x][y]+" ");
	}
	System.out.println();
      }
  }

   static final void output(int[][] box){
     for(int x = 0; x < box.length; x++){
	for(int y = 0; y < box[0].length; y++){
	  System.out.print(box[x][y]+" ");
	}
	System.out.println();
      }
   }
}
    
  
    
    
