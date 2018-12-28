package healpix.core.test;

import java.util.*;
import healpix.core.*;
import junit.framework.*;

public class TestFaceNeigh extends TestCase
{

  public static void testFaceNeigh() throws Exception
  {
	  int nside = 4;
	  System.out.println("Start  ::" +new Date());
          NestedMap nm = new NestedMap(nside);
	  for (int i=0;i<12;i++) {
		NestedFace faces[]= nm.getFace(i).adjacentFaces();
		System.out.print ("Neighbours of "+i+" are ");
		for (int j = 0; j < faces.length; j++) {
		 System.out.print (" "+faces[j].faceNum());
		}
                System.out.println();
          }
	  System.out.println("Finished  ::" +new Date());
  }

}

