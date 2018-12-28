package healpix.core.test;

import java.util.*;
import healpix.core.*;
import junit.framework.*;

public class TestFaceCorn extends TestCase
{

  public  void testFaceCorn() throws Exception
  {
	  int nside = 4;
	  int face = 0;

	NestedMap nm = new NestedMap(nside);
	for (int i=0;i<12;i++) doFace(i,nm);
	System.out.println("Finished  ::" +new Date());
  }

	 void show(NestedFace faces[], String comment) {
		System.out.print(comment+ " :");
		for (int j = 0; j < faces.length; j++) {
			 System.out.print (" "+faces[j].faceNum());
		}
		System.out.println();
	}
	 void show(NestedFace face, String comment) throws Exception{
		System.out.print(comment+ " :");
		 System.out.print (" "+face.faceNum());
		System.out.println();
	}

	 void doFace(int face, NestedMap nm) throws Exception{
	
		NestedFace theFace = nm.getFace(face);
		NestedFace[] faces= theFace.sCornerNeigh();
		show(faces,"South Corner "+face)	;
		faces= theFace.wCornerNeigh();
		show(faces,"West Corner  "+face)	;
		faces= theFace.nCornerNeigh();
		show(faces,"North Corner "+face)	;
		faces= theFace.eCornerNeigh();
		show(faces,"East Corner  "+face)	;

		NestedFace nf;
		nf= theFace.sFaceNeigh();
		show(nf,"South Side   "+face)	;
		nf= theFace.wFaceNeigh();
		show(nf,"West Side    "+face)	;
		nf= theFace.nFaceNeigh();
		show(nf,"North Side    "+face)	;
		nf= theFace.eFaceNeigh();
		show(nf,"East Side     "+face)	;
	
	}
}

