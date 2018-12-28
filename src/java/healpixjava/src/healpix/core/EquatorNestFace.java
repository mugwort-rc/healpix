// Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.core/EquatorNestFace.java

package healpix.core;


public class EquatorNestFace extends NestedFace {
	
	/**
	 */
	public EquatorNestFace(int facenum, NestedMap map) {
		super(facenum,map);
    }
	
	/**
	 */
	public NestedFace[] nCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),2,3);		
    }
	
	/**
	 */
	public NestedFace[] sCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),0,5);		
    }
	
	/**
	 */
	public NestedFace[] wCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),0,1,2);		
    }
	
	/**
	 */
	public NestedFace[] eCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),3,4,5);		
    }
	
	/**
	 */
	public NestedFace nFaceNeigh() throws Exception {
		return adjacentFaces()[3];
    }
	
	/**
	 */
	public NestedFace sFaceNeigh() throws Exception {
		return adjacentFaces()[0];
    }
	
	/**
	 */
	public NestedFace wFaceNeigh() throws Exception {
		return adjacentFaces()[2];
    }
	
	/**
	 */
	public NestedFace eFaceNeigh() throws Exception {
		return adjacentFaces()[5];
    }
	
	/**
	 */
	protected void findAdjacentFaces() throws Exception {
		adjacent[0] =  map.getFace(8 + ((face-1)%4));
		adjacent[1] =  map.getFace(4 + ((face -1)%4));
		adjacent[2] =  map.getFace((face-1)%4);
		adjacent[3] =  map.getFace(face%4);
		adjacent[4] =  map.getFace(4 + ((face +1)%4));
		adjacent[5] =  map.getFace(8 + (face %4));
    }
	
	/**
	 */
	public int[] pixWcorner() throws Exception {
        NestedFace faces[] = wCornerNeigh();
        int ind =0,faceind=0;
        int pix[] = new int[8];
        // south side
        pix[ind++]=faces[faceind].xy2pix(nside-1,nside-2);
        pix[ind++]=faces[faceind++].xy2pix(nside-1,nside-1);

		// west tip
        pix[ind++] = faces[faceind++].xy2pix(nside-1,0);
        // north side
        pix[ind++]=faces[faceind].xy2pix(0,0);
        pix[ind++]=faces[faceind].xy2pix(1,0);
        // now the internal pixels
        int pixnum= xy2pix(0,nside-1);
        pix[ind++] = pixnum+1;
        pix[ind++] = pixnum-1;
        pix[ind++] = pixnum-2;

        return pix;
    }
	
	/**
	 */
	public int[] pixEcorner() throws Exception {
        NestedFace faces[] = eCornerNeigh();
        int ind =0,faceind=0;
        int pix[] = new int[8];
        pix[ind++]=faces[2].xy2pix(nside-2,nside-1);
        // first the internal pixels
        int pixnum= xy2pix(nside-1,0);
        pix[ind++] = pixnum-1;
        pix[ind++] = pixnum+1;
        pix[ind++] = pixnum+2;
        // north side
        pix[ind++]=faces[0].xy2pix(0,1);
        pix[ind++]=faces[0].xy2pix(0,0);
        // east side
        pix[ind++]=faces[1].xy2pix(0,nside-1);
		// back to the se side
        pix[ind++]=faces[2].xy2pix(nside-1,nside-1);
        return pix;
    }
	
	/**
	 */
	public NestedFace sTipNeigh() throws Exception {
		throw new Exception(" No south tip neighbour for equitorial faces");
    }
	
	/**
	 */
	public NestedFace wTipNeigh() throws Exception {
		return adjacentFaces()[1];
    }
	
	/**
	 */
	public NestedFace eTipNeigh() throws Exception {
		return adjacentFaces()[4];
    }
	
	/**
	 */
	public NestedFace nTipNeigh() throws Exception {
		throw new Exception(" No north tip neighbour for equitorial faces");
    }
	
	/**
	   No north corner here so fill with -1..
	 */
	protected void fillBoxNCorner(int[][] theBox, int xd, int yd) throws Exception {
		int bx = theBox.length -xd ;
		fillBoxVal(theBox,-1,bx,0,xd,yd);
	}
	
	/**
	   No south corner on the equator fill with -1.
	 */
	protected void fillBoxSCorner(int[][] theBox, int xd, int yd) throws Exception {
		int by =theBox[0].length -yd;
		fillBoxVal(theBox,-1,0,by,xd,yd);
	}
}
