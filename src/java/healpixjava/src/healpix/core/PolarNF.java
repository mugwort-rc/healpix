// Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.core/PolarNF.java

package healpix.core;


/**
   Abstraction of Polar nested face - must be either north or south
 */
public abstract class PolarNF extends NestedFace {
	
	/**
	   Constructor - can not be called directly use getFace
	 */
	public PolarNF(int facenum, NestedMap map) {
		super(facenum,map);
    }
	
	/**
	 */
	public NestedFace nFaceNeigh() throws Exception {
		return adjacentFaces()[4];
    }
	
	/**
	 */
	public NestedFace sFaceNeigh() throws Exception {
		return adjacentFaces()[1];
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
	public NestedFace[] nCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),2,3,4);
    }
	
	/**
	 */
	public NestedFace[] sCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),0,1,5);
    }
	
	/**
	 */
	public NestedFace[] wCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),1,2);
    }
	
	/**
	 */
	public NestedFace[] eCornerNeigh() throws Exception {
		return assemble(adjacentFaces(),4,5);
    }
	
	/**
	 */
	public NestedFace sTipNeigh() throws Exception {
		return adjacentFaces()[0];
	}
	
	/**
	 */
	public NestedFace wTipNeigh() throws Exception {
		return adjacentFaces()[2];
	}
	
	/**
	 */
	public NestedFace eTipNeigh() throws Exception {
		return adjacentFaces()[4];
	}
	
	/**
	 */
	public NestedFace nTipNeigh() throws Exception {
		return adjacentFaces()[3];
	}
	
	/**
	   Rip is now in east west for North and south faces so we fill with blancks here.
	 */
	protected void fillBoxECorner(int[][] theBox, int xd, int yd) throws Exception {
		int bx = theBox.length - xd ;
		int by = theBox[0].length -yd;
		fillBoxVal(theBox,-1,bx,by,xd,yd);
	}
	
	/**
	   blanks go in here (ripping of faces North and south)
	 */
	protected void fillBoxWCorner(int[][] theBox, int xd, int yd) throws Exception {
		fillBoxVal(theBox,-1,0,0,xd,yd);
	}
}
