package healpix.core;

import healpix.tools.*;
import java.awt.Point;

/**
   A class representing a face in Healpix - this is abstract faces are either Northern, Southern or Equatorial. This represents the Nested Scheme.
 */
public abstract class NestedFace {
	protected int face;
	protected int nside;
	protected int npface;
	
	/**
	   offset number of pixels before this faces e.g. (face*(nside*nside))-1 this has to be added to any internal pixel number to find its actual pixel number.
	 */
	protected int offset;
	protected healpix.core.NestedFace adjacent[];
	NestedMap map;
	
	/**
	   Constructor - can not be called directly use getFace
	 */
	protected NestedFace(int facenum, NestedMap map) {
		this.face = facenum;
		this.nside = map.nside();
		this.npface = map.nopixface();
		this.offset = (npface*facenum);
		//if (offset > 0) offset--;
		this.map = map;
	
    }
	
	/**
	   Faces are constructed according to number, this get method is a factory for the correct type of face
	 */
	public static NestedFace getFace(int facenum, NestedMap imap) throws Exception {
		NestedFace newface = null;
		if (facenum < 0 || facenum > 11) 
			throw new Exception ("Face out of range"+facenum);	
		if (facenum < 4) { // North Polar region
			newface = new NorthNestFace(facenum,imap);
		}
		if (facenum > 7) { // South Polar region
			newface = new SouthNestFace(facenum,imap);
		}
		if ( facenum <8  && facenum > 3) { //Equatorial region
			newface = new EquatorNestFace(facenum,imap);
		}
		
		return newface;
	
    }
	
	/**
	   Faces which are adjacent to the given face -  There are 6. They are returned in clockwise order starting south
	   In each face the first pixel is in the lowest corner of the diamond
	   the faces are                    (x,y) coordinate on each face
	   .   .   .   .   <--- North Pole
	   / \ / \ / \ / \                          ^        ^
	   . 0 . 1 . 2 . 3 . <--- z = 2/3             \      /
	   \ / \ / \ / \ /                        y   \    /  x
	   4 . 5 . 6 . 7 . 4 <--- equator               \  /
	   / \ / \ / \ / \                              \/
	   . 8 . 9 .10 .11 . <--- z = -2/3              (0,0) : lowest corner
	   \ / \ / \ / \ /
	   .   .   .   .   <--- South Pole
	   This reurns adjacent if it is already filled otherwise it used the abstract method findAdjacentFaces to have it filed first.
	 */
	public NestedFace[] adjacentFaces() throws Exception {
		if (adjacent == null) {
			adjacent = new NestedFace[6];
			findAdjacentFaces();
		}
		return adjacent;
    }
	
	/**
	   Faces which are adjacent to the given point in the given face - this is none if its in side a face but you should not be calling this in that case
	 */
	public NestedFace[] adjacentFaces(Point p) throws Exception {
		if (!(p.x==0 || p.y==0|| p.x==nside-1 || p.y==nside-1)) throw new Exception(" Internal point ");

		// south corner
		if (p.x==0 && p.y==0)	{
			return sCornerNeigh();
		}
		// north corner
		if (p.x==nside-1 && p.y==nside-1)	{
			return nCornerNeigh();
		}
		// east corner
		if (p.x==nside-1 && p.y==0)	{
			return eCornerNeigh();
		}
		// west  corner
		if (p.x==nside-1 && p.y==0)	{
			return wCornerNeigh();
		}

		// sides 
		NestedFace[] sideFace = new NestedFace[1];
		if (p.y == 0) sideFace[0]= eFaceNeigh();
		if (p.y == nside-1) sideFace[0] = wFaceNeigh();
		if (p.x == nside-1) sideFace[0] = nFaceNeigh();
		if (p.x == 0) sideFace[0] = sFaceNeigh();

		return sideFace;
    }
	
	/**
	   Pixels which are adjacent to the given point inside this face returns an array of pixel numbers
	 */
	public int[] internalNeigh(Point p) throws Exception {
		if (p.x==0 || p.y==0|| p.x==nside-1 || p.y==nside-1) 
			throw new Exception(" Not Internal point "+p);
		int[] neigh = new int[8];
		neigh[0] = xy2pix(p.x-1,p.y-1);
		neigh[1] = xy2pix(p.x-1,p.y);
		neigh[2] = xy2pix(p.x-1,p.y+1);
		neigh[3] = xy2pix(p.x,p.y+1);
		neigh[4] = xy2pix(p.x+1,p.y+1);
		neigh[5] = xy2pix(p.x+1,p.y);
		neigh[6] = xy2pix(p.x+1,p.y-1);
		neigh[7] = xy2pix(p.x,p.y-1);
		return neigh;
	
    }
	
	/**
	   Pixels which are adjacent to the given point in this face returns an array of pixel numbers
	 */
	public int[] neighbours(Point p) throws Exception {
		// easy case first
		if (!(p.x==0 || p.y==0|| p.x==nside-1 || p.y==nside-1))
			return internalNeigh(p);
		// now we have a boundary pixel to deal with 
		// Corners
		if (p.x==0 && p.y==0) return pixScorner();
		if (p.x==0 && p.y==nside-1) return pixWcorner();
		if (p.x==nside-1 && p.y==nside-1) return pixNcorner();
		if (p.x==nside-1 && p.y==0) return pixEcorner();

		// sides
		if (p.x==0) return pixSside(p.y); // guess you couls say SW etc.
		if (p.y==0)  return pixEside(p.x);
		if (p.x==nside-1) return pixNside(p.y);
		if (p.y==nside-1)  return pixWside(p.x);
		throw new Exception("Should Not Happen: neaighbours "+p);
    }
	
	/**
	   Conveinence function to select items from an array
	 */
	protected static NestedFace[] assemble(NestedFace[] inp, int i1, int i2) {
		NestedFace[] retarray = new NestedFace[2];
		retarray[0] = inp[i1];
		retarray[1] = inp[i2];
		return retarray;
	
    }
	
	/**
	   Conveinence function to select items from an array
	 */
	protected static NestedFace[] assemble(NestedFace[] inp, int i1, int i2, int i3) {
		NestedFace[] retarray = new NestedFace[3];
		retarray[0] = inp[i1];
		retarray[1] = inp[i2];
		retarray[2] = inp[i3];
		return retarray;
	
    }
	
	/**
	 */
	public abstract NestedFace[] nCornerNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace[] sCornerNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace[] wCornerNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace[] eCornerNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace nFaceNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace sFaceNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace wFaceNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace eFaceNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace sTipNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace wTipNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace eTipNeigh() throws Exception;
	
	/**
	 */
	public abstract NestedFace nTipNeigh() throws Exception;
	
	/**
	   actually find the faces if - will be called from adjacentfaces as required, this should populate adjacent.
	 */
	protected abstract void findAdjacentFaces() throws Exception;
	
	/**
	 */
	public int faceNum() {
		return face;
    }
	
	/**
	   Convert an x and y to a pix number within the face
	 */
	public int xy2pix(Point p) throws Exception {
		return( xy2pix(p.x,p.y));
    }
	
	/**
	   Convert a pix number to x,y within the face
	 */
	public Point pix2xy(int pix) throws Exception {
		int tpix = (pix >= npface) ? (pix -offset) : pix;
		return(Healpix.pix2xy_nest(tpix));
    }
	
	/**
	   Convert an x and y to a pix number within the face
	 */
	public int xy2pix(int x, int y) throws Exception {
		int pix = Healpix.xy2pix_nest(x,y);
		return (pix+offset);
    }
	
	/**
	 */
	public int[] neighbours(int pix) throws Exception {
		return neighbours(pix2xy(pix));
    }
	
	/**
	   Neighbours of the pixel in the south corner
	 */
	public int[] pixScorner() throws Exception {
		NestedFace faces[] = sCornerNeigh();
		int ind =0,faceind=0;;
		int pix[] = null;

		if (faces.length==3) { // eight neighbours
			pix = new int[8];
			pix[ind++] = faces[faceind++].xy2pix(nside-1,nside-1);
			
		} else {
			pix = new int[7];
		}
		// west side
		pix[ind++]=faces[faceind].xy2pix(nside-1,0); 
		pix[ind++]=faces[faceind++].xy2pix(nside-1,1);

		// now the internal pixels
		pix[ind++] = offset+2; 
		pix[ind++] = offset+3;
		pix[ind++] = offset+1;
		
		// east side
		pix[ind++]=faces[faceind].xy2pix(1,nside-1);
		pix[ind++]=faces[faceind].xy2pix(0,nside-1);
		return pix;
    }
	
	/**
	   Neighbours of the pixel in the west corner
	 */
	public abstract int[] pixWcorner() throws Exception;
	
	/**
	   Neighbours of the pixel in the north corner
	 */
	public int[] pixNcorner() throws Exception {
		NestedFace faces[] = nCornerNeigh();
		int ind =0,faceind=0;;
		int pix[] = null;

		if (faces.length==3) { // eight neighbours
			pix = new int[8];
		} else {
			pix = new int[7];
		}
		// south internal 
		pix[ind++] = offset+npface-4; 
		pix[ind++] = offset+npface-2; 
		// west side
		pix[ind++]=faces[faceind].xy2pix(nside-2,0); 
		pix[ind++]=faces[faceind++].xy2pix(nside-1,0);

		if (faces.length==3) { // eight neighbours
			//north if appropriate
			pix[ind++]=faces[faceind++].xy2pix(0,0);
		}
		// east side
		pix[ind++]=faces[faceind].xy2pix(0,nside-1);
		pix[ind++]=faces[faceind].xy2pix(0,nside-2);
		// internal east pix
		pix[ind++] = offset+npface-3; 
		return pix;
    }
	
	/**
	   Neighbours of the pixel in the east corner
	 */
	public abstract int[] pixEcorner() throws Exception;
	
	/**
	   Pixels borderinng the pixel on the south(west) face.
	 */
	public int[] pixSside(int y) throws Exception {
		if (y==0 || y == nside-1) throw new Exception("Corner Not side pix");
		int[] neigh = new int[8];
		NestedFace sf = sFaceNeigh();
		neigh[0] = sf.xy2pix(nside-1,y-1);
		neigh[1] = sf.xy2pix(nside-1,y);
		neigh[2] = sf.xy2pix(nside-1,y+1);
		// internal pixels
		neigh[3] = xy2pix(0,y+1);
		neigh[4] = xy2pix(0+1,y+1);
		neigh[5] = xy2pix(0+1,y);
		neigh[6] = xy2pix(0+1,y-1);
		neigh[7] = xy2pix(0,y-1);
		return neigh;
    }
	
	/**
	   Pixels borderinng the pixel on the nortth(east )face.
	 */
	public int[] pixNside(int y) throws Exception {
		if (y==0 || y == nside-1) throw new Exception("Corner Not side pix");
		int[] neigh = new int[8];
		// internal
		neigh[0] = xy2pix(nside-2,y-1);
		neigh[1] = xy2pix(nside-2,y);
		neigh[2] = xy2pix(nside-2,y+1);
		neigh[3] = xy2pix(nside-1,y+1);
		NestedFace sf = nFaceNeigh();
		neigh[4] = sf.xy2pix(0,y+1);
		neigh[5] = sf.xy2pix(0,y);
		neigh[6] = sf.xy2pix(0,y-1);
		// internal pixels
		neigh[7] = xy2pix(nside-1,y-1);
		return neigh;
    }
	
	/**
	   Pixels borderinng the pixel on the (north)west face.
	 */
	public int[] pixWside(int x) throws Exception {
		if (x==0 || x == nside-1) throw new Exception("Corner Not side pix");
		int[] neigh = new int[8];
		int ind=0;
		// internal
		neigh[ind++] = xy2pix(x-1,nside-2);
		neigh[ind++] = xy2pix(x-1,nside-1);
		NestedFace f = wFaceNeigh();
		neigh[ind++] = f.xy2pix(x-1,0);
		neigh[ind++] = f.xy2pix(x,0);
		neigh[ind++] = f.xy2pix(x+1,0);
		// internal pixels
		neigh[ind++] = xy2pix(x+1,nside-1);
		neigh[ind++] = xy2pix(x+1,nside-2);
		neigh[ind++] = xy2pix(x,nside-2);
		return neigh;
    }
	
	/**
	   Pixels borderinng the pixel on the (south)east face.
	 */
	public int[] pixEside(int x) throws Exception {
		if (x==0 || x == nside-1) throw new Exception("Corner Not side pix");
		int[] neigh = new int[8];
		NestedFace f = eFaceNeigh();
		int ind=0;
		neigh[ind++] = f.xy2pix(x-1,nside-1);
		// internal
		neigh[ind++] = xy2pix(x-1,0);
		neigh[ind++] = xy2pix(x-1,1);
		neigh[ind++] = xy2pix(x,1);
		neigh[ind++] = xy2pix(x+1,1);
		neigh[ind++] = xy2pix(x+1,0);
		// other eside external
		neigh[ind++] = f.xy2pix(x+1,nside-1);
		neigh[ind++] = f.xy2pix(x,nside-1);
		return neigh;
    }
	
	/**
	 */

	
	
	
	/**
	   Put givin value in a square of the 2'd array starting at xmin,yim for  xd,yd elements.
	 */
	public void fillBoxVal(int[][] theBox, int val, int xmin, int ymin, int xd, int yd) {
		for (int x=0; x< xd; x++) {
			for (int y=0; y < yd; y++ ){
				theBox[x+xmin][y+ymin] = val;
			}
		}
	}
	
	/**
	   Get pixelnumbers of the neighours of a pixel in this face.
	 */
	public int[][] box(int pix, int nspix, int wepix) throws Exception {
		// allocate result 2-d array
		int xdim = nspix*2 +1;	
		int ydim = wepix*2 +1;	
		int[][] theBox = new int[xdim][ydim]; 
		// given pixnum is at center (zero offset).
		theBox[nspix][wepix] = pix;
		Point pixp = pix2xy(pix);
		// get internal part of the face and calculate how far we need
		// to go in other directions outside this face.
		int lowx = pixp.x - nspix;
		int sw = 0, se=0, ne=0,nw=0;
		if (lowx< 0 ) {
			sw= Math.abs(lowx); // number pixels we need to go south west  
			lowx=0;
		}
		int highx = pixp.x + nspix;
		if (highx >= nside ) {
			ne= highx- (nside -1); // number pixels  north east
			highx=nside-1;
		}
		int highy = pixp.y + wepix;
		if (highy >= nside ) {
			nw= highy- (nside -1); // number pixels north west
			highy=nside-1;
		}
		int lowy = pixp.y - wepix;
		if (lowy< 0 ) {
			se= Math.abs(lowy); // number pixels south east
			lowy=0;
		}

		// set numbers in theBox
		// the array offset in y is upside down compared to healpix
		// so we need to invert it here 	
		int ypos = theBox[0].length - se -1 ;
		for (int x=0; x<= highx-lowx; x++) {
			for (int y=0; y <= highy-lowy; y++ ){
				// this xy2pix returns the full pixelnumer including its offset
				theBox[x+sw][ypos-y] = xy2pix(lowx+x,lowy+y);
			}
		}

		// now we have to start looking in adjacent faces as necessary ...
		// going clockwise starting at the bottom of the face.
		int xpos;
		if (sw>0 && se>0 ) { //got a corner to deal with
			fillBoxSCorner(theBox,sw,se);
		}
		if (sw>0) { // need to go south west
			fillBoxSW(theBox,lowy,highy,sw,ypos);
		}
		if (sw>0 && nw>0 ) { //got a corner to deal with
			fillBoxWCorner(theBox,sw,nw);
		}
		if (nw>0) { // need to go south west
			fillBoxNW(theBox,lowx,highx,sw,nw);
		}
		if (nw>0 && ne>0 ) { //got a corner to deal with
			fillBoxNCorner(theBox,ne,nw);
		}
		if (ne>0) { // need to go north east
		    xpos= (highx-lowx) +sw + 1;
			fillBoxNE(theBox,lowy,highy,xpos,ypos);
		}
		if (se>0 && ne>0 ) { //got a corner to deal with
			fillBoxECorner(theBox,ne,se);
		}
		if (se>0) { // need to go south east
		    xpos= (highx-lowx) + 1;
			fillBoxSE(theBox,lowx,highx,sw,se);
		}

		return theBox;
    }
	
	/**
	   Take a chunk of pixel numebrs from a face on the North East. ymin and ymax are the Healpix Y values to look for.  Size of the box gives the limit on X, bx and by give the location in the Box to fill in from . For north faces this is
	   filled with -1.
	 */
	protected void fillBoxNE(int[][] theBox, int ymin, int ymax, int bx, int by) throws Exception {
		// get the north east neighbour 
		NestedFace f = nFaceNeigh();
		for (int x=0; x< theBox.length -bx; x++) {
			for (int y=0; y <= ymax-ymin; y++ ){
				theBox[x+bx][by-y] = f.xy2pix(x,ymin+y);
			}
		}
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on the South East.ximn and xmax are the Healpix X values to look for.  Size of the box gives the limit on Y, bx and by give the location in the Box to fill in from.
	   SouthNestFace fills this with -1 by overridding this method.
	 */
	protected void fillBoxSE(int[][] theBox, int xmin, int xmax, int bx, int by) throws Exception {
		// get the south east neighbour 
		NestedFace f = eFaceNeigh();
		int ypos=theBox[0].length -by;
		for (int x=0; x<= xmax-xmin; x++) {
			for (int y=0; y < by; y++ ){
				theBox[x+bx][ypos+y] = f.xy2pix(xmin+x,(nside -1)-y);
			}
		}
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on the North West .ximn and xmax are the Healpix X values to look for.  Size of the box gives the limit on Y, bx and by give the location in the Box to fill in from.
	   NorthNestFace fills this with -1 by overridding this method.
	 */
	protected void fillBoxNW(int[][] theBox, int xmin, int xmax, int bx, int by) throws Exception {
		// get the NW neighbour 
		NestedFace f = wFaceNeigh();
		for (int x=0; x<= xmax-xmin; x++) {
			for (int y=0; y < by; y++ ){
				theBox[x+bx][y] = f.xy2pix(xmin+x,(by-1)-y);

			}
		}
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on the South West . ymin and ymax are the Healpix Y values to look for.   bx limits x by gives the y offset in the Box to fill in from . For south faces this is filled with -1 (SouthNestFace overrides this method)..
	 */
	protected void fillBoxSW(int[][] theBox, int ymin, int ymax, int bx, int by) throws Exception {
		// get the south west neighbour 
		NestedFace f = sFaceNeigh();
		int lowx = nside-bx;
		for (int x=0; x<bx; x++) {
			for (int y=0; y <= ymax-ymin; y++ ){
				theBox[x][by-y] = f.xy2pix(lowx+x,ymin+y);
			}
		}
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on East corner. xd and yd are the Healpix X an Y number of pixels to go into the corner Face. The rest can be worked out from theBox dimensions .
	 */
	protected void fillBoxECorner(int[][] theBox, int xd, int yd) throws Exception {
		int bx = theBox.length - xd ;
		int by = theBox[0].length -yd;
		// get the east corner 
		NestedFace f = eTipNeigh();
		for (int x=0; x< xd; x++) {
			for (int y=0; y < yd; y++ ){
				theBox[x+bx][by+y] = f.xy2pix(x,(nside -1)-y);
			}
		}
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on West corner. xd and yd are the Healpix X an Y number of pixels to go into the corner Face. The rest can be worked out from theBox dimensions .
	 */
	protected void fillBoxWCorner(int[][] theBox, int xd, int yd) throws Exception {
		// get the west corner 
		NestedFace f = wTipNeigh();
		for (int x=0; x< xd; x++) {
			for (int y=0; y < yd; y++ ){
				theBox[(xd-1)-x][(yd-1)-y] = f.xy2pix((nside -1)-x,y);
			}
		}
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on North corner. xd and yd are the Healpix X an Y number of pixels to go into the corner Face. The rest can be worked out from theBox dimensions .
	 */
	protected void fillBoxNCorner(int[][] theBox, int xd, int yd) throws Exception {
		int by = yd -1;
		int bx = theBox.length -xd ;
		// get the north corner 
		NestedFace f = nTipNeigh();
		for (int x=0; x< xd; x++) {
			for (int y=0; y < yd; y++ ){
				theBox[bx+x][by-y] = f.xy2pix(x,y);
			}
		}
	}
	
	/**
	   Take a chunk of pixel numbers from a face on South corner. xd and yd are the Healpix X an Y number of pixels to go into the corner Face. The rest can be worked out from theBox dimensions .
	 */
	protected void fillBoxSCorner(int[][] theBox, int xd, int yd) throws Exception {
		int by = theBox[0].length -yd;
		// get the south corner 
		int xpos = nside - xd;
		NestedFace f = sTipNeigh();
		for (int x=0; x< xd; x++) {
			for (int y=0; y < yd; y++ ){
				theBox[x][by+y] = f.xy2pix(xpos+x,(nside -1)-y);
			}
		}
	}

}
