// Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.core/NorthNestFace.java

package healpix.core;


/**
   A Norther Face of a nest jhealpix map.
 */
public class NorthNestFace extends PolarNF {
	
	/**
	   Constructor - can not be called directly use getFace
	 */
	public NorthNestFace(int facenum, NestedMap map) {
		super(facenum,map);
    }
	
	/**
	 */
	protected void findAdjacentFaces() throws Exception {
		adjacent[0] =  map.getFace(face +8);
		adjacent[1] =  map.getFace(face +4);
		adjacent[2] =  map.getFace(((face -1 ) +4)%4);
		adjacent[3] =  map.getFace((face +2)%4);
		adjacent[4] =  map.getFace((face +1)%4);
		adjacent[5] =  map.getFace(4 + ((face+1)%4));
    }
	
	/**
	   Neighbours of the pixel in the north corner
	 */
	public int[] pixNcorner() throws Exception {
		NestedFace faces[] = nCornerNeigh();
        int ind =0,faceind=0;;
        int pix[] = new int[8];
        // south internal
        pix[ind++] = offset+npface-4;
        pix[ind++] = offset+npface-2;
        // west side
        pix[ind++]=faces[faceind].xy2pix(nside-1,nside-2);
        pix[ind++]=faces[faceind++].xy2pix(nside-1,nside-1);
        pix[ind++]=faces[faceind++].xy2pix(nside-1,nside-1);

        // east side
        pix[ind++]=faces[faceind].xy2pix(nside-1,nside-1);
        pix[ind++]=faces[faceind].xy2pix(nside-2,nside-1);
        // internal east pix
        pix[ind++] = offset+npface-3;
        return pix;
    }
	
	/**
	 */
	public int[] pixWcorner() throws Exception {
        NestedFace faces[] = wCornerNeigh();
        int ind =0,faceind=0;
        int pix[] = new int[7];
        // south side
        pix[ind++]=faces[faceind].xy2pix(nside-1,nside-2);
        pix[ind++]=faces[faceind++].xy2pix(nside-1,nside-1);
        // north side
        pix[ind++]=faces[faceind].xy2pix(nside-1,0);
        pix[ind++]=faces[faceind].xy2pix(nside-1,1);
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
        int pix[] = new int[7];
        // south side
        pix[ind++]=faces[1].xy2pix(nside-2,nside-1);
        // now the internal pixels
        int pixnum= xy2pix(nside-1,0);
        pix[ind++] = pixnum-1;
        pix[ind++] = pixnum+1;
        pix[ind++] = pixnum+2;
        // east side
        pix[ind++]=faces[0].xy2pix(1,nside-1);
        pix[ind++]=faces[0].xy2pix(0,nside-1);
        pix[ind++]=faces[1].xy2pix(nside-1,nside-1);
        return pix;
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
        neigh[4] = sf.xy2pix(y+1,nside-1);
        neigh[5] = sf.xy2pix(y,nside-1);
        neigh[6] = sf.xy2pix(y-1,nside-1);
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
        NestedFace sf = wFaceNeigh();
        neigh[ind++] = sf.xy2pix(nside-1,x-1);
        neigh[ind++] = sf.xy2pix(nside-1,x);
        neigh[ind++] = sf.xy2pix(nside-1,x+1);
        // internal pixels
        neigh[ind++] = xy2pix(x+1,nside-1);
        neigh[ind++] = xy2pix(x+1,nside-2);
        neigh[ind++] = xy2pix(x,nside-2);
		return neigh;
    }
	
	/**
	   Skew the easttip neighbour
	 */
	protected void fillBoxNE(int[][] theBox, int ymin, int ymax, int bx, int by) throws Exception {
		int xd=theBox.length -bx;
		int yd = ymax - ymin +1;
		int xmax = theBox.length-1;
		int xpos = nside -xd ;

        // get the north west face 
        NestedFace f = nFaceNeigh();
        for (int x=0; x< xd; x++) {
           for (int y=0; y < yd; y++ ){

              theBox[xmax-x][by-y] = f.xy2pix(ymin+y,xpos+x);

           }
        }
	}
	
	/**
	   Skew west tip neighbour
	 */
	protected void fillBoxNW(int[][] theBox, int xmin, int xmax, int bx, int by) throws Exception {
		int xd = xmax-xmin +1;
		int yd = by;
		int ypos = nside -yd;
		//int xpos = nside -xd ;
		int xpos = xmin;
	//System.out.println(" NW xpos:"+xpos+" ypos:"+ypos+" xmin:"+xmin+" xmax:"+xmax+" xd:"+xd+" yd:"+yd+" bx:"+bx+" by:"+by);
        // get the north west face 
        NestedFace f = wFaceNeigh();
        for (int x=0; x< xd; x++) {
           for (int y=0; y < yd; y++ ){
              theBox[x+bx][y] = f.xy2pix(ypos+y,xpos+x);
//	System.out.println(" NW x:"+(x+bx)+" y:"+(y)+" pixx:"+(ypos+y)+" pixy:"+(xpos+x+" val:"+theBox[x+bx][y]));
           }
        }
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on North corner. xd and yd are the Healpix X an Y number of pixels to go into the corner Face. The rest can be worked out from theBox dimensions . The north neighbour here is turned around.
	 */
	protected void fillBoxNCorner(int[][] theBox, int xd, int yd) throws Exception {
        int bx = theBox.length -1 ;
		int xpos = nside -xd;
		int ypos = nside -yd;
        // get the north corner 
		
        NestedFace f = nTipNeigh();
        for (int x=0; x< xd; x++) {
           for (int y=0; y < yd; y++ ){
              theBox[bx-x][y] = f.xy2pix(xpos+x,ypos+y);
           }
        }

	}
}
