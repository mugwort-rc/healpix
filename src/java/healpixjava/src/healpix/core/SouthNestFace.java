// Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.core/SouthNestFace.java

package healpix.core;


public class SouthNestFace extends PolarNF {
	
	/**
	   Constructor - can not be called directly use getFace
	 */
	public SouthNestFace(int facenum, NestedMap map) {
		super (facenum,map);
    }
	
	/**
	 */
	protected void findAdjacentFaces() throws Exception {
		adjacent[0] =  map.getFace(8 + ((face +2)%4));
		adjacent[1] =  map.getFace(8 + ((face +3) %4));
		adjacent[2] =  map.getFace(face -4);
		adjacent[3] =  map.getFace(face -8);
		adjacent[4] =  map.getFace(4 + ((face +1)%4));
		adjacent[5] =  map.getFace(8 + ((face +1)%4));

    }
	
	/**
	   Neighbours of the pixel in the south corner
	 */
	public int[] pixScorner() throws Exception {
        NestedFace faces[] = sCornerNeigh();
        int ind =0,faceind=0;;
        int pix[] = new int[8];

        pix[ind++] = faces[faceind++].xy2pix(0,0);

        // west side
        pix[ind++]=faces[faceind].xy2pix(0,0);
        pix[ind++]=faces[faceind++].xy2pix(1,0);

        int toff=offset;
        // now the internal pixels
        pix[ind++] = toff+2;
        pix[ind++] = toff+3;
        pix[ind++] = toff+1;

        // east side
        pix[ind++]=faces[faceind].xy2pix(0,1);
        pix[ind++]=faces[faceind].xy2pix(0,0);
        return pix;
    }
	
	/**
	 */
	public int[] pixWcorner() throws Exception {
        NestedFace faces[] = wCornerNeigh();
        int ind =0,faceind=0;
        int pix[] = new int[7];
        // south side
        pix[ind++]=faces[faceind].xy2pix(nside-2,0);
        pix[ind++]=faces[faceind++].xy2pix(nside-1,0);
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
        int pix[] = new int[7];
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
        pix[ind++]=faces[1].xy2pix(0,nside-2);
        return pix;
    }
	
	/**
	   Pixels borderinng the pixel on the south(west) face.
	 */
	public int[] pixSside(int y) throws Exception {
        if (y==0 || y == nside-1) throw new Exception("Corner Not side pix");
        int[] neigh = new int[8];
        NestedFace sf = sFaceNeigh();
        neigh[0] = sf.xy2pix(y-1,0);
        neigh[1] = sf.xy2pix(y,0);
        neigh[2] = sf.xy2pix(y+1,0);
        // internal pixels
        neigh[3] = xy2pix(0,y+1);
        neigh[4] = xy2pix(0+1,y+1);
        neigh[5] = xy2pix(0+1,y);
        neigh[6] = xy2pix(0+1,y-1);
        neigh[7] = xy2pix(0,y-1);
        return neigh;
    }
	
	/**
	   Pixels borderinng the pixel on the (south)east face.
	 */
	public int[] pixEside(int x) throws Exception {
        if (x==0 || x == nside-1) throw new Exception("Corner Not side pix");
        int[] neigh = new int[8];
        int ind=0;
        // internal
        neigh[ind++] = xy2pix(x-1,0);
        neigh[ind++] = xy2pix(x-1,1);
        neigh[ind++] = xy2pix(x,1);
        neigh[ind++] = xy2pix(x+1,1);
        neigh[ind++] = xy2pix(x+1,0);
        // other eside external
        NestedFace f = eFaceNeigh();
        neigh[ind++] = f.xy2pix(0,x+1);
        neigh[ind++] = f.xy2pix(0,x);
        neigh[ind++] = f.xy2pix(0,x-1);
		return neigh;
    }
	
	/**
	   Rotate the corner face and make it adjacent for box puposes
	 */
	protected void fillBoxSE(int[][] theBox, int xmin, int xmax, int bx, int by) throws Exception {
		int yoff=theBox[0].length -by;
        int xd = xmax-xmin +1;
		int ypos= nside - xd;
        // get the south east face
        NestedFace f = eFaceNeigh();
    //System.out.println(" SE xmin:"+xmin+" xmax:"+xmax+" ypos:"+ypos+" yoff:"+yoff+" xd:"+xd+" bx:"+bx+" by:"+by+" looking in face "+f.faceNum());
        for (int x=0; x< xd; x++) {
           for (int y=0; y < by; y++ ){
              theBox[x+bx][yoff+y] = f.xy2pix(y,xmin+x);
  // System.out.println(" SE x:"+(x+bx)+" y:"+(y+yoff)+" pixx:"+(y)+" pixy:"+(ypos+x)+" val:"+theBox[x+bx][yoff+y]);
           }
        }

	}
	
	/**
	   Fill in with the corner face
	 */
	protected void fillBoxSW(int[][] theBox, int ymin, int ymax, int bx, int by) throws Exception {
		int yd=ymax-ymin +1;
		int yoff = by-yd + 1;
        int xpos = bx -1 ;
        // get the north west face 
        NestedFace f = sFaceNeigh();
        //System.out.println(" SW xpos:"+xpos+" ymin:"+ymin+" ymax:"+ymax+" yd:"+yd+" bx:"+bx+" by:"+by+" lloking in face "+f.faceNum());
        for (int x=0; x< bx; x++) {
           for (int y=0; y < yd; y++ ){
              theBox[x][yoff+y] = f.xy2pix(ymax-y,xpos-x);
        //System.out.println(" SW x:"+(x)+" y:"+(yoff+y)+" pixx:"+(ymax-y)+" pixy:"+(xpos-x)+" val:"+theBox[x][y+yoff]);
           }
        }
	
	}
	
	/**
	 */
	public NestedFace wTipNeigh() throws Exception {
	     return adjacentFaces()[1];
	}
	
	/**
	 */
	public NestedFace eTipNeigh() throws Exception {
	     return adjacentFaces()[5];
	}
	
	/**
	   Take a chunk of pixel numebrs from a face on South corner. xd and yd are the Healpix X an Y number of pixels to go into the corner Face. The rest can be worked out from theBox dimensions .The south neighbour here is upside down.
	 */
	protected void fillBoxSCorner(int[][] theBox, int xd, int yd) throws Exception {
		int by = theBox[0].length -yd;
		// get the south corner 
		NestedFace f = sTipNeigh();
		for (int x=0; x< xd; x++) {
			for (int y=0; y < yd; y++ ){
				theBox[(xd-1)-x][by+y] = f.xy2pix(x,y);
			}
		}

	}
}
