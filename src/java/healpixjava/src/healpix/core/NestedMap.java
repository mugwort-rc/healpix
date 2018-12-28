// Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.core/NestedMap.java

package healpix.core;

import healpix.tools.*;
import java.awt.Point;

/**
   A class representing the geometry of the Helapix nestedMap just a smart container for faces so we do not keep constructing them.
 */
public class NestedMap {
	protected int nside;
	protected int nopix;
	protected int nopixface;
	protected NestedFace faces[] = new NestedFace [12] ;
	
	/**
	 */
	public NestedMap(int nside) {
		this.nside = nside;
		this.nopixface = (nside*nside);
		this.nopix = nopixface*12;
	
    }
	
	/**
	 */
	public NestedFace getFace(int face) throws Exception {
		if (faces[face] == null) {
			faces[face] = NestedFace.getFace(face,this);
		}
		return faces[face];
	
    }
	
	/**
	 */
	public int nopix() {
     return nopix;
    }
	
	/**
	 */
	public int nopixface() {
     return nopixface;
    }
	
	/**
	 */
	public int nside() {
     return nside;
    }
	
	
	/**
	   Return the numbers of the pixels in an area around the given pixel. The parameter nspix is the number of pixels to gather north and south of pixnum (e.g. like a raduis) likewise wepix is for west and east of pixnum. Hence we could get a rectangle.
	   The returned 2-d array will have pixnum at the center and will have dimensions nspix*2 +1 and wepix*2+1.
	   The box must be smaller than 2*nside for now - this seems fair as a box bigger than 2*nside would be covering a very large area and then one should be using a map.
	   Thge resulting box is in healpix oprientation e.g. rotated 45 degrees - it does not
	   seem sensible to do anything else.
	   Note that holes will appear in returned array because Healpix does not map on to a flat 2d structure. Such holes will be assigned -1 values.  These occur in all places where a pixel has 7 neighbours e.g. North and south of the Equatorial Faces. Here the faces are treated as if the sphere were cut along the face borders and laid out flat. North south wraping is taken into acount. The orientation of the North and South Polar faces change depending on whether the search is about a pixel in a Polar or Equitorial face - the faces are oriented to give most neighbours. this gives the best Gausian results.
	   Hence searchin for 2 pixels around pixel 21 of an NSIDE=4 Map (North Pole) would give a 5,5 array of pixel numbers with some -1 values where the east corner   neighbours are cut away because of flattening the sphere e.g. face 3 and 7 are ripped apart there. . Also you would need to rotate it 45degrees for Healpix orientation to visiualise it. SO it looks like:
	   025 028 029 046 044
	   019 022 023 043 041
	   017 020 021 042 040
	   107 110 111  -1  -1
	   105 108 109  -1  -1
	   
	   However due to the orirntation change searching aroun 111 in the Equatorial will produce
	   019 022 023  -1  -1
	   017 020 021  -1  -1
	   107 110 111 042 043
	   105 108 109 040 041
	   099 102 103 034 035
	   
	   Here (since we are in an eqitorial face) the rip is made between faces 1 and 2.
	   
	   
	   wil Jan 6 2000
	 */
	public int[][] box(int pix, int nspix, int wepix) throws Exception {
		// it seems fair to restrict to nspix< nside e.g. a box can not 
		// be biger than a face - this will make life easier !
		if (nspix >= nside || nspix >= nside) 
			throw new Exception("Box dimensions too large for given geometry");
		// which face is pix in.
	 	NestedFace face = getFace(pix/nopixface);  
		// and away we go
		return face.box(pix,nspix,wepix);
	}
	
	
	/**
	   Count blancks(-1) in theBox.
	 */
	public int countBlancks(int[][] theBox) throws Exception {
		int count=0;
		for (int x = 0; x < theBox.length;x++) {
			for (int y = 0; y < theBox[0].length;y++) {
				if (theBox[x][y] < 0) count++;
			}
		}
		return count;
	}
}
