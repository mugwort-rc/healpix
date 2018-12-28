// Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.core/FacePoint.java

package healpix.core;

import java.awt.Point;


public class FacePoint extends Point {
    public int face;
    
    /**
     */
    public int getFace() {
		return face;
    }
    
    /**
     */
    public FacePoint(int x, int y, int face) {
		super(x,y);
		this.face = face;
    }
}
