// Source file: /usr4/users/ggiardin/PlanckDev/planck/healpixj/java/planck/dmci/nonpersist/AngularPositionImp.java

package healpix.core;

import java.text.*;



/**
   An angular posiont theta phi
 */
public class AngularPosition {
    protected double theta = 0;
    protected double phi = 0;
    
    /**
       Default constructor 
       @roseuid 37C4EA240780
     */
    public AngularPosition() {
    }
    
    /**
       Simple constructor init both values
       @roseuid 37C4EA2443F8
     */
    public AngularPosition(double theta, double phi) {
    	
		this.theta=theta;
		this.phi=phi;
    }
    
    /**
       @roseuid 37C4EA24E691
     */
    public double theta() {
    return theta;
    }
    
    /**
       @roseuid 37C4EA24E7A8
     */
    public double phi() {
    return phi;
    }
    
    /**
       @roseuid 37C4EA24E7A9
     */
    public void setTheta(double val) {
      this.theta=val;
    }
    
    /**
       @roseuid 37C4EA24E7AB
     */
    public void setPhi(double val) {
		this.phi=val;
	
    }
    
    /**
       @roseuid 37C4EA24E7AD
     */
    public String toString() {
		DecimalFormat form = new DecimalFormat(" 0.000");
		return ("theta:"+form.format(theta)+" phi:"+form.format(phi));
	
    }
    
    /**
       @roseuid 3810539264F0
     */
    public void init(double t, double phi) {
    }
}
/*
AngularPositionImp.setTheata(double){
		this.theta=val;
	
    }

*/
