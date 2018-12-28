package healpix.tools;

import healpix.core.AngularPosition;

public final class CoordTransform {

    protected static final double pi = 3.14159265;
    protected static final double twoPi=2.0*pi;
    protected static final double fourPi=4.0*pi;
    protected static final double degToRad=180.0/pi;

    protected static final double[][] psi = {
	{0.57595865315,4.92619181360,0.00000000000,0.00000000000,0.11129056012,4.70053728340},
	{0.57477043300,4.93682924650,0.00000000000,0.00000000000,0.11142137093,4.71279419371}
    };
    protected static final double[][] stheta = {
	{0.88781538514,-0.88781538514, 0.39788119938,-0.39788119938, 0.86766174755,-0.86766174755},
	{0.88998808748,-0.88998808748, 0.39777715593,-0.39777715593, 0.86766622025,-0.86766622025}
    };

    protected static final double[][] ctheta = {
	{0.46019978478,0.46019978478,0.91743694670,0.91743694670,0.49715499774,0.49715499774},
	{0.45598377618,0.45598377618,0.91748206207,0.91748206207,0.49714719172,0.49714719172}
    };

    protected static final double[][] phi = {
	{4.92619181360,0.57595865315,0.00000000000,0.00000000000,4.70053728340,0.11129056012},
	{4.93682924650,0.57477043300,0.00000000000,0.00000000000,4.71279419371,0.11142137093}
    };
    
    public static final int EQ2GAL = 0; //RA-Dec (2000) -> Galactic
    public static final int GAL2EQ = 1; //Galactic      -> RA-Dec
    public static final int EQ2ECL = 2; //RA-Dec        -> Ecliptic
    public static final int ECL2EQ = 3; //Ecliptic      -> RA-Dec
    public static final int ECL2GAL = 4; //Ecliptic      -> Galactic
    public static final int GAL2ECL = 5; //Galactic      -> Ecliptic  

    /**Transforms an angular position in radians in a given coordinate system to a position
       in an other coordinate system, also in radians. RA-Dec position are intended in 
       Equinox J2000
       */
    public static final AngularPosition transform(AngularPosition pos, int trType) 
	throws Exception {

	double ao,bo,a,b,sb,cb,cbsa;
	int J2000 = 1;
	//by setting J2000 = 0, RA-Dec are intended in Equinox 1950.

	a= pos.phi() - phi[J2000][trType];
	b= pos.theta();
	sb=Math.sin(b);
	cb=Math.cos(b);
	cbsa=cb*Math.sin(a);
	b=-stheta[J2000][trType] * cbsa + ctheta[J2000][trType]*sb;
	b=Math.max(-1.0,Math.min(b,1.0));
	bo=Math.asin(b);
	
	a=Math.atan2(ctheta[J2000][trType] * cbsa+ stheta[J2000][trType]*sb,cb*Math.cos(a));
	ao=(a+psi[J2000][trType]+fourPi)%twoPi;

	return new AngularPosition(bo, ao);
				      
    }

    /**Transforms an angular position in degrees in a given coordinate system to a position
       in an other coordinate systems, also in degrees. RA-Dec position are intended in 
       Equinox J2000
       */
    public static final AngularPosition transformInDeg(AngularPosition pos, int trType) 
	throws Exception {

	double ao,bo,a,b,sb,cb,cbsa;
	int J2000 = 1;
	//by setting J2000 = 0, RA-Dec are intended in Equinox 1950.

	a= pos.phi()/degToRad - phi[J2000][trType];
	b= pos.theta()/degToRad;
	sb=Math.sin(b);
	cb=Math.cos(b);
	cbsa=cb*Math.sin(a);
	b=-stheta[J2000][trType] * cbsa + ctheta[J2000][trType]*sb;
	b=Math.max(-1.0,Math.min(b,1.0));
	bo=Math.asin(b)*degToRad;
	
	a=Math.atan2(ctheta[J2000][trType] * cbsa+ stheta[J2000][trType]*sb,cb*Math.cos(a));
	ao= ((a+psi[J2000][trType]+fourPi)%twoPi)*degToRad;

	return new AngularPosition(bo, ao);
				      
    }
}
    
