/*
 * HEALPix Java code supported by the Gaia project.
 * Copyright (C) 2006-2011 Gaia Data Processing and Analysis Consortium
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
package healpix.plot3d.gui.healpix3d;

import healpix.plot3d.canvas3d.Group3DAxis;
import healpix.plot3d.canvas3d.Group3DCircle;

import javax.media.j3d.Alpha;
import javax.media.j3d.AmbientLight;
import javax.media.j3d.Appearance;
import javax.media.j3d.BoundingSphere;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.ColoringAttributes;
import javax.media.j3d.DirectionalLight;
import javax.media.j3d.LineAttributes;
import javax.media.j3d.RotationInterpolator;
import javax.media.j3d.TransformGroup;
import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3f;

import com.sun.j3d.utils.behaviors.mouse.MouseRotate;
import com.sun.j3d.utils.behaviors.mouse.MouseTranslate;
import com.sun.j3d.utils.behaviors.mouse.MouseZoom;
import com.sun.j3d.utils.universe.SimpleUniverse;


public class HealCanvas extends Canvas3D implements RotateAble {

    private static final long serialVersionUID = 1L;
    protected TransformGroup objTrans;
    protected TransformGroup objScale;
    protected BranchGroup	scene;
    protected BranchGroup	axisGroup;
    protected BranchGroup	equatorGroup;
    protected BranchGroup	zoneGroup;
    protected BranchGroup	ringGroup;
    protected BranchGroup	ringsGroup;
    protected BranchGroup	solidFaceGroup;
    protected BranchGroup	nestGroup;
    protected BranchGroup	faceGroup;
    protected SimpleUniverse uni=null;
    protected MouseZoom zoom;
    protected RotationInterpolator rotator;
    protected boolean zoneView=true;
    protected boolean nestView=true;
    protected boolean faceView=false;
    protected boolean ringsView=false;
    protected boolean ringView=true;
    protected boolean axisView=true;
    protected boolean solidFaceView = false;
    protected int nside = 8;
    protected int zoneNo = 3;
    protected int faceNo = 0;
    protected int ringNo = 5;
    

	public RotationInterpolator getRotationInterpolator (){
		return rotator;
	}
	
    public BranchGroup createSceneGraph() {

	// Create the root of the branch graph
	BranchGroup objRoot = new BranchGroup();
	objRoot.setCapability(BranchGroup.ALLOW_DETACH);

	axisGroup = new Group3DAxis();
	equatorGroup = new BranchGroup();        // The sphere equator


	// Create the transform group node and initialize it to the
	// identity.  Enable the TRANSFORM_WRITE capability so that
	// our behavior code can modify it at runtime.  Add it to the
	// root of the subgraph.
	
	objScale = new TransformGroup();
	objTrans = new TransformGroup();
	axisGroup.setCapability(BranchGroup.ALLOW_DETACH);
	equatorGroup.setCapability(BranchGroup.ALLOW_DETACH);

	objTrans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
	objTrans.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
	objTrans.setCapability(TransformGroup.ALLOW_CHILDREN_READ);
	objTrans.setCapability(TransformGroup.ALLOW_CHILDREN_WRITE);
	objTrans.setCapability(TransformGroup.ALLOW_CHILDREN_EXTEND);
	objScale.addChild(objTrans);
	objRoot.addChild(objScale);

	makeGrids(nside);// nside==8 - different grids

	BoundingSphere bounds =
	    new BoundingSphere(new Point3d(0.0,0.0,0.0), 100.0);

        // ----------------
        // Generate Equator Branch
        // ----------------

        Group3DCircle equator = new Group3DCircle(1.1,40);
        Appearance equatorApp = equator.getAppearance();
        ColoringAttributes equatorColor = new ColoringAttributes();
        LineAttributes eqLine = new LineAttributes();
        eqLine.setLineWidth(3.0f);
        equatorColor.setColor(1.0f, 1.0f, 0.0f);  // yellow
        equatorApp.setColoringAttributes( equatorColor );
        equatorApp.setLineAttributes( eqLine );
        equatorGroup.addChild(equator);

		objTrans.addChild(equatorGroup);
		objTrans.addChild(axisGroup);
	
        // Set up the background
//        Color3f bgColor = new Color3f(0.99f, 0.99f, 1.0f);
//        Background bgNode = new Background(bgColor);
//        bgNode.setApplicationBounds(bounds);
//        objRoot.addChild(bgNode);

        // Set up the ambient light
        AmbientLight ambientLightNode = new AmbientLight();
        ambientLightNode.setInfluencingBounds(bounds);
        objRoot.addChild(ambientLightNode);

        // Set up the directional lights
        Color3f light1Color = new Color3f(0.8f, 0.8f, 0.8f);
        Vector3f light1Direction  = new Vector3f(4.0f, -7.0f, -12.0f);
        Color3f light2Color = new Color3f(0.3f, 0.3f, 0.3f);
        Vector3f light2Direction  = new Vector3f(-6.0f, -2.0f, -1.0f);

        DirectionalLight light1
            = new DirectionalLight(light1Color,light1Direction);
        light1.setInfluencingBounds(bounds);
        objRoot.addChild(light1);
        DirectionalLight light2
            = new DirectionalLight(light2Color, light2Direction);
        light2.setInfluencingBounds(bounds);
        objRoot.addChild(light2);

          // Create the rotate behavior node
          //MouseRotateY behavior = new MouseRotateY();
          MouseRotate behavior = new MouseRotate();
          behavior.setTransformGroup(objTrans);
          objTrans.addChild(behavior);
          behavior.setSchedulingBounds(bounds);

          // Create the zoom behavior node
          zoom = new MouseZoom();
          zoom.setTransformGroup(objTrans);
          objTrans.addChild(zoom);
          zoom.setSchedulingBounds(bounds);

          // Create the translate behavior node
          MouseTranslate behavior3 = new MouseTranslate();
          behavior3.setTransformGroup(objTrans);
          objTrans.addChild(behavior3);
          behavior3.setSchedulingBounds(bounds);


		// Auto rotator 
		Alpha mover = new Alpha (0,4000);
		rotator = new RotationInterpolator(mover,objTrans);
        rotator.setSchedulingBounds(bounds);
        objTrans.addChild(rotator);
		//rotator.setAlpha(null);



	return objRoot;
    }

	protected void makeGrids(int nside) {
		if (zoneView) setZone(zoneNo); else zoneGroup=null;
		if (nestView) setFace(faceNo); else nestGroup=null;
		if (ringView) setRing(ringNo); else ringGroup=null;
		if (faceView) setFaces(); else faceGroup=null;
		if (ringsView) setRings(); else ringsGroup=null;
		if (solidFaceView) setSolidFace(); else solidFaceGroup=null;
	}

    public HealCanvas() {
		super(SimpleUniverse.getPreferredConfiguration());
    }

	public void setupScene() {
	    if (scene != null) {
			scene.detach();
			scene=null;
		}	
		scene = createSceneGraph();
		if (uni == null) uni = new SimpleUniverse(this);
	}

	public void showScene() {
        scene.compile();
		uni.addBranchGraph(scene);
        // This will move the ViewPlatform back a bit so the
        // objects in the scene can be viewed.
        uni.getViewingPlatform().setNominalViewingTransform();
		uni.getViewer().getView().setFieldOfView(1.2);
		uni.getViewer().getView().setFrontClipDistance(0.01);
	}

	public void setNside(int nside) {
			this.nside = nside;
			makeGrids(nside);
	}

	public void setViewAxis(boolean b) {
		if (b) {
			if (!axisView) objTrans.addChild(axisGroup);
		} else {
			if(axisView) axisGroup.detach();
		}
		axisView=b;
	}
	public void setViewZone(boolean b) {
		if (b) {
			if (zoneGroup==null) setZone(zoneNo);	
			if (!zoneView) objTrans.addChild(zoneGroup);
		} else {
			if (zoneGroup!=null && zoneView) zoneGroup.detach();
		}
		zoneView=b;
	}

	public void setViewRings(boolean b) {
		if (b) {
			if (ringsGroup==null) setRings();	
			if (!ringsView) objTrans.addChild(ringsGroup);
		} else {
			if (ringsGroup!=null && ringsView) ringsGroup.detach();
		}
		ringsView=b;
	}
    
        public void setViewSolidFace(boolean b) {
		if (b) {
			if (solidFaceGroup==null) setSolidFace();	
			if (!solidFaceView) objTrans.addChild(solidFaceGroup);
		} else {
			if (solidFaceGroup!=null && solidFaceView) solidFaceGroup.detach();
		}
		solidFaceView=b;
	}

	public void setViewFaces(boolean b) {
		if (b) {
			if (faceGroup==null) setFaces();	
			if (!faceView) objTrans.addChild(faceGroup);
		} else {
			if (faceGroup!=null && faceView) faceGroup.detach();
		}
		faceView=b;
	}

	public void setViewNest(boolean b) {
		if (b) {
			if (nestGroup==null) setFace(faceNo);	
			if (!nestView) objTrans.addChild(nestGroup);
		} else {
			if (nestView&& nestGroup!=null) nestGroup.detach();
		}
		nestView=b;
		
	}
	public void setViewRing(boolean b) {
		if (b) {
			if (ringGroup==null) setRing(ringNo);	
			if (!ringView) objTrans.addChild(ringGroup);
		} else {
			if (ringView && ringGroup!=null) ringGroup.detach();
		}
		ringView=b;
	}

	public void setZone(int z) {
		zoneNo = z;
		if (zoneView && zoneGroup!=null) zoneGroup.detach();
		zoneGroup = new BranchGroup();        // Show one zone
		HealSphere zone = new ZoneSphere(nside,z);
    	zoneGroup.addChild(zone);
		zoneGroup.setCapability(BranchGroup.ALLOW_DETACH);
		objTrans.addChild(zoneGroup);
		if (!zoneView) zoneGroup.detach();
	}

	public void setFace(int f) {
		faceNo = f;
		if (nestView && nestGroup!=null) nestGroup.detach();
		nestGroup = new BranchGroup();        // Show one zone
		HealSphere face = new NestSphere(nside,f);
    	nestGroup.addChild(face);
		nestGroup.setCapability(BranchGroup.ALLOW_DETACH);
		objTrans.addChild(nestGroup);
		if (!nestView) nestGroup.detach();
	}

	public void setRing(int r) {
		ringNo = r;
		if (ringView && ringGroup!=null) ringGroup.detach();
		ringGroup = new BranchGroup();        // Show one zone
		HealSphere s = new RingSphere(nside,ringNo);
    	ringGroup.addChild(s);
		ringGroup.setCapability(BranchGroup.ALLOW_DETACH);
		objTrans.addChild(ringGroup);
		if (!ringView) ringGroup.detach();
	}
 
	public void setFaces() {
		if (faceView&& faceGroup!=null) faceGroup.detach();
		faceGroup = new BranchGroup();        // Show one zone
		//int[] faces = {0,1,4,5,6,7,8};
		int[] faces = {0,1,5,8,9};
		float[][] col = { {1.0f,0.0f,0.0f},{0.0f,1.0f,0.0f},
							{0.0f,0.0f,1.0f},{0.8f,0.6f,0.6f},
							{0.6f,0.8f,0.6f},{1.0f,0.0f,0.0f}};
		for (int f = 0; f < faces.length; f++) {
			HealSphere face = new NestSphere(nside,faces[f]);
			(face.getAppearance().getColoringAttributes()).setColor(col[f][0],col[f][1],col[f][2]);
    		faceGroup.addChild(face);
		}
		faceGroup.setCapability(BranchGroup.ALLOW_DETACH);
		objTrans.addChild(faceGroup);
		if (!faceView) faceGroup.detach();
	}

	public void setRings() {
	    if (ringsView&& ringsGroup!=null) ringsGroup.detach();
	    ringsGroup = new BranchGroup();        // Show one zone
	    //double[][] col = { {0.9f,0.1f,0.1f}, {0.1f,0.1f,0.9f}, {0.4f,0.6f,0.4f}};
	    float[][] col = { {0.6f,0.4f,0.4f}, {0.1f,0.1f,0.9f}};
	    for (int f = 0; f < (int)(nside*1.5); f+=3) {
		HealSphere ring = new RingSphere(nside,f);
		int colo = f%col.length;
		(ring.getAppearance().getColoringAttributes()).setColor(col[colo][0],col[colo][1],col[colo][2]);
    		ringsGroup.addChild(ring);
	    }
	    ringsGroup.setCapability(BranchGroup.ALLOW_DETACH);
	    objTrans.addChild(ringsGroup);
	    if (!ringsView) ringsGroup.detach();
	}


    public void setSolidFace() {

	if (solidFaceView && solidFaceGroup !=null) solidFaceGroup.detach();

	solidFaceGroup = new BranchGroup();
	int nfaces =12;
	float[] color = new float[3];
	color[0] = 1.0f; // red
	color[1] = 0.1f; // green
		
	for (int f = 0; f < nfaces; f++) {
	    FaceSphere face = new FaceSphere(nside,f);
	    float b = (float) f / (float) nfaces; //farction of blue

	    (face.getAppearance().getColoringAttributes()).setColor((1-b)*color[0],color[1],b);
	    solidFaceGroup.addChild(face);
	}

	solidFaceGroup.setCapability(BranchGroup.ALLOW_DETACH);
	objTrans.addChild(solidFaceGroup);
	if (!solidFaceView) solidFaceGroup.detach();
    }

     /*public void setSolidFace() {
		if (solidFaceView&& solidFaceGroup !=null) solidFaceGroup.detach();
		solidFaceGroup = new BranchGroup();        // Show one zone
		HealSphere face = new FaceSphere(nside,5);
		solidFaceGroup.addChild(face);
		solidFaceGroup.setCapability(BranchGroup.ALLOW_DETACH);
		objTrans.addChild(solidFaceGroup);
		if (!solidFaceView) solidFaceGroup.detach();
		}*/
    

}
