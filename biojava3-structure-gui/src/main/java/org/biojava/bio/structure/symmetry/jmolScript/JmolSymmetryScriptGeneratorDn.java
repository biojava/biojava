/**
 * 
 */
package org.biojava.bio.structure.symmetry.jmolScript;

import org.biojava.bio.structure.symmetry.core.RotationAxisAligner;
import org.biojava.bio.structure.symmetry.geometry.Prism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorDn extends JmolSymmetryScriptGeneratorPointGroup {

	public JmolSymmetryScriptGeneratorDn(RotationAxisAligner axisTransformation, String name) {
		super(axisTransformation, name);
		int fold = axisTransformation.getRotationGroup().getRotation(0).getFold();
		
		// special case for D2. Since there is no 2-fold prism, draw a 4-fold
		// prism that encases the D2 structure
		if (axisTransformation.getRotationGroup().getPointGroup().equals("D2")) {
			fold = 4;
		}
		
		Prism p = new Prism(fold);
		p.setHeight(axisTransformation.getDimension().z*2);
		p.setInscribedRadius(axisTransformation.getRadius());
		setPolyhedron(p);
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		double polyhedronExtension = Math.max(getPolyhedron().getCirumscribedRadius(), getAxisTransformation().getDimension().z);

		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
	public int getOrientationCount() {
		// for Dn point groups the last view is redundant due to symmetry.
		return getPolyhedron().getViewCount()-1;
	}
	
	/**
	 * Returns the name of a specific orientation
	 * @param index orientation index
	 * @return name of orientation
	 */
	public String getOrientationName(int index) {	
		if (index == 0 && getAxisTransformation().getRotationGroup().getPointGroup().equals("D2")) {
			return "Front C2 axis";
		} else {
			return getPolyhedron().getViewName(index);
		}
	}
}
