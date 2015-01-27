/**
 * 
 */
package org.biojava.bio.structure.symmetry.jmolScript;

import org.biojava.bio.structure.symmetry.core.RotationAxisAligner;
import org.biojava.bio.structure.symmetry.geometry.Prism;
import org.biojava.bio.structure.symmetry.geometry.RectangularPrism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorCn extends JmolSymmetryScriptGeneratorPointGroup {

	public JmolSymmetryScriptGeneratorCn(RotationAxisAligner axisTransformation, String name) {
		super(axisTransformation, name);
		if (axisTransformation.getRotationGroup().getPointGroup().equals("C2")) {
			setPolyhedron(new RectangularPrism(axisTransformation.getDimension().z*2, axisTransformation.getDimension().x*2, axisTransformation.getDimension().y*2));
		} else {
			Prism p = new Prism(axisTransformation.getRotationGroup().getRotation(0).getFold());
			p.setHeight(axisTransformation.getDimension().z*2);
			p.setInscribedRadius(axisTransformation.getRadius());
			setPolyhedron(p);
		}
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		RotationAxisAligner at = getAxisTransformation();
		double polyhedronExtension = Math.max(getPolyhedron().getCirumscribedRadius(), at.getDimension().z);
		
		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
	public int getOrientationCount() {
		//  the last two views (top, bottom) are not that interesting.
		if (getAxisTransformation().getRotationGroup().getPointGroup().equals("C2")) {
		    return getPolyhedron().getViewCount()-2;
		}
		return getPolyhedron().getViewCount();
	}
	
	/**
	 * Returns the name of a specific orientation
	 * @param index orientation index
	 * @return name of orientation
	 */
	public String getOrientationName(int index) {	
		if (getAxisTransformation().getRotationGroup().getPointGroup().equals("C2")) {
			if (index == 0) {
				return "Front C2 axis";
			} else if (index == 2) {
				return "Back C2 axis";
			}
		} 
		return getPolyhedron().getViewName(index);
	}
}
