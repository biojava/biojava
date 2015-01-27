/**
 * 
 */
package org.biojava.bio.structure.symmetry.jmolScript;

import org.biojava.bio.structure.symmetry.core.RotationAxisAligner;
import org.biojava.bio.structure.symmetry.geometry.RectangularPrism;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorC1 extends JmolSymmetryScriptGeneratorPointGroup {

	public JmolSymmetryScriptGeneratorC1(RotationAxisAligner axisTransformation, String name) {
		super(axisTransformation, name);
		setPolyhedron(new RectangularPrism(axisTransformation.getDimension().z*2, axisTransformation.getDimension().x*2, axisTransformation.getDimension().y*2));
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		RotationAxisAligner at = getAxisTransformation();
		double polyhedronExtension = Math.max(at.getDimension().x, at.getDimension().y);
		
		polyhedronExtension = Math.max(at.getDimension().z, polyhedronExtension);
		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
	public int getOrientationCount() {
		// the last two views (top, bottom) are not that interesting.
		return getPolyhedron().getViewCount()-2;
	}

}
