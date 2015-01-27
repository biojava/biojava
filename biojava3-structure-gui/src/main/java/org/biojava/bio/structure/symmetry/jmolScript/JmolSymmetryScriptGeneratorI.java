/**
 * 
 */
package org.biojava.bio.structure.symmetry.jmolScript;

import org.biojava.bio.structure.symmetry.core.RotationAxisAligner;
import org.biojava.bio.structure.symmetry.geometry.Icosahedron;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorI extends JmolSymmetryScriptGeneratorPointGroup {

	public JmolSymmetryScriptGeneratorI(RotationAxisAligner axisTransformation, String name) {
		super(axisTransformation, name);
		double radius = Math.max(axisTransformation.getDimension().z, axisTransformation.getRadius());
		Icosahedron i = new Icosahedron();
		i.setMidRadius(radius);
		setPolyhedron(i);
	}
	
	public int getZoom() {
		// find maximum extension of structure
		double maxExtension = getMaxExtension();
		// find maximum extension of polyhedron
		double polyhedronExtension = getPolyhedron().getCirumscribedRadius();
		
		int zoom = Math.round((float)(maxExtension/polyhedronExtension * 110));
		if (zoom > 100) {
			zoom = 100;
		}
		return zoom;
	}
	
}
