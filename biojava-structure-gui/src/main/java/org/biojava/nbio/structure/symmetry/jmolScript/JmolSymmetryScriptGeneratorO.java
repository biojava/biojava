/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
/**
 *
 */
package org.biojava.nbio.structure.symmetry.jmolScript;

import org.biojava.nbio.structure.symmetry.axis.RotationAxisAligner;
import org.biojava.nbio.structure.symmetry.geometry.Octahedron;


/**
 * @author Peter
 *
 */
public class JmolSymmetryScriptGeneratorO extends JmolSymmetryScriptGeneratorPointGroup {

	public JmolSymmetryScriptGeneratorO(RotationAxisAligner axisTransformation, String name) {
		super(axisTransformation, name);
		Octahedron o = new Octahedron();
		double radius = Math.max(axisTransformation.getDimension().z, axisTransformation.getRadius());
		o.setMidRadius(radius);
		setPolyhedron(o);
	}

	@Override
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
