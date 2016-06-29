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
import org.biojava.nbio.structure.symmetry.geometry.Prism;


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

	@Override
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

	@Override
	public int getOrientationCount() {
		// for Dn point groups the last view is redundant due to symmetry.
		return getPolyhedron().getViewCount()-1;
	}

	/**
	 * Returns the name of a specific orientation
	 * @param index orientation index
	 * @return name of orientation
	 */
	@Override
	public String getOrientationName(int index) {
		if (index == 0 && getAxisTransformation().getRotationGroup().getPointGroup().equals("D2")) {
			return "Front C2 axis";
		} else {
			return getPolyhedron().getViewName(index);
		}
	}
}
