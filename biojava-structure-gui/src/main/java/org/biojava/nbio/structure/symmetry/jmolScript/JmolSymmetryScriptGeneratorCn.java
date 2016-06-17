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
import org.biojava.nbio.structure.symmetry.geometry.RectangularPrism;


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

	@Override
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

	@Override
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
	@Override
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
