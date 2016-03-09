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
package org.biojava.nbio.structure.symmetry.geometry;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import java.util.List;

/**
 * @author Peter
 *
 */
public interface Polyhedron {

	public Point3d[] getVertices();
	public List<int[]> getLineLoops();
	public double getCirumscribedRadius();
	public int getViewCount();
	public String getViewName(int index);
	public Matrix3d getViewMatrix(int index);
}
