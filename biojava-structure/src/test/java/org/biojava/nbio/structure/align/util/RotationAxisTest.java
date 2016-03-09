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
 * Created on Oct 7, 2013
 * Author: blivens
 *
 */

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
 * Created on Oct 7, 2013
 * Author: blivens
 *
 */
package org.biojava.nbio.structure.align.util;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * @author blivens
 *
 */
public class RotationAxisTest {

	@Test
	public void testProjection() throws Exception{
		RotationAxis axis;
		Atom dir,pos,projected;

		dir = new AtomImpl();
		pos = new AtomImpl();

		// 180 around z
		dir.setCoords( new double[] {0,0,1});
		pos.setCoords( new double[] {0,0,0} );
		axis = new RotationAxis(dir, pos, Math.PI);

		pos.setCoords( new double[] {1,2,3});
		projected = axis.getProjectedPoint(pos);
		assertArrayEquals(new double[] {0,0,3},projected.getCoords(),1e-14);

		double dist = axis.getProjectedDistance(pos);
		assertEquals(Math.sqrt(5),dist,1e-14);


		// main diagonal through (1,1,0)
		dir.setCoords( new double[] {2,2,2});
		pos.setCoords( new double[] {1,1,0});
		axis = new RotationAxis(dir, pos, Math.PI);

		pos.setCoords( new double[] {1,1,0});
		projected = axis.getProjectedPoint(pos);
		assertArrayEquals(new double[] {1,1,0},projected.getCoords(),1e-14);

		pos.setCoords( new double[] {0,0,-1});
		projected = axis.getProjectedPoint(pos);
		assertArrayEquals(new double[] {0,0,-1},projected.getCoords(),1e-14);

		pos.setCoords( new double[] {-.5,-.5,0});
		projected = axis.getProjectedPoint(pos);
		assertArrayEquals(new double[] {0,0,-1},projected.getCoords(),1e-14);

		dist = axis.getProjectedDistance(pos);
		assertEquals(Math.sqrt(3/2.),dist,1e-14);

		pos.setCoords( new double[] {0,0,0});
		projected = axis.getProjectedPoint(pos);
		assertArrayEquals(new double[] {1/3.,1/3.,-2/3.},projected.getCoords(),1e-14);

	}

}
