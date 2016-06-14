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
package org.biojava.nbio.structure.symmetry.internal;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.junit.Test;

/**
 * Run two easy cases of internal symmetry and test that results are significant
 * and order is correct.
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class TestCeSymm {

	@Test
	public void testEasyCases() throws IOException, StructureException {

		String[] names = new String[] { "1hiv.A", "4i4q", "1n0r.A"};
		int[] orders = new int[] { 2, 3, 4 };

		for (int i = 0; i < names.length; i++) {

			Structure s = StructureTools.getStructure(names[i]);
			Atom[] atoms = StructureTools.getRepresentativeAtomArray(s);

			CeSymmResult result = CeSymm.analyze(atoms);

			assertTrue(result.isSignificant());
			assertEquals(result.getNumRepeats(), orders[i]);
		}
	}
}
