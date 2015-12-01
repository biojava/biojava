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

import java.io.IOException;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.junit.Test;

/**
 * Test for the old bug of non-independent CeSymm runs.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 * 
 */
public class TestCeSymm {

	@Test
	public void testIndependence() throws IOException, StructureException {

		// Only instantiate one CeSymm class
		CeSymm ce = new CeSymm();
		String name;
		Atom[] atoms;

		name = "1MER.A";
		atoms = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		ce.analyze(atoms);

		name = "1ijq.A:377-641";
		atoms = StructureTools.getRepresentativeAtomArray(StructureTools
				.getStructure(name));
		// This was causing an assertion error if runs are dependent
		ce.analyze(atoms);
	}

}
