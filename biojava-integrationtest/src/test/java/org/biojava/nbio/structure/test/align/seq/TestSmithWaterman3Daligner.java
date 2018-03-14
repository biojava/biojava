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
package org.biojava.nbio.structure.test.align.seq;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.seq.SmithWaterman3DParameters;
import org.biojava.nbio.structure.align.seq.SmithWaterman3Daligner;
import org.junit.Test;

/**
 * Test the superposition based on a sequence alignment on different cases.
 * 
 * @author Aleix Lafita
 *
 */
public class TestSmithWaterman3Daligner {

	/**
	 * Test the changes in the alignment by the max RMSD parameter. Use the
	 * chain A and B of hemoglobin, because they have very similar sequences,
	 * but some columns of the alignment have to be removed for a better
	 * superposition (lower RMSD).
	 */
	@Test
	public void testMaxRMSD() throws StructureException, IOException {

		Structure s1 = StructureTools.getStructure("1A3N.A");
		Structure s2 = StructureTools.getStructure("1A3N.B");

		Atom[] ca1 = StructureTools.getRepresentativeAtomArray(s1);
		Atom[] ca2 = StructureTools.getRepresentativeAtomArray(s2);

		SmithWaterman3Daligner aligner = new SmithWaterman3Daligner();
		SmithWaterman3DParameters params = new SmithWaterman3DParameters();

		// Use no restriction on the RMSD
		params.setMaxRmsd(99.0);

		AFPChain afpChain = aligner.align(ca1, ca2, params);

		assertEquals("RMSD is wrong", 1.39, afpChain.getTotalRmsdOpt(), 0.005);
		assertEquals("Length is wrong", 137, afpChain.getOptLength());
		
		// Restrict it to 1A RMSD (18 columns have to be dropped)
		params.setMaxRmsd(1.0);

		afpChain = aligner.align(ca1, ca2, params);

		assertTrue("RMSD is above the threshold", afpChain.getTotalRmsdOpt() < 1.0);
		assertEquals("Length is wrong", 119, afpChain.getOptLength());
		
		// Restrict it to 0A RMSD (the minimum length is relevant )
		params.setMaxRmsd(0.0);
		params.setMinLen(30);

		afpChain = aligner.align(ca1, ca2, params);

		assertEquals("Length shoild be the minimum possible", 30, afpChain.getOptLength());

	}
}
