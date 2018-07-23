/*
 *                  BioJava development code
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
 * Created on Jun 8, 2007
 *
 */
package org.biojava.nbio.structure.test;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.StructurePairAligner;
import org.biojava.nbio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.InputStream;

public class StructurePairAlignerTest {

	private Structure structure1;
	private Structure structure2;

	@Before
	public void setUp() throws Exception
	{
		InputStream inStream = this.getClass().getResourceAsStream("/5pti.pdb");
		Assert.assertNotNull(inStream);
		InputStream inStream2 = this.getClass().getResourceAsStream("/1tap.pdb");
		Assert.assertNotNull(inStream2);

		PDBFileParser pdbpars = new PDBFileParser();

		structure1 = pdbpars.parsePDBFile(inStream) ;
		structure2 = pdbpars.parsePDBFile(inStream2);


		Assert.assertNotNull(structure1);
		Assert.assertNotNull(structure2);
		Assert.assertEquals("structure does not contain one chain ", 1, structure1.size());

	}


	@Test
	public void testAlignStructureStructure() {
		StructurePairAligner aligner = new StructurePairAligner();

		boolean allFine = true;
		String msg = "";
		try {
			aligner.align(structure1,structure2);

			AlternativeAlignment[] aligs = aligner.getAlignments();
			Assert.assertEquals("the number of obtained alternative alignments is not correct", 20, aligs.length);
			AlternativeAlignment a = aligs[0];

			Assert.assertNotNull(a);

			Assert.assertEquals("the expected nr of eq. residues is not correct.", 47, a.getEqr());

			// they are v. close, but not identical
			Assert.assertTrue(a.getRmsd() < 4);
			Assert.assertTrue(a.getRmsd() > 3);
			Assert.assertTrue(a.getPercId() > 9);
			Assert.assertTrue(a.getScore() > 140);

		} catch (Exception e){
			msg = e.getMessage();
			allFine = false;
			e.printStackTrace();
		}
		Assert.assertTrue(allFine);
		Assert.assertEquals("an error occured", "", msg);



	}



}
