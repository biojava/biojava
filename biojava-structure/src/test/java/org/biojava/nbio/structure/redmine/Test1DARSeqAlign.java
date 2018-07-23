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
package org.biojava.nbio.structure.redmine;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.IOException;

/** test for https://redmine.open-bio.org/issues/3282
 *
 * @author Andreas Prlic
 *
 */
public class Test1DARSeqAlign {

	@Test
	public void test1DAR() throws StructureException, IOException {
		AtomCache cache = new AtomCache();
		FileParsingParameters orig = cache.getFileParsingParams();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);

		cache.setFileParsingParams(params);

		boolean usingReducedChemCompProvider = false;

		ChemCompProvider ccp =ChemCompGroupFactory.getChemCompProvider();
		if (ccp.getClass().getName().contains("ReducedChemCompProvider") ) {
			usingReducedChemCompProvider = true;

			ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
		}

		Structure struc = cache.getStructure("1DAR");

		Chain a = struc.getPolyChainByPDB("A");
		Chain b = struc.getNonPolyChainsByPDB("A").get(0);

		//System.out.println(c.getSeqResGroups());

		Group g = b.getGroupByPDB(ResidueNumber.fromString("692"));
		//System.out.println(g);
		//System.out.println(FileConvert.toPDB(g.getAtom(0)));

		Group g3 = a.getGroupByPDB(ResidueNumber.fromString("689"));
		//System.out.println(g3);
		//System.out.println(FileConvert.toPDB(g3.getAtom(0)));

		assertTrue(! a.getSeqResGroups().contains(g));

		assertTrue( g instanceof NucleotideImpl);

		assertTrue(g.getType().equals(GroupType.NUCLEOTIDE));

		assertTrue( g3.getPDBName().equals("LYS"));
		assertTrue( a.getSeqResGroups().contains(g3));

		assertTrue( g3 instanceof AminoAcid);




		if (usingReducedChemCompProvider)
			ChemCompGroupFactory.setChemCompProvider(ccp);

		cache.setFileParsingParams(orig);
	}
}
