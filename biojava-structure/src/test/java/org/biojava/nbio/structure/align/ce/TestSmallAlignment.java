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
 * Created on May 26, 2010
 *
 */
package org.biojava.nbio.structure.align.ce;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.xml.AFPChainFlipper;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;

public class TestSmallAlignment {


	@Test
	public void test1a4w() throws IOException, StructureException {

		String name1 = "1a4w.I";
		String name2 = "1a4w.H";
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		StructureAlignment ce =StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

		AFPChain afpChain = ce.align(ca1, ca2);

		Assert.assertTrue(afpChain != null);

		afpChain.setName1(name1);
		afpChain.setName2(name2);

		Assert.assertTrue(afpChain.getNrEQR() == 0);

		String xml = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);

		AFPChain newChain = AFPChainXMLParser.fromXML(xml, ca1, ca2);

		Assert.assertTrue(newChain != null);

		String xml2 = AFPChainXMLConverter.toXML(newChain,ca1,ca2);

		Assert.assertEquals(xml, xml2);

		// test flipping

		AFPChain flipped = AFPChainFlipper.flipChain(afpChain);

		String xmlFlipped = AFPChainXMLConverter.toXML(flipped, ca2, ca1);

		AFPChain flippedRecreated = AFPChainXMLParser.fromXML(xmlFlipped, ca2, ca1);

		AFPChain reverted = AFPChainFlipper.flipChain(flippedRecreated);

		String revertedXML = AFPChainXMLConverter.toXML(reverted, ca1, ca2);
		Assert.assertEquals(xml, revertedXML);

	}
}
