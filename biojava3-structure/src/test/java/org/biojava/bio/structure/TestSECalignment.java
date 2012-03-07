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
 * Created on Jan 25, 2010
 *
 */
package org.biojava.bio.structure;

import java.io.InputStream;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.util.StringManipulationTestsHelper;
import org.biojava3.core.util.StringManipulationHelper;

/** This test makes sure that the new representation of selenocysteins as SEC amino acids does not
 * affect the structure alignment results.
 * 
 * @author andreas
 *
 */
public class TestSECalignment extends  TestCase {

	public void testOldSecOutput() throws Exception {

		InputStream inStream = this.getClass().getResourceAsStream("/ce_1fdo.A_2iv2.X.out");
		assertNotNull(inStream);
		String xml = StringManipulationHelper.convertStreamToString(inStream);

		AtomCache cache = new AtomCache();
		String name1="1FDO.A";
		String name2="2IV2.X";
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		AFPChain afpChainOrig = AFPChainXMLParser.fromXML(xml, ca1, ca2);

		// calc time is hardware dependent.... overwrite...
		afpChainOrig.setCalculationTime(-1);

		//String ce1 = afpChainOrig.toFatcat(ca1, ca2);

		String xmlComp =  AFPChainXMLConverter.toXML(afpChainOrig, ca1, ca2);

		//			assertEquals( xml, xmlComp);
		StringManipulationTestsHelper.assertEqualsIgnoreEndline(xml, xmlComp);
		StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

		AFPChain afpChainNew = ce.align(ca1,ca2);
		afpChainNew.setCalculationTime(-1);
		afpChainNew.setName1(name1);
		afpChainNew.setName2(name2);

		String xmlNew = AFPChainXMLConverter.toXML(afpChainNew,ca1,ca2);
		//String ce2 = afpChainNew.toFatcat(ca1, ca2);
		// FIXME version number, new line character difference?
		//			assertEquals(xml,xmlNew);
		StringManipulationTestsHelper.assertEqualsIgnoreEndline(xml,xmlNew);

		//assertEquals(ce1,ce2);

	}
}
