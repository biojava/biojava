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
package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.TmpAtomCache;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;

import junit.framework.TestCase;

public class TestSmallAlignment extends TestCase{


	public void test1a4w(){

		String name1 = "1a4w.I";
		String name2 = "1a4w.H";
		try {
			Atom[] ca1 = TmpAtomCache.cache.getAtoms(name1);
			Atom[] ca2 = TmpAtomCache.cache.getAtoms(name2);

			StructureAlignment ce =StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

			AFPChain afpChain = ce.align(ca1, ca2);

			assertTrue(afpChain != null);

			assertTrue(afpChain.getNrEQR() ==0 );
			
			String xml = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
			
			AFPChain newChain = AFPChainXMLParser.fromXML(xml, ca1, ca2);
			assertTrue(newChain != null);
			
			String xml2 = AFPChainXMLConverter.toXML(newChain,ca1,ca2);
			
			assertEquals(xml,xml2);
			
			
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

	}
}
