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
 * created at 20 Mar 2014
 * Author: ap3 
 */

package org.biojava.nbio.structure.test.io;


import junit.framework.TestCase;
import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.StructureIO;

public class TestSeqResParsing extends TestCase {

	public void test11GS(){

		String pdbID = "11GS";

		Structure s;
		try {
			AtomCache cache = new AtomCache();
			cache.getFileParsingParams().setAlignSeqRes(true);
			
			StructureIO.setAtomCache(cache);
			
			s = StructureIO.getStructure(pdbID);
			assertNotNull(s);
			assertTrue(s.getChains().size() > 0);
			Chain c = s.getChain(0);

			assertTrue(c.getSeqResGroups().size() > 2);
			
			Group first  = c.getSeqResGroup(0);
			Group second = c.getSeqResGroup(1);
			Group third  = c.getSeqResGroup(2);

			assertTrue(first instanceof AminoAcid);
			assertTrue(second instanceof AminoAcid);
			assertTrue(third instanceof AminoAcid);
			
			AminoAcid aafirst = (AminoAcid) first;
			AminoAcid aasecond = (AminoAcid)second;
			AminoAcid aathird = (AminoAcid) third;
			
			assertTrue(aafirst.getRecordType().equals(AminoAcid.SEQRESRECORD));
			assertTrue(aasecond.getRecordType().equals(AminoAcid.SEQRESRECORD));
			assertTrue(aathird.getRecordType().equals(AminoAcid.ATOMRECORD));
		

		} catch (Exception e) {

			e.printStackTrace();
			fail(e.getMessage());
		}
		

	}
}
