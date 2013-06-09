/**
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
 * Created on 2013-05-28
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.io;

import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava3.core.sequence.ProteinSequence;
import org.junit.Before;
import org.junit.Test;


/**
 * A test for {@link FastaAFPChainConverter}.
 * @author dmyersturnbull
 *
 */
public class FastaAFPChainConverterTest {

	private AtomCache cache;
	
	@Before
	public void setUp() {
		cache = new AtomCache();
	}
	
//	@Test
	public void testIncomplete() throws IOException, StructureException {
		Structure s1 = cache.getStructure("1w0p");
		Structure s2 = cache.getStructure("1qdm");
		ProteinSequence seq1 = new ProteinSequence("GWGG----SEL--YRRNTSLNS--QQDW-------QSNAKIRIVDGAA-----NQIQ");
		ProteinSequence seq2 = new ProteinSequence("WMQNQLAQNKT--QDLILDYVNQLCNRL---PSPMESAV----DCGSLGSMPDIEFT");
		FastaAFPChainConverter conv = new FastaAFPChainConverter(true);
		AFPChain afpChain = conv.fastaToAfpChain(seq1, seq2, s1, s2);
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
	}

}
