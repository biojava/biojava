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
package org.biojava.nbio.core.sequence.loader;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;


/**
 * Testing example for issue #834
 *
 * @author Jacek Grzebyta
 * @author Paolo Pavan
 * @see InfoTask
 */
public class SimpleGenbankProxySequenceReaderTest {


	private final static Logger logger = LoggerFactory.getLogger(SimpleGenbankProxySequenceReaderTest.class);

	@Test(expected = IOException.class)
	public void testWrongSequence() throws Exception {
		logger.info("test wrong sequence");

		String wrongGi = "34567";

        GenbankProxySequenceReader<AminoAcidCompound> genbankReader
				= new GenbankProxySequenceReader<>(System.getProperty("java.io.tmpdir"),
				wrongGi,
                AminoAcidCompoundSet.aminoAcidCompoundSet);

		ProteinSequence seq = new ProteinSequence(genbankReader);
	}
}
