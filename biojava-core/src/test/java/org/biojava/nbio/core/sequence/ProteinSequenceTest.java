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
package org.biojava.nbio.core.sequence;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.loader.StringProxySequenceReader;
import org.junit.jupiter.api.Test;

public class ProteinSequenceTest {

    @Test
    void basicTest() throws Exception {
        ProteinSequence proteinSequence = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");
        assertNotNull(proteinSequence.toString());
        assertEquals(24, proteinSequence.getLength());

        StringProxySequenceReader<AminoAcidCompound> sequenceStringProxyLoader = new StringProxySequenceReader<AminoAcidCompound>(
                "XRNDCEQGHILKMFPSTWYVBZJA", AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence proteinSequenceFromProxy = new ProteinSequence(sequenceStringProxyLoader);
        assertNotNull(proteinSequenceFromProxy.toString());
        assertEquals(24, proteinSequence.getLength());
    }
}
