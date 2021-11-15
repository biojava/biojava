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
		assertNotNull( proteinSequence.toString());
        assertEquals(24, proteinSequence.getLength());

		StringProxySequenceReader<AminoAcidCompound> sequenceStringProxyLoader = new StringProxySequenceReader<AminoAcidCompound>("XRNDCEQGHILKMFPSTWYVBZJA", AminoAcidCompoundSet.getAminoAcidCompoundSet());
		ProteinSequence proteinSequenceFromProxy = new ProteinSequence(sequenceStringProxyLoader);
		assertNotNull( proteinSequenceFromProxy.toString());
        assertEquals(24, proteinSequence.getLength());
	}  
}
