package org.biojava.nbio.core.alignment.matrices;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class DefaultAAIndexProviderTest {
   
    private static final String BENS940102 = "BENS940102";

    DefaultAAIndexProvider provider = new DefaultAAIndexProvider();
    @Test
    void newAAIndexProviderReturnsNullIfNotExists(){ 
        assertNull(provider.getMatrix("unknown"));
    }

    @Test
    void aaIndexProviderGetByName(){ 
        SubstitutionMatrix<AminoAcidCompound> matrix =  provider.getMatrix(BENS940102);
        assertNotNull(matrix);
        assertEquals(BENS940102, matrix.getName());
    }
}