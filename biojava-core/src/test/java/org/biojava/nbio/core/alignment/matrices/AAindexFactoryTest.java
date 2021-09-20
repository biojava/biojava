package org.biojava.nbio.core.alignment.matrices;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class AAindexFactoryTest {
   

    DefaultAAIndexProvider provider = new DefaultAAIndexProvider();
   
    @Test
    void aaProviderIsSingleton(){ 
        AAIndexProvider provider = AAindexFactory.getAAIndexProvider();
        assertNotNull(provider);
        AAIndexProvider provider2 = AAindexFactory.getAAIndexProvider();
        assertTrue(provider == provider2);
    }

    @Test
    void cannotSetProviderToNull(){
        AAindexFactory.setAAIndexProvider(null);
        assertNotNull(AAindexFactory.getAAIndexProvider());
    }

}