package org.biojava.nbio.core;

import junit.framework.TestCase;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.RNASequenceCreator;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.DNAToRNATranslator;
import org.junit.Test;

/**
 * A Test case for https://github.com/biojava/biojava/issues/344
 *
 * Created by andreas on 12/4/15.
 */

public class TestAmbiguityCompoundSet extends TestCase{

    @Test
    public void testCompountSet(){
        try {

            CompoundSet<NucleotideCompound> dnaSet = AmbiguityDNACompoundSet.getDNACompoundSet();
            CompoundSet<NucleotideCompound> rnaSet = AmbiguityRNACompoundSet.getRNACompoundSet();

            DNASequence dna=new DNASequence("AGTCS", dnaSet);

            assertEquals("AGTCS",dna.toString());

            RNASequence rna = dna.getRNASequence();

            rna = new RNASequence(dna.getSequenceAsString().replaceAll("T", "U"), AmbiguityRNACompoundSet.getRNACompoundSet()); //fails with missing compound S

            assertEquals("AGUCS",rna.toString());

            /** now, do the translation also using the underlying API (should not be needed for a user)
             *
             */
            DNAToRNATranslator translator = new DNAToRNATranslator(new RNASequenceCreator(rnaSet
                    ),dnaSet,rnaSet,false);

            Sequence<NucleotideCompound> translated = translator.createSequence(dna);

            assertEquals("AGUCS", translated.toString());
            
        } catch (CompoundNotFoundException e) {
            e.printStackTrace();
            fail(e.getMessage());
        }
    }
}
