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
 * Created on November 21, 2010
 * Author: Mark Chapman
 */

package org.biojava3.core.sequence;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.Before;
import org.junit.Test;

public class MultipleSequenceAlignmentTest {

    private MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msaProteins;
    private MultipleSequenceAlignment<DNASequence,NucleotideCompound> msaDNA;

    @Before
    public void setup() {
        msaProteins = new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();
        for (int i = 0; i < 8; i++) {
            msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        }
        msaDNA = new MultipleSequenceAlignment<DNASequence, NucleotideCompound>();
        for (int i = 0; i < 7; i++) {
            msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        }
    }

    @Test
    public void testGetCompoundsAt() {
        AminoAcidCompound aminoAcid = AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("N");
        List<AminoAcidCompound> colProteins = new ArrayList<AminoAcidCompound>();
        for (int i = 0; i < 8; i++) {
            colProteins.add(aminoAcid);
        }
        assertEquals(msaProteins.getCompoundsAt(3), colProteins);
        NucleotideCompound nucleotide = DNACompoundSet.getDNACompoundSet().getCompoundForString("C");
        List<NucleotideCompound> colDNA = new ArrayList<NucleotideCompound>();
        for (int i = 0; i < 7; i++) {
            colDNA.add(nucleotide);
        }
        assertEquals(msaDNA.getCompoundsAt(3), colDNA);
    }

}
