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
 * Created on June 17, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class NeedlemanWunschQuadraticTest {

    private static final double PRECISION = 0.00000001;

    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private NeedlemanWunschQuadratic<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() throws CompoundNotFoundException {
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        gaps = new SimpleGapPenalty((short) 10, (short) 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        alignment = new NeedlemanWunschQuadratic<ProteinSequence, AminoAcidCompound>(query, target, gaps, blosum62);
        self = new NeedlemanWunschQuadratic<ProteinSequence, AminoAcidCompound>(query, query, gaps, blosum62);
    }

    @Test
    public void testNeedlemanWunsch() {
        NeedlemanWunschQuadratic<ProteinSequence, AminoAcidCompound> nw =
                new NeedlemanWunschQuadratic<ProteinSequence, AminoAcidCompound>();
        nw.setQuery(query);
        nw.setTarget(target);
        nw.setGapPenalty(gaps);
        nw.setSubstitutionMatrix(blosum62);
        assertEquals(nw.getScore(), alignment.getScore(), PRECISION);
    }

    @Test
    public void testAnchoredDNAAlignment() throws CompoundNotFoundException {
        DNASequence query = new DNASequence(  "ACGTACCGGTTTT", DNACompoundSet.getDNACompoundSet());
        DNASequence target = new DNASequence("TACGTCCGGTTACGTACGTT", DNACompoundSet.getDNACompoundSet());
        NeedlemanWunschQuadratic<DNASequence, NucleotideCompound> aligner = new NeedlemanWunschQuadratic<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)5, (short)2), SubstitutionMatrixHelper.getNuc4_4());
        assertEquals(String.format("-ACGTACCGGTT-------TT%nTACGT-CCGGTTACGTACGTT%n"), aligner.getPair().toString());
    }
}
