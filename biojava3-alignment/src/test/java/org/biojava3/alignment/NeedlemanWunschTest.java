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

import static org.junit.Assert.*;

import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

public class NeedlemanWunschTest {

    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private NeedlemanWunsch<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() {
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        gaps = new SimpleGapPenalty((short) 10, (short) 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        alignment = new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>(query, target, gaps, blosum62);
        self = new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>(query, query, gaps, blosum62);
    }

    @Test
    public void testNeedlemanWunsch() {
        NeedlemanWunsch<ProteinSequence, AminoAcidCompound> nw =
                new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>();
        nw.setQuery(query);
        nw.setTarget(target);
        nw.setGapPenalty(gaps);
        nw.setSubstitutionMatrix(blosum62);
        assertEquals(nw.getScore(), alignment.getScore());
    }

    @Test
    public void testGetQuery() {
        assertEquals(alignment.getQuery(), query);
        assertEquals(self.getQuery(), query);
    }

    @Test
    public void testGetTarget() {
        assertEquals(alignment.getTarget(), target);
        assertEquals(self.getTarget(), query);
    }

    @Test
    public void testGetGapPenalty() {
        assertEquals(alignment.getGapPenalty(), gaps);
        assertEquals(self.getGapPenalty(), gaps);
    }

    @Test
    public void testGetSubstitutionMatrix() {
        assertEquals(alignment.getSubstitutionMatrix(), blosum62);
        assertEquals(self.getSubstitutionMatrix(), blosum62);
    }

    @Test
    public void testIsStoringScoreMatrix() {
        assertFalse(alignment.isStoringScoreMatrix());
        assertFalse(self.isStoringScoreMatrix());
    }

    @Test
    public void testGetScoreMatrix() {
        short[][] scores = alignment.getScoreMatrix();
        assertEquals(scores[2][1], -6);
        scores = self.getScoreMatrix();
        assertEquals(scores[3][4], 4);
    }

    @Test
    public void testGetScoreMatrixAsString() {
        assertEquals(alignment.getScoreMatrixAsString(), String.format(
                "        R   D   G%n" +
                "    0 -11 -12 -13%n" +
                "A -11  -1 -12 -12%n" +
                "R -12  -6  -3 -14%n" +
                "N -13 -12  -5  -3%n" +
                "D -14 -14  -6  -6%n"));
        assertEquals(self.getScoreMatrixAsString(), String.format(
                "        A   R   N   D%n" +
                "    0 -11 -12 -13 -14%n" +
                "A -11   4  -7  -8  -9%n" +
                "R -12  -7   9  -2  -3%n" +
                "N -13  -8  -2  15   4%n" +
                "D -14  -9  -3   4  21%n"));
    }

    @Test
    public void testGetScoreMatrixAt() {
        assertEquals(alignment.getScoreMatrixAt(1,1), -1);
        assertEquals(self.getScoreMatrixAt(1,1), 4);
    }

    @Test
    public void testGetComputationTime() {
        assertTrue(alignment.getComputationTime() > 0);
        assertTrue(self.getComputationTime() > 0);
    }

    @Test
    public void testGetProfile() {
        assertEquals(alignment.getProfile().toString(), String.format("ARND%n-RDG%n"));
        assertEquals(self.getProfile().toString(), String.format("ARND%nARND%n"));
    }

    @Test
    public void testGetMaxScore() {
        assertEquals(alignment.getMaxScore(), 21);
        assertEquals(self.getMaxScore(), 21);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(alignment.getMinScore(), -27);
        assertEquals(self.getMinScore(), -28);
    }

    @Test
    public void testGetScore() {
        assertEquals(alignment.getScore(), -6);
        assertEquals(self.getScore(), 21);
    }

    @Test
    public void testGetPair() {
        assertEquals(alignment.getPair().toString(), String.format("ARND%n-RDG%n"));
        assertEquals(self.getPair().toString(), String.format("ARND%nARND%n"));
    }

}
