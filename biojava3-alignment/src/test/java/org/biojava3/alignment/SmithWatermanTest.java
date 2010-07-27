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
 * Created on June 29, 2010
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

public class SmithWatermanTest {

    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private SmithWaterman<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() {
        query = new ProteinSequence("AERNDKK");
        target = new ProteinSequence("ERDNKGFPS");
        gaps = new SimpleGapPenalty((short) 2, (short) 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        alignment = new SmithWaterman<ProteinSequence, AminoAcidCompound>(query, target, gaps, blosum62);
        self = new SmithWaterman<ProteinSequence, AminoAcidCompound>(query, query, gaps, blosum62);
    }

    @Test
    public void testSmithWaterman() {
        SmithWaterman<ProteinSequence, AminoAcidCompound> alig =
                new SmithWaterman<ProteinSequence, AminoAcidCompound>();
        alig.setQuery(query);
        alig.setTarget(target);
        alig.setGapPenalty(gaps);
        alig.setSubstitutionMatrix(blosum62);
        assertEquals(alig.getPair().toString(), String.format("ERNDKK%nER-DNK%n"));
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
    public void testSetStoringScoreMatrix() {
        assertFalse(alignment.isStoringScoreMatrix());
        alignment.setStoringScoreMatrix(true);
        assertTrue(alignment.isStoringScoreMatrix());
    }

    @Test
    public void testGetScoreMatrix() {
        short[][] scores = alignment.getScoreMatrix();
        assertEquals(scores[2][2], 2);
        assertEquals(scores[4][3], 11);
        scores = self.getScoreMatrix();
        assertEquals(scores[2][2], 9);
        assertEquals(scores[4][3], 11);
    }

    @Test
    public void testGetScoreMatrixAsString() {
        assertEquals(alignment.getScoreMatrixAsString(), String.format(
                "      E  R  D  N  K  G  F  P  S%n" +
                "   0  0  0  0  0  0  0  0  0  0%n" +
                "A  0  0  0  0  0  0  0  0  0  1%n" +
                "E  0  5  2  2  0  1  0  0  0  0%n" +
                "R  0  2 10  7  6  5  4  3  2  1%n" +
                "N  0  1  7 11 13 10  9  8  7  6%n" +
                "D  0  2  6 13 12 12  9  8  7  7%n" +
                "K  0  1  5 10 13 17 14 13 12 11%n" +
                "K  0  1  4  9 10 18 15 14 13 12%n"));
        assertEquals(self.getScoreMatrixAsString(), String.format(
                "      A  E  R  N  D  K  K%n" +
                "   0  0  0  0  0  0  0  0%n" +
                "A  0  4  1  0  0  0  0  0%n" +
                "E  0  1  9  6  5  4  3  2%n" +
                "R  0  0  6 14 11 10  9  8%n" +
                "N  0  0  5 11 20 17 16 15%n" +
                "D  0  0  4 10 17 26 23 22%n" +
                "K  0  0  3  9 16 23 31 28%n" +
                "K  0  0  2  8 15 22 28 36%n"));
    }

    @Test
    public void testGetScoreMatrixAt() {
        assertEquals(alignment.getScoreMatrixAt(5,3), 13);
        assertEquals(self.getScoreMatrixAt(5,4), 17);
    }

    @Test
    public void testGetComputationTime() {
        assertTrue(alignment.getComputationTime() > 0);
        assertTrue(self.getComputationTime() > 0);
    }

    @Test
    public void testGetProfile() {
        assertEquals(alignment.getProfile().toString(), String.format("ERNDKK%nER-DNK%n"));
        assertEquals(self.getProfile().toString(), String.format("AERNDKK%nAERNDKK%n"));
    }

    @Test
    public void testGetMaxScore() {
        assertEquals(alignment.getMaxScore(), 50);
        assertEquals(self.getMaxScore(), 36);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(alignment.getMinScore(), 0);
        assertEquals(self.getMinScore(), 0);
    }

    @Test
    public void testGetScore() {
        assertEquals(alignment.getScore(), 18);
        assertEquals(self.getScore(), 36);
    }

    @Test
    public void testGetPair() {
        assertEquals(alignment.getPair().toString(), String.format("ERNDKK%nER-DNK%n"));
        assertEquals(self.getPair().toString(), String.format("AERNDKK%nAERNDKK%n"));
    }

}
