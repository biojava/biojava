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

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class SmithWatermanTest {

	private static final double PRECISION = 0.00000001;
	
    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private SmithWaterman<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() throws CompoundNotFoundException { 
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
        int[][][] scores = alignment.getScoreMatrix();
        assertEquals(scores[2][2][2], 2);
        assertEquals(scores[4][3][0], 11);
        scores = self.getScoreMatrix();
        assertEquals(scores[2][2][0], 9);
        assertEquals(scores[4][3][1], 11);
    }

    @Test
    public void testGetScoreMatrixAsString() {
        assertEquals(alignment.getScoreMatrixAsString(), String.format(
                "Substitution%n" +
                "      E  R  D  N  K  G  F  P  S%n" +
                "   0  0  0  0  0  0  0  0  0  0%n" +
                "A  0  0  0  0  0  0  0  0  0  1%n" +
                "E  0  5  0  2  0  1  0  0  0  0%n" +
                "R  0  0 10  0  2  2  0  0  0  0%n" +
                "N  0  0  2 11 13  6  5  1  1  3%n" +
                "D  0  2  0 13 12 12  9  6  7  7%n" +
                "K  0  1  4  5 13 17 10  6  7  7%n" +
                "K  0  1  3  4 10 18 15 11 12 12%n" +
                "%nDeletion%n" +
                "      E  R  D  N  K  G  F  P  S%n" +
                "   0  0  0  0  0  0  0  0  0  0%n" +
                "A  0  0  0  0  0  0  0  0  0  0%n" +
                "E  0  0  0  0  0  0  0  0  0  0%n" +
                "R  0  2  0  0  0  0  0  0  0  0%n" +
                "N  0  1  7  0  0  0  0  0  0  0%n" +
                "D  0  0  6  8 10  3  2  0  0  0%n" +
                "K  0  0  5 10  9  9  6  3  4  4%n" +
                "K  0  0  4  9 10 14  7  3  4  4%n" +
                "%nInsertion%n" +
                "      E  R  D  N  K  G  F  P  S%n" +
                "   0  0  0  0  0  0  0  0  0  0%n" +
                "A  0  0  0  0  0  0  0  0  0  0%n" +
                "E  0  0  2  1  0  0  0  0  0  0%n" +
                "R  0  0  0  7  6  5  4  3  2  1%n" +
                "N  0  0  0  0  8 10  9  8  7  6%n" +
                "D  0  0  0  0 10  9  9  8  7  6%n" +
                "K  0  0  0  1  2 10 14 13 12 11%n" +
                "K  0  0  0  0  1  7 15 14 13 12%n"));
        assertEquals(self.getScoreMatrixAsString(), String.format(
                "Substitution%n" +
                "      A  E  R  N  D  K  K%n" +
                "   0  0  0  0  0  0  0  0%n" +
                "A  0  4  0  0  0  0  0  0%n" +
                "E  0  0  9  1  0  2  1  1%n" +
                "R  0  0  1 14  6  3  6  5%n" +
                "N  0  0  0  6 20 12 10  9%n" +
                "D  0  0  2  3 12 26 16 15%n" +
                "K  0  0  1  6 10 16 31 28%n" +
                "K  0  0  1  5  9 15 28 36%n" +
                "%nDeletion%n" +
                "      A  E  R  N  D  K  K%n" +
                "   0  0  0  0  0  0  0  0%n" +
                "A  0  0  0  0  0  0  0  0%n" +
                "E  0  1  0  0  0  0  0  0%n" +
                "R  0  0  6  0  0  0  0  0%n" +
                "N  0  0  5 11  3  0  3  2%n" +
                "D  0  0  4 10 17  9  7  6%n" +
                "K  0  0  3  9 16 23 13 12%n" +
                "K  0  0  2  8 15 22 28 25%n" +
                "%nInsertion%n" +
                "      A  E  R  N  D  K  K%n" +
                "   0  0  0  0  0  0  0  0%n" +
                "A  0  0  1  0  0  0  0  0%n" +
                "E  0  0  0  6  5  4  3  2%n" +
                "R  0  0  0  0 11 10  9  8%n" +
                "N  0  0  0  0  3 17 16 15%n" +
                "D  0  0  0  0  0  9 23 22%n" +
                "K  0  0  0  0  3  7 13 28%n" +
                "K  0  0  0  0  2  6 12 25%n"));
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
        assertEquals(alignment.getMaxScore(), 50, PRECISION);
        assertEquals(self.getMaxScore(), 36, PRECISION);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(alignment.getMinScore(), 0, PRECISION);
        assertEquals(self.getMinScore(), 0, PRECISION);
    }

    @Test
    public void testGetScore() {
        assertEquals(alignment.getScore(), 18, PRECISION);
        assertEquals(self.getScore(), 36, PRECISION);
    }

    @Test
    public void testGetPair() {
        assertEquals(alignment.getPair().toString(), String.format("ERNDKK%nER-DNK%n"));
        assertEquals(self.getPair().toString(), String.format("AERNDKK%nAERNDKK%n"));
    }

}
