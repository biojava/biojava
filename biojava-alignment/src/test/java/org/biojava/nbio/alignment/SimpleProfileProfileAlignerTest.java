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
 * Created on July 7, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.ProfilePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class SimpleProfileProfileAlignerTest {

	private static final double PRECISION = 0.00000001;
	
    private ProteinSequence protein1, protein2, protein3, protein4;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private Profile<ProteinSequence, AminoAcidCompound> prof1, prof2, prof3, prof4;
    private SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound> sppa1, sppa2, sppa3;
    private ProfilePair<ProteinSequence, AminoAcidCompound> pp1, pp2, all;

    @Before
    public void setup() throws CompoundNotFoundException { 
        protein1 = new ProteinSequence("ARND");
        protein2 = new ProteinSequence("ARND");
        protein3 = new ProteinSequence("HILK");
        protein4 = new ProteinSequence("ANDR");
        gaps = new SimpleGapPenalty(2, 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        prof1 = new SimpleProfile<ProteinSequence, AminoAcidCompound>(protein1);
        prof2 = new SimpleProfile<ProteinSequence, AminoAcidCompound>(protein2);
        prof3 = new SimpleProfile<ProteinSequence, AminoAcidCompound>(protein3);
        prof4 = new SimpleProfile<ProteinSequence, AminoAcidCompound>(protein4);
        sppa1 = new SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound>(prof1, prof2, gaps, blosum62);
        pp1 = sppa1.getPair();
        sppa2 = new SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound>(prof3, prof4, gaps, blosum62);
        pp2 = sppa2.getPair();
        sppa3 = new SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound>(pp1, pp2, gaps, blosum62);
        all = sppa3.getPair();
    }

    @Test
    public void testSimpleProfileProfileAligner() {
        SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound> alig =
                new SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound>();
        alig.setQuery(prof1);
        alig.setTarget(prof2);
        alig.setGapPenalty(gaps);
        alig.setSubstitutionMatrix(blosum62);
        assertEquals(alig.getScore(), sppa1.getScore(), PRECISION);
    }

    @Test
    public void testSimpleProfileProfileAlignerProfileOfSCProfileOfSCGapPenaltySubstitutionMatrixOfC() {
        assertNotNull(sppa1);
        assertNotNull(sppa2);
        assertNotNull(sppa3);
    }

    @Test
    public void testGetQuery() {
        assertEquals(prof1, sppa1.getQuery());
        assertEquals(prof3, sppa2.getQuery());
        assertEquals(pp1, sppa3.getQuery());
    }

    @Test
    public void testGetTarget() {
        assertEquals(prof2, sppa1.getTarget());
        assertEquals(prof4, sppa2.getTarget());
        assertEquals(pp2, sppa3.getTarget());
    }

    @Test
    public void testGetGapPenalty() {
        assertEquals(gaps, sppa1.getGapPenalty());
        assertEquals(gaps, sppa2.getGapPenalty());
        assertEquals(gaps, sppa3.getGapPenalty());
    }

    @Test
    public void testGetSubstitutionMatrix() {
        assertEquals(blosum62, sppa1.getSubstitutionMatrix());
        assertEquals(blosum62, sppa2.getSubstitutionMatrix());
        assertEquals(blosum62, sppa3.getSubstitutionMatrix());
    }

    @Test
    public void testIsStoringScoreMatrix() {
        assertFalse(sppa1.isStoringScoreMatrix());
        assertFalse(sppa2.isStoringScoreMatrix());
        assertFalse(sppa3.isStoringScoreMatrix());
    }

    @Test
    public void testGetScoreMatrix() {
        int[][][] scores = sppa1.getScoreMatrix();
        assertEquals(1, scores[2][1][1]);
        scores = sppa2.getScoreMatrix();
        assertEquals(-7, scores[3][4][0]);
        scores = sppa3.getScoreMatrix();
        assertEquals(1, scores[1][2][2]);
    }

    @Test
    public void testGetScoreMatrixAsString() {
        assertEquals(String.format(
                "Substitution%n" +
                "        A   R   N   D%n" +
                "    0  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "A  -\u221E   4  -4  -6  -7%n" +
                "R  -\u221E  -4   9   1  -2%n" +
                "N  -\u221E  -6   1  15   7%n" +
                "D  -\u221E  -7  -2   7  21%n" +
                "%nDeletion%n" +
                "        A   R   N   D%n" +
                "   -2  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "A  -3  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "R  -4   1  -7  -9 -10%n" +
                "N  -5   0   6  -2  -5%n" +
                "D  -6  -1   5  12   4%n" +
                "%nInsertion%n" +
                "        A   R   N   D%n" +
                "   -2  -3  -4  -5  -6%n" +
                "A  -\u221E  -\u221E   1   0  -1%n" +
                "R  -\u221E  -\u221E  -7   6   5%n" +
                "N  -\u221E  -\u221E  -9  -2  12%n" +
                "D  -\u221E  -\u221E -10  -5   4%n"),
                sppa1.getScoreMatrixAsString());
        assertEquals(String.format(
                "Substitution%n" +
                "        A   N   D   R%n" +
                "    0  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "H  -\u221E  -2  -2  -5  -5%n" +
                "I  -\u221E  -4  -5  -5  -8%n" +
                "L  -\u221E  -5  -7  -9  -7%n" +
                "K  -\u221E  -6  -5  -7  -6%n" +
                "%nDeletion%n" +
                "        A   N   D   R%n" +
                "   -2  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "H  -3  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "I  -4  -5  -5  -8  -8%n" +
                "L  -5  -6  -6  -8  -9%n" +
                "K  -6  -7  -7  -9 -10%n" +
                "%nInsertion%n" +
                "        A   N   D   R%n" +
                "   -2  -3  -4  -5  -6%n" +
                "H  -\u221E  -\u221E  -5  -5  -6%n" +
                "I  -\u221E  -\u221E  -7  -8  -8%n" +
                "L  -\u221E  -\u221E  -8  -9 -10%n" +
                "K  -\u221E  -\u221E  -9  -8  -9%n"),
                sppa2.getScoreMatrixAsString());
        assertEquals(String.format(
                "Substitution%n" +
                "        -   H   I   L   K%n" +
                "    0  -\u221E  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "A  -\u221E   4  -5  -5  -6  -7%n" +
                "R  -\u221E  -4   4  -1  -2   3%n" +
                "N  -\u221E  -6   5   3  -2   0%n" +
                "D  -\u221E  -7   0   7  -1   0%n" +
                "%nDeletion%n" +
                "        -   H   I   L   K%n" +
                "   -2  -\u221E  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "A  -3  -\u221E  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "R  -4   1  -8  -8  -9 -10%n" +
                "N  -5   0   1  -4  -5   0%n" +
                "D  -6  -1   2   0  -5  -1%n" +
                "%nInsertion%n" +
                "        -   H   I   L   K%n" +
                "   -2  -3  -4  -5  -6  -7%n" +
                "A  -\u221E  -\u221E   1   0  -1  -2%n" +
                "R  -\u221E  -\u221E  -7   1   0  -1%n" +
                "N  -\u221E  -\u221E  -9   2   1   0%n" +
                "D  -\u221E  -\u221E -10  -3   4   3%n"),
                sppa3.getScoreMatrixAsString());
    }

    @Test
    public void testGetComputationTime() {
        assertTrue(sppa1.getComputationTime() > 0);
        assertTrue(sppa2.getComputationTime() > 0);
        assertTrue(sppa3.getComputationTime() > 0);
    }

    @Test
    public void testGetProfile() {
        assertEquals(sppa1.getProfile().toString(), String.format("ARND%nARND%n"));
        assertEquals(sppa2.getProfile().toString(), String.format("-HILK%nAND-R%n"));
        assertEquals(sppa3.getProfile().toString(), String.format("ARND--%nARND--%n--HILK%nA-ND-R%n"));
    }

    @Test
    public void testGetMaxScore() {
        assertEquals(21, sppa1.getMaxScore(), PRECISION);
        assertEquals(21, sppa2.getMaxScore(), PRECISION);
        assertEquals(21, sppa3.getMaxScore(), PRECISION);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(-12, sppa1.getMinScore(), PRECISION);
        assertEquals(-12, sppa2.getMinScore(), PRECISION);
        assertEquals(-13, sppa3.getMinScore(), PRECISION);
    }

    @Test
    public void testGetScore() {
        assertEquals(21, sppa1.getScore(), PRECISION);
        assertEquals(-6, sppa2.getScore(), PRECISION);
        assertEquals(3, sppa3.getScore(), PRECISION);
    }

    @Test
    public void testGetPair() {
        assertEquals(String.format("ARND%nARND%n"), pp1.toString());
        assertEquals(String.format("-HILK%nAND-R%n"), pp2.toString());
        assertEquals(String.format("ARND--%nARND--%n--HILK%nA-ND-R%n"), all.toString());
    }

}
