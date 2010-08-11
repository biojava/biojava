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

package org.biojava3.alignment;

import static org.junit.Assert.*;

import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.ProfilePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

import org.junit.Before;
import org.junit.Test;

public class SimpleProfileProfileAlignerTest {

    private ProteinSequence protein1, protein2, protein3, protein4;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private Profile<ProteinSequence, AminoAcidCompound> prof1, prof2, prof3, prof4;
    private SimpleProfileProfileAligner<ProteinSequence, AminoAcidCompound> sppa1, sppa2, sppa3;
    private ProfilePair<ProteinSequence, AminoAcidCompound> pp1, pp2, all;

    @Before
    public void setup() {
        protein1 = new ProteinSequence("ARND");
        protein2 = new ProteinSequence("ARND");
        protein3 = new ProteinSequence("HILK");
        protein4 = new ProteinSequence("ANDR");
        gaps = new SimpleGapPenalty((short) 2, (short) 1);
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
        assertEquals(alig.getScore(), sppa1.getScore());
    }

    @Test
    public void testSimpleProfileProfileAlignerProfileOfSCProfileOfSCGapPenaltySubstitutionMatrixOfC() {
        assertNotNull(sppa1);
        assertNotNull(sppa2);
        assertNotNull(sppa3);
    }

    @Test
    public void testGetQuery() {
        assertEquals(sppa1.getQuery(), prof1);
        assertEquals(sppa2.getQuery(), prof3);
        assertEquals(sppa3.getQuery(), pp1);
    }

    @Test
    public void testGetTarget() {
        assertEquals(sppa1.getTarget(), prof2);
        assertEquals(sppa2.getTarget(), prof4);
        assertEquals(sppa3.getTarget(), pp2);
    }

    @Test
    public void testGetGapPenalty() {
        assertEquals(sppa1.getGapPenalty(), gaps);
        assertEquals(sppa2.getGapPenalty(), gaps);
        assertEquals(sppa3.getGapPenalty(), gaps);
    }

    @Test
    public void testGetSubstitutionMatrix() {
        assertEquals(sppa1.getSubstitutionMatrix(), blosum62);
        assertEquals(sppa2.getSubstitutionMatrix(), blosum62);
        assertEquals(sppa3.getSubstitutionMatrix(), blosum62);
    }

    @Test
    public void testIsStoringScoreMatrix() {
        assertFalse(sppa1.isStoringScoreMatrix());
        assertFalse(sppa2.isStoringScoreMatrix());
        assertFalse(sppa3.isStoringScoreMatrix());
    }

    @Test
    public void testGetScoreMatrix() {
        short[][][] scores = sppa1.getScoreMatrix();
        assertEquals(scores[2][1][1], 1);
        scores = sppa2.getScoreMatrix();
        assertEquals(scores[3][4][0], -7);
        scores = sppa3.getScoreMatrix();
        assertEquals(scores[1][2][2], 1);
    }

    @Test
    public void testGetScoreMatrixAsString() {
        assertEquals(sppa1.getScoreMatrixAsString(), String.format(
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
                "D  -\u221E  -\u221E -10  -5   4%n"));
        assertEquals(sppa2.getScoreMatrixAsString(), String.format(
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
                "K  -\u221E  -\u221E  -9  -8  -9%n"));
        assertEquals(sppa3.getScoreMatrixAsString(), String.format(
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
                "D  -\u221E  -\u221E -10  -3   4   3%n"));
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
        assertEquals(sppa1.getMaxScore(), 21);
        assertEquals(sppa2.getMaxScore(), 21);
        assertEquals(sppa3.getMaxScore(), 21);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(sppa1.getMinScore(), -12);
        assertEquals(sppa2.getMinScore(), -12);
        assertEquals(sppa3.getMinScore(), -13);
    }

    @Test
    public void testGetScore() {
        assertEquals(sppa1.getScore(), 21);
        assertEquals(sppa2.getScore(), -6);
        assertEquals(sppa3.getScore(), 3);
    }

    @Test
    public void testGetPair() {
        assertEquals(pp1.toString(), String.format("ARND%nARND%n"));
        assertEquals(pp2.toString(), String.format("-HILK%nAND-R%n"));
        assertEquals(all.toString(), String.format("ARND--%nARND--%n--HILK%nA-ND-R%n"));
    }

}
