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

import java.util.Arrays;
import java.util.List;

import org.biojava3.alignment.Alignments.PairwiseScorer;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

import org.junit.Before;
import org.junit.Test;

public class GuideTreeTest {

    List<ProteinSequence> proteins;
    GapPenalty gaps;
    SubstitutionMatrix<AminoAcidCompound> blosum62;
    private GuideTree<ProteinSequence, AminoAcidCompound> tree;

    @Before
    public void setup() {
        proteins = Arrays.asList(new ProteinSequence[] {new ProteinSequence("ARND"), new ProteinSequence("ARND"),
                new ProteinSequence("HILK"), new ProteinSequence("ANDR")});
        gaps = new SimpleGapPenalty((short) 2, (short) 1);
        blosum62 = new SimpleSubstitutionMatrix<AminoAcidCompound>();
        tree = new GuideTree<ProteinSequence, AminoAcidCompound>(proteins, PairwiseScorer.GLOBAL_IDENTITIES, gaps,
                blosum62);
    }

    @Test
    public void testGuideTree() {
        assertNotNull(tree);
    }

    @Test
    public void testGetAllPairsScores() {
        assertArrayEquals(tree.getAllPairsScores(), new int[] {4, 0, 3, 0, 3, 0});
    }

    @Test
    public void testGetDistanceMatrix() {
        assertArrayEquals(tree.getDistanceMatrix(), new double[][] {
                {0.0,                 0.0, 1.0, 0.19999999999999996},
                {0.0,                 0.0, 1.0, 0.4},
                {1.0,                 1.0, 0.0, 1.0},
                {0.19999999999999996, 0.4, 1.0, 0.0}});
    }

    @Test
    public void testGetRoot() {
        assertEquals(tree.getRoot().getProfile().toString(), String.format("%s%n%s%n%s%n%s%n",
                "--ARND-",
                "--ARND-",
                "HILK---",
                "--A-NDR"));
    }

    @Test
    public void testGetScoreMatrix() {
        assertArrayEquals(tree.getScoreMatrix(), new int[][] {
                {4, 4, 0, 3},
                {4, 6, 0, 3},
                {0, 0, 5, 0},
                {3, 3, 0, 6}});
    }

    @Test
    public void testToString() {
        assertEquals(tree.toString(),
                "(((1:0.0,2:0.0):0.19999999999999996,3:0.8):0.09999999999999998,4:0.09999999999999998)");
    }

}
