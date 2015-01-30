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

import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.Alignments.ProfileProfileAlignerType;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.GuideTreeNode;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class GuideTreeTest {

    private List<ProteinSequence> proteins;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private GuideTree<ProteinSequence, AminoAcidCompound> tree;

    @Before
    public void setup() throws CompoundNotFoundException { 
        proteins = Arrays.asList(new ProteinSequence[] {new ProteinSequence("ARND"), new ProteinSequence("ARND"),
                new ProteinSequence("HILK"), new ProteinSequence("ANDR")});
        gaps = new SimpleGapPenalty((short) 2, (short) 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        tree = new GuideTree<ProteinSequence, AminoAcidCompound>(proteins, Alignments.getAllPairsScorers(proteins,
                PairwiseSequenceScorerType.GLOBAL_IDENTITIES, gaps, blosum62));
    }

    @Test
    public void testGuideTree() {
        assertNotNull(tree);
    }

    @Test
    public void testGetAllPairsScores() {
        assertArrayEquals(tree.getAllPairsScores(), new double[] {4, 0, 3, 0, 3, 0}, 0.00001);
    }

    @Test
    public void testGetDistanceMatrix() {
        assertArrayEquals(tree.getDistanceMatrix(), new double[][] {
                {0.0, 0.0, 1.0, 0.4},
                {0.0, 0.0, 1.0, 0.4},
                {1.0, 1.0, 0.0, 1.0},
                {0.4, 0.4, 1.0, 0.0}});
    }

    @Test
    public void testGetRoot() {
        assertEquals(Alignments.getProgressiveAlignment(tree, ProfileProfileAlignerType.GLOBAL, gaps,
                blosum62).toString(), String.format("%s%n%s%n%s%n%s%n",
                "--ARND-",
                "--ARND-",
                "HILK---",
                "--A-NDR"));
    }

    @Test
    public void testGetScoreMatrix() {
        assertArrayEquals(tree.getScoreMatrix(), new double[][] {
                {4, 4, 0, 3},
                {4, 6, 0, 3},
                {0, 0, 5, 0},
                {3, 3, 0, 6}});
    }

    @Test
    public void testGetSequences() {
        List<ProteinSequence> list = tree.getSequences();
        assertEquals(list.size(), 4);
        assertEquals(list.get(0), proteins.get(0));
        assertEquals(list.get(1), proteins.get(1));
        assertEquals(list.get(2), proteins.get(2));
        assertEquals(list.get(3), proteins.get(3));
    }

    @Test
    public void testIterator() {
        int i = 0;
        for (GuideTreeNode<ProteinSequence, AminoAcidCompound> n : tree) {
            switch (i++) {
            case 0: assertEquals(n.getName(), "1"); break;
            case 1: assertEquals(n.getName(), "2"); break;
            case 2: assertEquals(n.getName(), ""); break;
            case 3: assertEquals(n.getName(), "3"); break;
            case 4: assertEquals(n.getName(), ""); break;
            case 5: assertEquals(n.getName(), "4"); break;
            case 6: assertEquals(n.getName(), ""); break;
            }
        }
    }

    @Test
    public void testToString() {
        assertEquals(tree.toString(),
                "(((1:0.0,2:0.0):-1.4,3:0.8999999999999999):-0.7,4:-0.7)");
    }

}
