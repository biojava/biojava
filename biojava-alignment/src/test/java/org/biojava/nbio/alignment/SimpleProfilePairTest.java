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

import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.ProfilePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class SimpleProfilePairTest {

    private ProteinSequence protein1, protein2, protein3, protein4;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private Profile<ProteinSequence, AminoAcidCompound> pair1, pair2;

    @Before
    public void setup() throws CompoundNotFoundException { 
        protein1 = new ProteinSequence("ARND");
        protein2 = new ProteinSequence("ARND");
        protein3 = new ProteinSequence("HILK");
        protein4 = new ProteinSequence("ANDR");
        gaps = new SimpleGapPenalty((short) 2, (short) 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        pair1 = new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>(protein1, protein2, gaps, blosum62).getPair();
        pair2 = new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>(protein3, protein4, gaps, blosum62).getPair();
    }

    @Test
    public void testSimpleProfilePair() {
        ProfilePair<ProteinSequence, AminoAcidCompound> all =
                new SimpleProfilePair<ProteinSequence, AminoAcidCompound>(pair1, pair2, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.GAP}), Arrays.asList(
                new Step[] {Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND}));
        assertEquals(all.toString(), String.format("ARND--%nARND--%n--HILK%nA-ND-R%n"));
    }

    @Test
    public void testGetQuery() {
        ProfilePair<ProteinSequence, AminoAcidCompound> all =
                new SimpleProfilePair<ProteinSequence, AminoAcidCompound>(pair1, pair2, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.GAP}), Arrays.asList(
                new Step[] {Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND}));
        assertEquals(all.getQuery(), pair1);
    }

    @Test
    public void testGetTarget() {
        ProfilePair<ProteinSequence, AminoAcidCompound> all =
                new SimpleProfilePair<ProteinSequence, AminoAcidCompound>(pair1, pair2, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.GAP}), Arrays.asList(
                new Step[] {Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND}));
        assertEquals(all.getTarget(), pair2);
    }

}
