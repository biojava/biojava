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
 * Created on June 15, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class SimpleProfileTest {

    private ProteinSequence query, target;
    private Profile<ProteinSequence, AminoAcidCompound> profile;

    @Before
    public void setup() {
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        profile = new SimpleProfile<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}), Arrays.asList(new Step[] {
                Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND}));
    }

    @Test(expected=IllegalArgumentException.class)
    public void testSimpleProfile() {
        new SimpleProfile<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}), Arrays.asList(new Step[] {
                Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND}));
    }

    @Test
    public void testGetAlignedSequenceInt() {
        assertEquals(profile.getAlignedSequence(1).toString(), "ARND-");
        assertEquals(profile.getAlignedSequence(2).toString(), "-R-DG");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds() {
        profile.getAlignedSequence(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds2() {
        profile.getAlignedSequence(3);
    }

    @Test
    public void testGetAlignedSequenceS() {
        assertEquals(profile.getAlignedSequence(query).toString(), "ARND-");
        assertEquals(profile.getAlignedSequence(target).toString(), "-R-DG");
        assertNull(profile.getAlignedSequence(new ProteinSequence("AR")));
    }

    @Test
    public void testGetAlignedSequences() {
        List<AlignedSequence<AminoAcidCompound>> list = profile.getAlignedSequences();
        assertEquals(list.size(), 2);
        assertEquals(list.get(0).toString(), "ARND-");
        assertEquals(list.get(1).toString(), "-R-DG");
    }

    @Test
    public void testGetAlignedSequencesIntArray() {
        List<AlignedSequence<AminoAcidCompound>> list = profile.getAlignedSequences(2, 1, 2);
        assertEquals(list.size(), 3);
        assertEquals(list.get(0).toString(), "-R-DG");
        assertEquals(list.get(1).toString(), "ARND-");
        assertEquals(list.get(2).toString(), "-R-DG");
    }

    @Test
    public void testGetAlignedSequencesSArray() {
        List<AlignedSequence<AminoAcidCompound>> list = profile.getAlignedSequences(query, query, target);
        assertEquals(list.size(), 3);
        assertEquals(list.get(0).toString(), "ARND-");
        assertEquals(list.get(1).toString(), "ARND-");
        assertEquals(list.get(2).toString(), "-R-DG");
    }

    @Test
    public void testGetCompoundAtIntInt() {
        assertEquals(profile.getCompoundAt(1, 4).getShortName(), "D");
        assertEquals(profile.getCompoundAt(2, 3).getShortName(), "-");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds() {
        profile.getCompoundAt(0, 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds2() {
        profile.getCompoundAt(3, 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds3() {
        profile.getCompoundAt(1, 0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds4() {
        profile.getCompoundAt(2, 6);
    }

    @Test
    public void testGetCompoundAtSInt() {
        assertEquals(profile.getCompoundAt(query, 2).getShortName(), "R");
        assertEquals(profile.getCompoundAt(target, 5).getShortName(), "G");
        assertNull(profile.getCompoundAt(new ProteinSequence("AR"), 3));
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds() {
        profile.getCompoundAt(query, 0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds2() {
        profile.getCompoundAt(target, 6);
    }

    @Test
    public void testGetCompoundSet() {
        assertEquals(profile.getCompoundSet(), AminoAcidCompoundSet.getAminoAcidCompoundSet());
    }

    @Test
    public void testGetCompoundsAt() {
        List<AminoAcidCompound> column = profile.getCompoundsAt(5);
        assertEquals(column.size(), 2);
        assertEquals(column.get(0).getShortName(), "-");
        assertEquals(column.get(1).getShortName(), "G");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds() {
        profile.getCompoundsAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds2() {
        profile.getCompoundsAt(6);
    }

    @Test
    public void testGetIndexOf() {
        AminoAcidCompoundSet cs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
        assertEquals(profile.getIndexOf(cs.getCompoundForString("A")), 1);
        assertEquals(profile.getIndexOf(cs.getCompoundForString("R")), 2);
        assertEquals(profile.getIndexOf(cs.getCompoundForString("N")), 3);
        assertEquals(profile.getIndexOf(cs.getCompoundForString("D")), 4);
        assertEquals(profile.getIndexOf(cs.getCompoundForString("G")), 5);
        assertEquals(profile.getIndexOf(cs.getCompoundForString("-")), 1);
        assertEquals(profile.getIndexOf(cs.getCompoundForString("E")), -1);
    }

    @Test
    public void testGetIndicesAt() {
        assertArrayEquals(profile.getIndicesAt(1), new int[] {1, 1});
        assertArrayEquals(profile.getIndicesAt(2), new int[] {2, 1});
        assertArrayEquals(profile.getIndicesAt(3), new int[] {3, 1});
        assertArrayEquals(profile.getIndicesAt(4), new int[] {4, 2});
        assertArrayEquals(profile.getIndicesAt(5), new int[] {4, 3});
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds() {
        profile.getIndicesAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds2() {
        profile.getIndicesAt(6);
    }

    @Test
    public void testGetLastIndexOf() {
        AminoAcidCompoundSet cs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("A")), 1);
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("R")), 2);
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("N")), 3);
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("D")), 4);
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("G")), 5);
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("-")), 5);
        assertEquals(profile.getLastIndexOf(cs.getCompoundForString("E")), -1);
    }

    @Test
    public void testGetLength() {
        assertEquals(profile.getLength(), 5);
    }

    @Test
    public void testGetOriginalSequences() {
        List<ProteinSequence> list = profile.getOriginalSequences();
        assertEquals(list.size(), 2);
        assertEquals(list.get(0), query);
        assertEquals(list.get(1), target);
    }

    @Test
    public void testGetSize() {
        assertEquals(profile.getSize(), 2);
    }

    @Ignore // TODO SimpleProfile.getSubProfile(Location)
    @Test
    public void testGetSubProfile() {
        fail("Not yet implemented");
    }

    @Test
    public void testIsCircular() {
        assertFalse(profile.isCircular());
    }

    @Ignore // TODO SimpleProfile.toString(int)
    @Test
    public void testToStringInt() {
        fail("Not yet implemented");
    }

    @Test
    public void testToString() {
        String newLine = System.getProperty("line.separator");
        assertEquals(profile.toString(), "ARND-" + newLine + "-R-DG" + newLine);
    }

    @Test
    public void testIterator() {
        for (AlignedSequence<AminoAcidCompound> s : profile) {
            assertEquals(s.toString().length(), 5);
        }
    }

}
