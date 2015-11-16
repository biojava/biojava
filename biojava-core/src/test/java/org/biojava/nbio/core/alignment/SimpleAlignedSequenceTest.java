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

package org.biojava.nbio.core.alignment;

import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.location.SimpleLocation;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.*;

public class SimpleAlignedSequenceTest {

    private ProteinSequence go, lo;
    private AlignedSequence<ProteinSequence, AminoAcidCompound> global, local, local2;
    private AminoAcidCompoundSet cs;

    @Before
    public void setup() throws CompoundNotFoundException { 
        go = new ProteinSequence("ARND");
        lo = new ProteinSequence("CEQGHILKM");
        global = new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>(go, Arrays.asList(new Step[] {
                Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.GAP}));
        local = new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>(lo, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND,
                Step.COMPOUND}), 1, 3);
        local2 = new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>(go, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND}), 1, 0);
        cs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
    }

    @Test(expected=IllegalArgumentException.class)
    public void testSimpleAlignedSequenceLocal() {
        new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>(lo, Arrays.asList(new Step[] {Step.COMPOUND,
                Step.COMPOUND, Step.GAP, Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND}));
    }

    @Test(expected=IllegalArgumentException.class)
    public void testSimpleAlignedSequenceLong() {
        new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>(go, Arrays.asList(new Step[] {Step.GAP,
                Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}));
    }

    @Test(expected=IllegalArgumentException.class)
    public void testSimpleAlignedSequenceShort() {
        new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>(go, Arrays.asList(new Step[] {Step.GAP,
                Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.GAP}));
    }

    @Test
    public void testGetAlignmentIndexAt() {
        assertEquals(global.getAlignmentIndexAt(1), 2);
        assertEquals(global.getAlignmentIndexAt(2), 3);
        assertEquals(global.getAlignmentIndexAt(3), 5);
        assertEquals(global.getAlignmentIndexAt(4), 6);
        assertEquals(local.getAlignmentIndexAt(1), 1);
        assertEquals(local.getAlignmentIndexAt(2), 1);
        assertEquals(local.getAlignmentIndexAt(3), 2);
        assertEquals(local.getAlignmentIndexAt(4), 5);
        assertEquals(local.getAlignmentIndexAt(5), 7);
        assertEquals(local.getAlignmentIndexAt(6), 8);
        assertEquals(local.getAlignmentIndexAt(7), 8);
        assertEquals(local.getAlignmentIndexAt(8), 8);
        assertEquals(local.getAlignmentIndexAt(9), 8);
        assertEquals(local2.getAlignmentIndexAt(1), 1);
        assertEquals(local2.getAlignmentIndexAt(2), 1);
        assertEquals(local2.getAlignmentIndexAt(3), 2);
        assertEquals(local2.getAlignmentIndexAt(4), 3);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignmentIndexAtOutOfBounds() {
        global.getAlignmentIndexAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignmentIndexAtOutOfBounds2() {
        global.getAlignmentIndexAt(5);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignmentIndexAtOutOfBounds3() {
        local.getAlignmentIndexAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignmentIndexAtOutOfBounds4() {
        local.getAlignmentIndexAt(10);
    }

    @Test
    public void testGetEnd() {
        assertEquals(global.getEnd().getPosition(), Integer.valueOf(6));
        assertEquals(local.getEnd().getPosition(), Integer.valueOf(8));
        assertEquals(local2.getEnd().getPosition(), Integer.valueOf(3));
    }

    @Test
    public void testGetLocationInAlignment() {
        assertEquals(global.getLocationInAlignment(), new SimpleLocation(2, 6, Strand.UNDEFINED,
                new SimpleLocation(2, 3, Strand.UNDEFINED), new SimpleLocation(5, 6, Strand.UNDEFINED)));
        assertEquals(local.getLocationInAlignment(), new SimpleLocation(1, 8, Strand.UNDEFINED,
                new SimpleLocation(1, 2, Strand.UNDEFINED), new SimpleLocation(5, 5, Strand.UNDEFINED),
                new SimpleLocation(7, 8, Strand.UNDEFINED)));
        assertEquals(local2.getLocationInAlignment(), new SimpleLocation(1, 3, Strand.UNDEFINED));
    }

    @Test
    public void testGetNumGaps() {
        assertEquals(global.getNumGaps(), 3);
        assertEquals(local.getNumGaps(), 2);
        assertEquals(local2.getNumGaps(), 0);
    }

    @Test
    public void testGetOriginalSequence() {
        assertEquals(global.getOriginalSequence(), go);
        assertEquals(local.getOriginalSequence(), lo);
        assertEquals(local2.getOriginalSequence(), go);
    }

    @Test
    public void testGetOverlapCount() {
        assertEquals(global.getOverlapCount(), 1);
        assertEquals(local.getOverlapCount(), 1);
        assertEquals(local2.getOverlapCount(), 1);
    }

    @Test
    public void testGetSequenceIndexAt() {
        assertEquals(global.getSequenceIndexAt(1), 1);
        assertEquals(global.getSequenceIndexAt(2), 1);
        assertEquals(global.getSequenceIndexAt(3), 2);
        assertEquals(global.getSequenceIndexAt(4), 2);
        assertEquals(global.getSequenceIndexAt(5), 3);
        assertEquals(global.getSequenceIndexAt(6), 4);
        assertEquals(global.getSequenceIndexAt(7), 4);
        assertEquals(local.getSequenceIndexAt(1), 2);
        assertEquals(local.getSequenceIndexAt(2), 3);
        assertEquals(local.getSequenceIndexAt(3), 3);
        assertEquals(local.getSequenceIndexAt(4), 3);
        assertEquals(local.getSequenceIndexAt(5), 4);
        assertEquals(local.getSequenceIndexAt(6), 4);
        assertEquals(local.getSequenceIndexAt(7), 5);
        assertEquals(local.getSequenceIndexAt(8), 6);
        assertEquals(local2.getSequenceIndexAt(1), 2);
        assertEquals(local2.getSequenceIndexAt(2), 3);
        assertEquals(local2.getSequenceIndexAt(3), 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetSequenceIndexAtOutOfBounds() {
        global.getSequenceIndexAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetSequenceIndexAtOutOfBounds2() {
        global.getSequenceIndexAt(8);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetSequenceIndexAtOutOfBounds3() {
        local.getSequenceIndexAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetSequenceIndexAtOutOfBounds4() {
        local.getSequenceIndexAt(9);
    }

    @Test
    public void testGetStart() {
        assertEquals(global.getStart().getPosition(), Integer.valueOf(2));
        assertEquals(local.getStart().getPosition(), Integer.valueOf(1));
        assertEquals(local2.getStart().getPosition(), Integer.valueOf(1));
    }

    @Test
    public void testIsCircular() {
        assertFalse(global.isCircular());
        assertFalse(local.isCircular());
        assertFalse(local2.isCircular());
    }

    @Test
    public void testCountCompounds() {
        assertEquals(global.countCompounds(cs.getCompoundForString("A"), cs.getCompoundForString("N"),
                cs.getCompoundForString("A"), cs.getCompoundForString("E"), cs.getCompoundForString("D")), 3);
        assertEquals(local.countCompounds(cs.getCompoundForString("A"), cs.getCompoundForString("N"),
                cs.getCompoundForString("A"), cs.getCompoundForString("E"), cs.getCompoundForString("D")), 1);
        assertEquals(local2.countCompounds(cs.getCompoundForString("A"), cs.getCompoundForString("N"),
                cs.getCompoundForString("A"), cs.getCompoundForString("E"), cs.getCompoundForString("D")), 2);
    }

    @Test
    public void testGetAccession() {
        assertNull(global.getAccession());
        assertNull(local.getAccession());
        assertNull(local2.getAccession());
    }

    @Test
    public void testGetAsList() {
        assertArrayEquals(global.getAsList().toArray(new AminoAcidCompound[7]), new AminoAcidCompound[] {
            cs.getCompoundForString("-"), cs.getCompoundForString("A"), cs.getCompoundForString("R"),
            cs.getCompoundForString("-"), cs.getCompoundForString("N"), cs.getCompoundForString("D"),
            cs.getCompoundForString("-")});
        assertArrayEquals(local.getAsList().toArray(new AminoAcidCompound[8]), new AminoAcidCompound[] {
            cs.getCompoundForString("E"), cs.getCompoundForString("Q"), cs.getCompoundForString("-"),
            cs.getCompoundForString("-"), cs.getCompoundForString("G"), cs.getCompoundForString("-"),
            cs.getCompoundForString("H"), cs.getCompoundForString("I")});
        assertArrayEquals(local2.getAsList().toArray(new AminoAcidCompound[3]), new AminoAcidCompound[] {
            cs.getCompoundForString("R"), cs.getCompoundForString("N"), cs.getCompoundForString("D")});
    }

    @Test
    public void testGetCompoundAt() {
        assertEquals(global.getCompoundAt(1), cs.getCompoundForString("-"));
        assertEquals(global.getCompoundAt(2), cs.getCompoundForString("A"));
        assertEquals(global.getCompoundAt(3), cs.getCompoundForString("R"));
        assertEquals(global.getCompoundAt(4), cs.getCompoundForString("-"));
        assertEquals(global.getCompoundAt(5), cs.getCompoundForString("N"));
        assertEquals(global.getCompoundAt(6), cs.getCompoundForString("D"));
        assertEquals(global.getCompoundAt(7), cs.getCompoundForString("-"));
        assertEquals(global.getCompoundAt(1), cs.getCompoundForString("-"));
        assertEquals(local.getCompoundAt(1), cs.getCompoundForString("E"));
        assertEquals(local.getCompoundAt(2), cs.getCompoundForString("Q"));
        assertEquals(local.getCompoundAt(3), cs.getCompoundForString("-"));
        assertEquals(local.getCompoundAt(4), cs.getCompoundForString("-"));
        assertEquals(local.getCompoundAt(5), cs.getCompoundForString("G"));
        assertEquals(local.getCompoundAt(6), cs.getCompoundForString("-"));
        assertEquals(local.getCompoundAt(7), cs.getCompoundForString("H"));
        assertEquals(local.getCompoundAt(8), cs.getCompoundForString("I"));
        assertEquals(local2.getCompoundAt(1), cs.getCompoundForString("R"));
        assertEquals(local2.getCompoundAt(2), cs.getCompoundForString("N"));
        assertEquals(local2.getCompoundAt(3), cs.getCompoundForString("D"));
    }

    @Test
    public void testGetCompoundSet() {
        assertEquals(global.getCompoundSet(), cs);
        assertEquals(local.getCompoundSet(), cs);
        assertEquals(local2.getCompoundSet(), cs);
    }

    @Test
    public void testGetIndexOf() {
        assertEquals(global.getIndexOf(cs.getCompoundForString("R")), 3);
        assertEquals(global.getIndexOf(cs.getCompoundForString("-")), 1);
        assertEquals(local.getIndexOf(cs.getCompoundForString("G")), 5);
        assertEquals(local.getIndexOf(cs.getCompoundForString("-")), 3);
        assertEquals(local2.getIndexOf(cs.getCompoundForString("N")), 2);
        assertEquals(local2.getIndexOf(cs.getCompoundForString("-")), -1);
    }

    @Test
    public void testGetLastIndexOf() {
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("R")), 3);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("-")), 7);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("G")), 5);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("-")), 6);
        assertEquals(local2.getLastIndexOf(cs.getCompoundForString("N")), 2);
        assertEquals(local2.getLastIndexOf(cs.getCompoundForString("-")), -1);
    }

    @Test
    public void testGetLength() {
        assertEquals(global.getLength(), 7);
        assertEquals(local.getLength(), 8);
        assertEquals(local2.getLength(), 3);
    }

    @Test
    public void testGetSequenceAsString() {
        assertEquals(global.getSequenceAsString(), "-AR-ND-");
        assertEquals(local.getSequenceAsString(), "EQ--G-HI");
        assertEquals(local2.getSequenceAsString(), "RND");
    }

    @Test
    public void testGetSequenceAsStringIntegerIntegerStrand() {
        assertEquals(global.getSubSequence(2, 5).getSequenceAsString(), "AR-N");
        assertEquals(local.getSubSequence(2, 6).getSequenceAsString(), "Q--G-");
        assertEquals(local2.getSubSequence(2, 3).getSequenceAsString(), "ND");
    }

    @Ignore // TODO SimpleAlignedSequence.getSubSequence(Integer, Integer)
    @Test
    public void testGetSubSequence() {
        fail("Not yet implemented");
    }

    @Test
    public void testIterator() {
        for (AminoAcidCompound c : global) {
            assertNotNull(cs.getStringForCompound(c));
        }
        for (AminoAcidCompound c : local) {
            assertNotNull(cs.getStringForCompound(c));
        }
        for (AminoAcidCompound c : local2) {
            assertNotNull(cs.getStringForCompound(c));
        }
    }

    @Test
    public void testToString() {
        assertEquals(global.toString(), "-AR-ND-");
        assertEquals(local.toString(), "EQ--G-HI");
        assertEquals(local2.toString(), "RND");
    }

}
