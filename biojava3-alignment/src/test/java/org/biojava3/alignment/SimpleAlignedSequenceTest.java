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

import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.location.SimpleLocation;
import org.biojava3.core.sequence.location.template.Location;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class SimpleAlignedSequenceTest {

    private SimpleAlignedSequence<AminoAcidCompound> as;
    private AminoAcidCompoundSet cs;

    @Before
    public void setup() {
        Step[] steps = { Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.GAP };
        as = new SimpleAlignedSequence<AminoAcidCompound>(new ProteinSequence("ARND"), Arrays.asList(steps));
        cs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
    }

    @Test
    public void testSimpleAlignedSequence() {
        Step[] steps = { Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND, Step.GAP };
        SimpleAlignedSequence<AminoAcidCompound> as =
                new SimpleAlignedSequence<AminoAcidCompound>(new ProteinSequence("ARND"), Arrays.asList(steps));
        assertEquals(as.getLength(), 7);
    }

    @Test
    public void testGetAlignmentIndexAt() {
        assertEquals(as.getAlignmentIndexAt(2), 3);
        assertEquals(as.getAlignmentIndexAt(4), 6);
    }

    @Test
    public void testGetEnd() {
        assertEquals(as.getEnd(), 6);
    }

    @Test
    public void testGetLocationInAlignment() {
        Location[] sublocations = { new SimpleLocation(2, 3, Strand.UNDEFINED),
                new SimpleLocation(5, 6, Strand.UNDEFINED) };
        assertEquals(as.getLocationInAlignment(), new SimpleLocation(2, 6, Strand.UNDEFINED,
                Arrays.asList(sublocations)));
    }

    @Test
    public void testGetNumGaps() {
        assertEquals(as.getNumGaps(), 1);
    }

    @Ignore // TODO implement ProteinSequence.equals(Object)
    @Test
    public void testGetOriginalSequence() {
        assertEquals(as.getOriginalSequence(), new ProteinSequence("ARND"));
    }

    @Test
    public void testGetOverlapCount() {
        assertEquals(as.getOverlapCount(), 1);
    }

    @Test
    public void testGetSequenceIndexAt() {
        assertEquals(as.getSequenceIndexAt(1), 1);
        assertEquals(as.getSequenceIndexAt(2), 1);
        assertEquals(as.getSequenceIndexAt(3), 2);
        assertEquals(as.getSequenceIndexAt(4), 2);
        assertEquals(as.getSequenceIndexAt(5), 3);
        assertEquals(as.getSequenceIndexAt(6), 4);
        assertEquals(as.getSequenceIndexAt(7), 4);
    }

    @Test
    public void testGetStart() {
        assertEquals(as.getStart(), 2);
    }

    @Test
    public void testIsCircular() {
        assertFalse(as.isCircular());
    }

    @Test
    public void testCountCompounds() {
        assertEquals(as.countCompounds(cs.getCompoundForString("A"), cs.getCompoundForString("N")), 2);
    }

    @Test
    public void testGetAccession() {
        assertNull(as.getAccession());
    }

    @Ignore // TODO fix isGap()
    @Test
    public void testGetAsList() {
        AminoAcidCompound[] compounds = { cs.getCompoundForString("-"), cs.getCompoundForString("A"),
                cs.getCompoundForString("R"), cs.getCompoundForString("-"), cs.getCompoundForString("N"),
                cs.getCompoundForString("D"), cs.getCompoundForString("-") }, list = new AminoAcidCompound[7];
        assertArrayEquals(as.getAsList().toArray(list), compounds);
    }

    @Test
    public void testGetCompoundAt() {
        assertEquals(as.getCompoundAt(1), cs.getCompoundForString("-"));
        assertEquals(as.getCompoundAt(2), cs.getCompoundForString("A"));
        assertEquals(as.getCompoundAt(3), cs.getCompoundForString("R"));
        assertEquals(as.getCompoundAt(4), cs.getCompoundForString("-"));
        assertEquals(as.getCompoundAt(5), cs.getCompoundForString("N"));
        assertEquals(as.getCompoundAt(6), cs.getCompoundForString("D"));
        assertEquals(as.getCompoundAt(7), cs.getCompoundForString("-"));
    }

    @Test
    public void testGetCompoundSet() {
        assertEquals(as.getCompoundSet(), cs);
    }

    @Test
    public void testGetIndexOf() {
        assertEquals(as.getIndexOf(cs.getCompoundForString("R")), 3);
    }

    @Test
    public void testGetLastIndexOf() {
        assertEquals(as.getLastIndexOf(cs.getCompoundForString("R")), 3);
    }

    @Test
    public void testGetLength() {
        assertEquals(as.getLength(), 7);
    }

    @Ignore // TODO fix isGap()
    @Test
    public void testGetSequenceAsString() {
        assertEquals(as.getSequenceAsString(), "-AR-ND-");
    }

    @Ignore // TODO fix isGap()
    @Test
    public void testGetSequenceAsStringIntegerIntegerStrand() {
        assertEquals(as.getSequenceAsString(2, 5, Strand.UNDEFINED), "AR-N");
    }

    @Ignore // TODO SimpleAlignedSequence.getSubSequence(Integer, Integer)
    @Test
    public void testGetSubSequence() {
        fail("Not yet implemented");
    }

    @Test
    public void testIterator() {
        for (AminoAcidCompound c : as) {
            assertNotNull(cs.getStringForCompound(c));
        }
    }

    @Ignore // TODO fix isGap()
    @Test
    public void testToString() {
        assertEquals(as.toString(), "-AR-ND-");
    }

}
