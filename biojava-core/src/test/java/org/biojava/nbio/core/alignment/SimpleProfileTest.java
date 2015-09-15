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
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.Profile.StringFormat;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class SimpleProfileTest {

    private ProteinSequence query, target;
    private Profile<ProteinSequence, AminoAcidCompound> global, local, single;

    @Before
    public void setup() throws CompoundNotFoundException { 
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        query.setAccession(new AccessionID("Query"));
        target.setAccession(new AccessionID("Target"));
        global = new SimpleProfile<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}), 0, 0, Arrays.asList(
                new Step[] {Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND}), 0, 0);
        local = new SimpleProfile<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND}), 1, 0, Arrays.asList(new Step[] { Step.COMPOUND,
                Step.GAP, Step.COMPOUND}), 0, 1);
        single = new SimpleProfile<ProteinSequence, AminoAcidCompound>(query);
    }

    @Test(expected=IllegalArgumentException.class)
    public void testSimpleProfile() {
        new SimpleProfile<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}), 0, 0, Arrays.asList(
                new Step[] {Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND}), 0, 0);
    }

    @Test
    public void testGetAlignedSequenceInt() {
        assertEquals(global.getAlignedSequence(1).toString(), "ARND-");
        assertEquals(global.getAlignedSequence(2).toString(), "-R-DG");
        assertEquals(local.getAlignedSequence(1).toString(), "RND");
        assertEquals(local.getAlignedSequence(2).toString(), "R-D");
        assertEquals(single.getAlignedSequence(1).toString(), "ARND");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds() {
        global.getAlignedSequence(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds2() {
        global.getAlignedSequence(3);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds3() {
        local.getAlignedSequence(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds4() {
        local.getAlignedSequence(3);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds5() {
        single.getAlignedSequence(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetAlignedSequenceIntOutOfBounds6() {
        single.getAlignedSequence(2);
    }

    @Test
    public void testGetAlignedSequenceS() throws CompoundNotFoundException { 
        assertEquals(global.getAlignedSequence(query).toString(), "ARND-");
        assertEquals(global.getAlignedSequence(target).toString(), "-R-DG");
        assertNull(global.getAlignedSequence(new ProteinSequence("AR")));
        assertEquals(local.getAlignedSequence(query).toString(), "RND");
        assertEquals(local.getAlignedSequence(target).toString(), "R-D");
        assertNull(local.getAlignedSequence(new ProteinSequence("AR")));
        assertEquals(single.getAlignedSequence(query).toString(), "ARND");
        assertNull(single.getAlignedSequence(target));
    }

    @Test
    public void testGetAlignedSequences() {
        List<AlignedSequence<ProteinSequence, AminoAcidCompound>> list = global.getAlignedSequences();
        assertEquals(list.size(), 2);
        assertEquals(list.get(0).toString(), "ARND-");
        assertEquals(list.get(1).toString(), "-R-DG");
        list = local.getAlignedSequences();
        assertEquals(list.size(), 2);
        assertEquals(list.get(0).toString(), "RND");
        assertEquals(list.get(1).toString(), "R-D");
        list = single.getAlignedSequences();
        assertEquals(list.size(), 1);
        assertEquals(list.get(0).toString(), "ARND");
    }

    @Test
    public void testGetAlignedSequencesIntArray() {
        List<AlignedSequence<ProteinSequence, AminoAcidCompound>> list = global.getAlignedSequences(2, 1, 2);
        assertEquals(list.size(), 3);
        assertEquals(list.get(0).toString(), "-R-DG");
        assertEquals(list.get(1).toString(), "ARND-");
        assertEquals(list.get(2).toString(), "-R-DG");
        list = local.getAlignedSequences(2, 2, 1);
        assertEquals(list.size(), 3);
        assertEquals(list.get(0).toString(), "R-D");
        assertEquals(list.get(1).toString(), "R-D");
        assertEquals(list.get(2).toString(), "RND");
        list = single.getAlignedSequences(1, 1);
        assertEquals(list.size(), 2);
        assertEquals(list.get(0).toString(), "ARND");
        assertEquals(list.get(1).toString(), "ARND");
    }

    @Test
    public void testGetAlignedSequencesSArray() {
        List<AlignedSequence<ProteinSequence, AminoAcidCompound>> list = global.getAlignedSequences(query, query,
                target);
        assertEquals(list.size(), 3);
        assertEquals(list.get(0).toString(), "ARND-");
        assertEquals(list.get(1).toString(), "ARND-");
        assertEquals(list.get(2).toString(), "-R-DG");
        list = local.getAlignedSequences(target, query, target);
        assertEquals(list.size(), 3);
        assertEquals(list.get(0).toString(), "R-D");
        assertEquals(list.get(1).toString(), "RND");
        assertEquals(list.get(2).toString(), "R-D");
        list = single.getAlignedSequences(query, query);
        assertEquals(list.size(), 2);
        assertEquals(list.get(0).toString(), "ARND");
        assertEquals(list.get(1).toString(), "ARND");
    }

    @Test
    public void testGetCompoundAtIntInt() {
        assertEquals(global.getCompoundAt(1, 4).getShortName(), "D");
        assertEquals(global.getCompoundAt(2, 3).getShortName(), "-");
        assertEquals(local.getCompoundAt(1, 1).getShortName(), "R");
        assertEquals(local.getCompoundAt(2, 2).getShortName(), "-");
        assertEquals(single.getCompoundAt(1, 3).getShortName(), "N");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds() {
        global.getCompoundAt(0, 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds2() {
        global.getCompoundAt(3, 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds3() {
        global.getCompoundAt(1, 0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds4() {
        global.getCompoundAt(2, 6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds5() {
        local.getCompoundAt(0, 2);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds6() {
        local.getCompoundAt(3, 2);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds7() {
        local.getCompoundAt(1, 0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds8() {
        local.getCompoundAt(2, 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtIntIntOutOfBounds9() {
        single.getCompoundAt(1, 0);
    }

    @Test
    public void testGetCompoundAtSInt() throws CompoundNotFoundException { 
        assertEquals(global.getCompoundAt(query, 2).getShortName(), "R");
        assertEquals(global.getCompoundAt(target, 5).getShortName(), "G");
        assertNull(global.getCompoundAt(new ProteinSequence("AR"), 3));
        assertEquals(local.getCompoundAt(query, 2).getShortName(), "N");
        assertEquals(local.getCompoundAt(target, 3).getShortName(), "D");
        assertNull(local.getCompoundAt(new ProteinSequence("AR"), 3));
        assertEquals(single.getCompoundAt(query, 2).getShortName(), "R");
        assertNull(single.getCompoundAt(target, 3));
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds() {
        global.getCompoundAt(query, 0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds2() {
        global.getCompoundAt(target, 6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds3() {
        local.getCompoundAt(target, 0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds4() {
        local.getCompoundAt(query, 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundAtSIntOutOfBounds5() {
        single.getCompoundAt(query, 0);
    }

    @Test
    public void testGetCompoundSet() {
        assertEquals(global.getCompoundSet(), AminoAcidCompoundSet.getAminoAcidCompoundSet());
        assertEquals(local.getCompoundSet(), AminoAcidCompoundSet.getAminoAcidCompoundSet());
        assertEquals(single.getCompoundSet(), AminoAcidCompoundSet.getAminoAcidCompoundSet());
    }

    @Test
    public void testGetCompoundsAt() {
        List<AminoAcidCompound> column = global.getCompoundsAt(5);
        assertEquals(column.size(), 2);
        assertEquals(column.get(0).getShortName(), "-");
        assertEquals(column.get(1).getShortName(), "G");
        column = local.getCompoundsAt(2);
        assertEquals(column.size(), 2);
        assertEquals(column.get(0).getShortName(), "N");
        assertEquals(column.get(1).getShortName(), "-");
        column = single.getCompoundsAt(2);
        assertEquals(column.size(), 1);
        assertEquals(column.get(0).getShortName(), "R");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds() {
        global.getCompoundsAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds2() {
        global.getCompoundsAt(6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds3() {
        local.getCompoundsAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds4() {
        local.getCompoundsAt(4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds5() {
        single.getCompoundsAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundsAtOutOfBounds6() {
        single.getCompoundsAt(5);
    }

    @Test
    public void testGetIndexOf() {
        AminoAcidCompoundSet cs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
        assertEquals(global.getIndexOf(cs.getCompoundForString("A")), 1);
        assertEquals(global.getIndexOf(cs.getCompoundForString("R")), 2);
        assertEquals(global.getIndexOf(cs.getCompoundForString("N")), 3);
        assertEquals(global.getIndexOf(cs.getCompoundForString("D")), 4);
        assertEquals(global.getIndexOf(cs.getCompoundForString("G")), 5);
        assertEquals(global.getIndexOf(cs.getCompoundForString("-")), 1);
        assertEquals(global.getIndexOf(cs.getCompoundForString("E")), -1);
        assertEquals(local.getIndexOf(cs.getCompoundForString("R")), 1);
        assertEquals(local.getIndexOf(cs.getCompoundForString("N")), 2);
        assertEquals(local.getIndexOf(cs.getCompoundForString("D")), 3);
        assertEquals(local.getIndexOf(cs.getCompoundForString("-")), 2);
        assertEquals(local.getIndexOf(cs.getCompoundForString("K")), -1);
        assertEquals(single.getIndexOf(cs.getCompoundForString("A")), 1);
        assertEquals(single.getIndexOf(cs.getCompoundForString("R")), 2);
        assertEquals(single.getIndexOf(cs.getCompoundForString("N")), 3);
        assertEquals(single.getIndexOf(cs.getCompoundForString("D")), 4);
        assertEquals(single.getIndexOf(cs.getCompoundForString("G")), -1);
    }

    @Test
    public void testGetIndicesAt() {
        assertArrayEquals(global.getIndicesAt(1), new int[] {1, 1});
        assertArrayEquals(global.getIndicesAt(2), new int[] {2, 1});
        assertArrayEquals(global.getIndicesAt(3), new int[] {3, 1});
        assertArrayEquals(global.getIndicesAt(4), new int[] {4, 2});
        assertArrayEquals(global.getIndicesAt(5), new int[] {4, 3});
        assertArrayEquals(local.getIndicesAt(1), new int[] {2, 1});
        assertArrayEquals(local.getIndicesAt(2), new int[] {3, 1});
        assertArrayEquals(local.getIndicesAt(3), new int[] {4, 2});
        assertArrayEquals(single.getIndicesAt(1), new int[] {1});
        assertArrayEquals(single.getIndicesAt(2), new int[] {2});
        assertArrayEquals(single.getIndicesAt(3), new int[] {3});
        assertArrayEquals(single.getIndicesAt(4), new int[] {4});
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds() {
        global.getIndicesAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds2() {
        global.getIndicesAt(6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds3() {
        local.getIndicesAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds4() {
        local.getIndicesAt(4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds5() {
        single.getIndicesAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndicesAtOutOfBounds6() {
        single.getIndicesAt(5);
    }

    @Test
    public void testGetLastIndexOf() {
        AminoAcidCompoundSet cs = AminoAcidCompoundSet.getAminoAcidCompoundSet();
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("A")), 1);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("R")), 2);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("N")), 3);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("D")), 4);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("G")), 5);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("-")), 5);
        assertEquals(global.getLastIndexOf(cs.getCompoundForString("E")), -1);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("R")), 1);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("N")), 2);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("D")), 3);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("-")), 2);
        assertEquals(local.getLastIndexOf(cs.getCompoundForString("K")), -1);
        assertEquals(single.getLastIndexOf(cs.getCompoundForString("A")), 1);
        assertEquals(single.getLastIndexOf(cs.getCompoundForString("R")), 2);
        assertEquals(single.getLastIndexOf(cs.getCompoundForString("N")), 3);
        assertEquals(single.getLastIndexOf(cs.getCompoundForString("D")), 4);
        assertEquals(single.getLastIndexOf(cs.getCompoundForString("G")), -1);
    }

    @Test
    public void testGetLength() {
        assertEquals(global.getLength(), 5);
        assertEquals(local.getLength(), 3);
        assertEquals(single.getLength(), 4);
    }

    @Test
    public void testGetOriginalSequences() {
        List<ProteinSequence> list = global.getOriginalSequences();
        assertEquals(list.size(), 2);
        assertEquals(list.get(0), query);
        assertEquals(list.get(1), target);
        list = local.getOriginalSequences();
        assertEquals(list.size(), 2);
        assertEquals(list.get(0), query);
        assertEquals(list.get(1), target);
        list = single.getOriginalSequences();
        assertEquals(list.size(), 1);
        assertEquals(list.get(0), query);
    }

    @Test
    public void testGetSize() {
        assertEquals(global.getSize(), 2);
        assertEquals(local.getSize(), 2);
        assertEquals(single.getSize(), 1);
    }

    @Ignore // TODO SimpleProfile.getSubProfile(Location)
    @Test
    public void testGetSubProfile() {
        fail("Not yet implemented");
    }

    @Test
    public void testIsCircular() {
        assertFalse(global.isCircular());
        assertFalse(local.isCircular());
        assertFalse(single.isCircular());
    }

    @Test
    public void testToStringInt() {
    
    	
    
    	
        assertEquals(global.toString(3), String.format(
                "          1 3%n" +
                "Query   1 ARN 3%n" +                
                "           | %n"+
                "Target  1 -R- 1%n" +
                "%n" +
                "          4 5%n" +
                "Query   4 D- 4%n" +
                "          | %n"+
                "Target  2 DG 3%n"));
        assertEquals(local.toString(4), String.format(
                "          1 3%n" +
                "Query   2 RND 4%n" +
                "          | |%n"+
                "Target  1 R-D 2%n"));
        assertEquals(single.toString(4), String.format(
                "         1  4%n" +
                "Query  1 ARND 4%n" ));
    }

    @Test
    public void testToStringFormatted() {
       
        assertEquals(global.toString(StringFormat.ALN), String.format(
                "CLUSTAL W MSA from BioJava%n%n" +
                "Query     ARND- 4%n" +
                "           | | %n"+
                "Target    -R-DG 3%n"));
        assertEquals(local.toString(StringFormat.FASTA), String.format(
                ">Query%n" +
                "RND%n" +
                ">Target%n" +
                "R-D%n"));
        assertEquals(single.toString(StringFormat.MSF), String.format(
                "MSA from BioJava%n%n" +
                " MSF: 4  Type: P  Check: 735 ..%n%n" +
                " Name: Query  Len: 4  Check:  735  Weight: 1.0%n" +
                "%n//%n%n" +
                "Query ARND%n"));
    }

    @Test
    public void testToString() {
        assertEquals(global.toString(), String.format("ARND-%n-R-DG%n"));
        assertEquals(local.toString(), String.format("RND%nR-D%n"));
        assertEquals(single.toString(), String.format("ARND%n"));
    }

    @Test
    public void testIterator() {
        for (AlignedSequence<ProteinSequence, AminoAcidCompound> s : global) {
            assertEquals(s.toString().length(), 5);
        }
        for (AlignedSequence<ProteinSequence, AminoAcidCompound> s : local) {
            assertEquals(s.toString().length(), 3);
        }
        for (AlignedSequence<ProteinSequence, AminoAcidCompound> s : single) {
            assertEquals(s.toString().length(), 4);
        }
    }

}
