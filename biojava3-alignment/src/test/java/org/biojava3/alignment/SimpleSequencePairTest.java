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
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class SimpleSequencePairTest {

    private ProteinSequence query, target;
    private SequencePair<ProteinSequence, AminoAcidCompound> pair;

    @Before
    public void setup() {
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        pair = new SimpleSequencePair<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}), Arrays.asList(new Step[] {
                Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND}));
    }

    @Test
    public void testGetCompoundInQueryAt() {
        assertEquals(pair.getCompoundInQueryAt(1).getShortName(), "A");
        assertEquals(pair.getCompoundInQueryAt(2).getShortName(), "R");
        assertEquals(pair.getCompoundInQueryAt(3).getShortName(), "N");
        assertEquals(pair.getCompoundInQueryAt(4).getShortName(), "D");
        assertEquals(pair.getCompoundInQueryAt(5).getShortName(), "-");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInQueryAtOutOfBounds() {
        pair.getCompoundInQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInQueryAtOutOfBounds2() {
        pair.getCompoundInQueryAt(6);
    }

    @Test
    public void testGetCompoundInTargetAt() {
        assertEquals(pair.getCompoundInTargetAt(1).getShortName(), "-");
        assertEquals(pair.getCompoundInTargetAt(2).getShortName(), "R");
        assertEquals(pair.getCompoundInTargetAt(3).getShortName(), "-");
        assertEquals(pair.getCompoundInTargetAt(4).getShortName(), "D");
        assertEquals(pair.getCompoundInTargetAt(5).getShortName(), "G");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInTargetAtOutOfBounds() {
        pair.getCompoundInTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInTargetAtOutOfBounds2() {
        pair.getCompoundInTargetAt(6);
    }

    @Test
    public void testGetIndexInQueryAt() {
        assertEquals(pair.getIndexInQueryAt(1), 1);
        assertEquals(pair.getIndexInQueryAt(2), 2);
        assertEquals(pair.getIndexInQueryAt(3), 3);
        assertEquals(pair.getIndexInQueryAt(4), 4);
        assertEquals(pair.getIndexInQueryAt(5), 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryAtOutOfBounds() {
        pair.getIndexInQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryAtOutOfBounds2() {
        pair.getIndexInQueryAt(6);
    }

    @Test
    public void testGetIndexInQueryForTargetAt() {
        assertEquals(pair.getIndexInQueryForTargetAt(1), 2);
        assertEquals(pair.getIndexInQueryForTargetAt(2), 4);
        assertEquals(pair.getIndexInQueryForTargetAt(3), 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryForTargetAtOutOfBounds() {
        pair.getIndexInQueryForTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryForTargetAtOutOfBounds2() {
        pair.getIndexInQueryForTargetAt(4);
    }

    @Test
    public void testGetIndexInTargetAt() {
        assertEquals(pair.getIndexInTargetAt(1), 1);
        assertEquals(pair.getIndexInTargetAt(2), 1);
        assertEquals(pair.getIndexInTargetAt(3), 1);
        assertEquals(pair.getIndexInTargetAt(4), 2);
        assertEquals(pair.getIndexInTargetAt(5), 3);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetAtOutOfBounds() {
        pair.getIndexInTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetAtOutOfBounds2() {
        pair.getIndexInTargetAt(6);
    }

    @Test
    public void testGetIndexInTargetForQueryAt() {
        assertEquals(pair.getIndexInTargetForQueryAt(1), 1);
        assertEquals(pair.getIndexInTargetForQueryAt(2), 1);
        assertEquals(pair.getIndexInTargetForQueryAt(3), 1);
        assertEquals(pair.getIndexInTargetForQueryAt(4), 2);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetForQueryAtOutOfBounds() {
        pair.getIndexInTargetForQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetForQueryAtOutOfBounds2() {
        pair.getIndexInTargetForQueryAt(5);
    }

    @Test
    public void testGetNumIdenticals() {
        assertEquals(pair.getNumIdenticals(), 2);
    }

    @Ignore // TODO implement AminoAcidCompoundSet.compoundsEquivalent(C, C)
    @Test
    public void testGetNumSimilars() {
        assertEquals(pair.getNumSimilars(), 2);
    }

    @Test
    public void testGetQuery() {
        assertEquals(pair.getQuery().getOriginalSequence(), query);
    }

    @Test
    public void testGetTarget() {
        assertEquals(pair.getTarget().getOriginalSequence(), target);
    }

}
