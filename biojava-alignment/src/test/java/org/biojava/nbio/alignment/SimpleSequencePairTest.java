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

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class SimpleSequencePairTest {

    private ProteinSequence query, target;
    private SequencePair<ProteinSequence, AminoAcidCompound> global, local;

    @Before
    public void setup() throws CompoundNotFoundException { 
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        global = new SimpleSequencePair<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.COMPOUND, Step.GAP}), Arrays.asList(new Step[] {
                Step.GAP, Step.COMPOUND, Step.GAP, Step.COMPOUND, Step.COMPOUND}));
        local = new SimpleSequencePair<ProteinSequence, AminoAcidCompound>(query, target, Arrays.asList(new Step[] {
                Step.COMPOUND, Step.COMPOUND, Step.COMPOUND}), 1, 0, Arrays.asList(new Step[] { Step.COMPOUND,
                Step.GAP, Step.COMPOUND}), 0, 1);
    }

    @Test
    public void testGetCompoundInQueryAt() {
        assertEquals(global.getCompoundInQueryAt(1).getShortName(), "A");
        assertEquals(global.getCompoundInQueryAt(2).getShortName(), "R");
        assertEquals(global.getCompoundInQueryAt(3).getShortName(), "N");
        assertEquals(global.getCompoundInQueryAt(4).getShortName(), "D");
        assertEquals(global.getCompoundInQueryAt(5).getShortName(), "-");
        assertEquals(local.getCompoundInQueryAt(1).getShortName(), "R");
        assertEquals(local.getCompoundInQueryAt(2).getShortName(), "N");
        assertEquals(local.getCompoundInQueryAt(3).getShortName(), "D");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInQueryAtOutOfBounds() {
        global.getCompoundInQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInQueryAtOutOfBounds2() {
        global.getCompoundInQueryAt(6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInQueryAtOutOfBounds3() {
        local.getCompoundInQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInQueryAtOutOfBounds4() {
        local.getCompoundInQueryAt(4);
    }

    @Test
    public void testGetCompoundInTargetAt() {
        assertEquals(global.getCompoundInTargetAt(1).getShortName(), "-");
        assertEquals(global.getCompoundInTargetAt(2).getShortName(), "R");
        assertEquals(global.getCompoundInTargetAt(3).getShortName(), "-");
        assertEquals(global.getCompoundInTargetAt(4).getShortName(), "D");
        assertEquals(global.getCompoundInTargetAt(5).getShortName(), "G");
        assertEquals(local.getCompoundInTargetAt(1).getShortName(), "R");
        assertEquals(local.getCompoundInTargetAt(2).getShortName(), "-");
        assertEquals(local.getCompoundInTargetAt(3).getShortName(), "D");
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInTargetAtOutOfBounds() {
        global.getCompoundInTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInTargetAtOutOfBounds2() {
        global.getCompoundInTargetAt(6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInTargetAtOutOfBounds3() {
        local.getCompoundInTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetCompoundInTargetAtOutOfBounds4() {
        local.getCompoundInTargetAt(4);
    }

    @Test
    public void testGetIndexInQueryAt() {
        assertEquals(global.getIndexInQueryAt(1), 1);
        assertEquals(global.getIndexInQueryAt(2), 2);
        assertEquals(global.getIndexInQueryAt(3), 3);
        assertEquals(global.getIndexInQueryAt(4), 4);
        assertEquals(global.getIndexInQueryAt(5), 4);
        assertEquals(local.getIndexInQueryAt(1), 2);
        assertEquals(local.getIndexInQueryAt(2), 3);
        assertEquals(local.getIndexInQueryAt(3), 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryAtOutOfBounds() {
        global.getIndexInQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryAtOutOfBounds2() {
        global.getIndexInQueryAt(6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryAtOutOfBounds3() {
        local.getIndexInQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryAtOutOfBounds4() {
        local.getIndexInQueryAt(4);
    }

    @Test
    public void testGetIndexInQueryForTargetAt() {
        assertEquals(global.getIndexInQueryForTargetAt(1), 2);
        assertEquals(global.getIndexInQueryForTargetAt(2), 4);
        assertEquals(global.getIndexInQueryForTargetAt(3), 4);
        assertEquals(local.getIndexInQueryForTargetAt(1), 2);
        assertEquals(local.getIndexInQueryForTargetAt(2), 4);
        assertEquals(local.getIndexInQueryForTargetAt(3), 4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryForTargetAtOutOfBounds() {
        global.getIndexInQueryForTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryForTargetAtOutOfBounds2() {
        global.getIndexInQueryForTargetAt(4);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryForTargetAtOutOfBounds3() {
        local.getIndexInQueryForTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInQueryForTargetAtOutOfBounds4() {
        local.getIndexInQueryForTargetAt(4);
    }

    @Test
    public void testGetIndexInTargetAt() {
        assertEquals(global.getIndexInTargetAt(1), 1);
        assertEquals(global.getIndexInTargetAt(2), 1);
        assertEquals(global.getIndexInTargetAt(3), 1);
        assertEquals(global.getIndexInTargetAt(4), 2);
        assertEquals(global.getIndexInTargetAt(5), 3);
        assertEquals(local.getIndexInTargetAt(1), 1);
        assertEquals(local.getIndexInTargetAt(2), 1);
        assertEquals(local.getIndexInTargetAt(3), 2);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetAtOutOfBounds() {
        global.getIndexInTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetAtOutOfBounds2() {
        global.getIndexInTargetAt(6);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetAtOutOfBounds3() {
        local.getIndexInTargetAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetAtOutOfBounds4() {
        local.getIndexInTargetAt(4);
    }

    @Test
    public void testGetIndexInTargetForQueryAt() {
        assertEquals(global.getIndexInTargetForQueryAt(1), 1);
        assertEquals(global.getIndexInTargetForQueryAt(2), 1);
        assertEquals(global.getIndexInTargetForQueryAt(3), 1);
        assertEquals(global.getIndexInTargetForQueryAt(4), 2);
        assertEquals(local.getIndexInTargetForQueryAt(1), 1);
        assertEquals(local.getIndexInTargetForQueryAt(2), 1);
        assertEquals(local.getIndexInTargetForQueryAt(3), 1);
        assertEquals(local.getIndexInTargetForQueryAt(4), 2);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetForQueryAtOutOfBounds() {
        global.getIndexInTargetForQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetForQueryAtOutOfBounds2() {
        global.getIndexInTargetForQueryAt(5);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetForQueryAtOutOfBounds3() {
        local.getIndexInTargetForQueryAt(0);
    }

    @Test(expected=IndexOutOfBoundsException.class)
    public void testGetIndexInTargetForQueryAtOutOfBounds4() {
        local.getIndexInTargetForQueryAt(5);
    }

    @Test
    public void testGetNumIdenticals() {
        assertEquals(global.getNumIdenticals(), 2);
        assertEquals(local.getNumIdenticals(), 2);
    }

    @Test
    public void testGetNumSimilars() {
        assertEquals(global.getNumSimilars(), 2);
        assertEquals(local.getNumSimilars(), 2);
    }

    @Test
    public void testGetQuery() {
        assertEquals(global.getQuery().getOriginalSequence(), query);
        assertEquals(local.getQuery().getOriginalSequence(), query);
    }

    @Test
    public void testGetTarget() {
        assertEquals(global.getTarget().getOriginalSequence(), target);
        assertEquals(local.getTarget().getOriginalSequence(), target);
    }

}
