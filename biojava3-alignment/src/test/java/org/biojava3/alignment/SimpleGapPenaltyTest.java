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
 * Created on Jun 9, 2010
 * Author: Mark 
 *
 */

package org.biojava3.alignment;

import static org.junit.Assert.*;

import org.biojava3.alignment.template.GapPenalty;
import org.junit.Test;

public class SimpleGapPenaltyTest {

    @Test
    public void testSimpleGapPenalty() {
        GapPenalty gaps = new SimpleGapPenalty();
        short gop = SimpleGapPenalty.Defaults.getOpenPenalty();
        short gep = SimpleGapPenalty.Defaults.getExtensionPenalty();
        assertEquals(gaps.getOpenPenalty(), gop);
        assertEquals(gaps.getExtensionPenalty(), gep);
        assertEquals(gaps.getType(), new SimpleGapPenalty().getType());
    }

    @Test
    public void testSimpleGapPenaltyShortShort() {
        short gop = 10, gep = 4;
        GapPenalty gaps = new SimpleGapPenalty(gop, gep);
        assertEquals(gaps.getOpenPenalty(), gop);
        assertEquals(gaps.getExtensionPenalty(), gep);
        assertEquals(gaps.getType(), GapPenalty.Type.AFFINE);
    }

    @Test
    public void testExtensionPenalty() {
        GapPenalty gaps = new SimpleGapPenalty();
        short gep = 14;
        gaps.setExtensionPenalty(gep);
        assertEquals(gaps.getExtensionPenalty(), gep);
    }

    @Test
    public void testOpenPenalty() {
        GapPenalty gaps = new SimpleGapPenalty();
        short gop = 27;
        gaps.setOpenPenalty(gop);
        assertEquals(gaps.getOpenPenalty(), gop);
    }

    @Test
    public void testType() {
        assertEquals(new SimpleGapPenalty((short) 7, (short) 0).getType(), GapPenalty.Type.CONSTANT);
        assertEquals(new SimpleGapPenalty((short) 5, (short) 5).getType(), GapPenalty.Type.LINEAR);
        assertEquals(new SimpleGapPenalty((short) 8, (short) 3).getType(), GapPenalty.Type.AFFINE);
    }

    @Test
    public void testDefaults() {
        short gop = 5;
        short gep = 0;
        SimpleGapPenalty.Defaults.setOpenPenalty(gop);
        SimpleGapPenalty.Defaults.setExtensionPenalty(gep);
        GapPenalty gaps = new SimpleGapPenalty();
        assertEquals(gaps.getOpenPenalty(), gop);
        assertEquals(gaps.getExtensionPenalty(), gep);
        assertEquals(gaps.getType(), GapPenalty.Type.CONSTANT);
    }

}
