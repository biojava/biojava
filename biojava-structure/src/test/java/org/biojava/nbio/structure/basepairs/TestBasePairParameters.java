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
 *
 */
package org.biojava.nbio.structure.basepairs;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.basepairs.BasePairParameters;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * This class tests the implementations of the search for base pairs for different RCSB structures
 * and the tests uses 3DNA as a comparator program. (other programs such as CURVES and NEWHELIX work similarly but
 * this implementation is closest to that of 3DNA)
 * @author Luke Czapla
 * @since 5.0.0
 *
 */
public class TestBasePairParameters {

    @Test
    public void testBasePair() throws IOException, StructureException {

        Structure structure = StructureIO.getStructure("1KX5");

        BasePairParameters bp = new BasePairParameters(structure);
        bp.analyze();
        //String sequence = bp.getPairSequence();

        assertEquals(147, bp.getLength());
        // below all this set of comparator data was from an external program, 3DNA.
        // next three in degrees: buckle, propeller, opening
        assertEquals(bp.getBuckle(0), -3.796, 0.1);
        assertEquals(bp.getPropeller(0), 4.482, 0.1);
        assertEquals(bp.getOpening(0), -0.730, 0.1);
        // next three in Å: shear, stretch, stagger
        assertEquals(bp.getShear(0), -0.324, 0.02);
        assertEquals(bp.getStretch(0), -0.578, 0.02);
        assertEquals(bp.getStagger(0), -0.336, 0.02);
        // next three in degrees: tilt, roll, twist
        assertEquals(bp.getTilt(1), 2.354, 0.1);
        assertEquals(bp.getRoll(1), 0.785, 0.1);
        assertEquals(bp.getTwist(1), 32.522, 0.5);
        // next three in Å, shift, slide, rise
        assertEquals(bp.getShift(1), -0.873, 0.02);
        assertEquals(bp.getSlide(1), -0.607, 0.02);
        assertEquals(bp.getRise(1), 3.070, 0.02);


        structure = StructureIO.getStructure("3PHP");
        bp = new TertiaryBasePairParameters(structure, true, false).analyze();
        assertEquals(9, bp.getLength());

        double[][] pairs = bp.getPairingParameters();
        double[][] steps = bp.getStepParameters();

        // test against values given by 3DNA, just using the raw arrays
        assertEquals(pairs[4][0], 0.060, 0.1);
        assertEquals(pairs[4][1], -9.323, 0.1);
        assertEquals(pairs[4][2], -5.109, 0.1);
        // next three in Å: shear, stretch, stagger
        assertEquals(pairs[4][3], 0.126, 0.02);
        assertEquals(pairs[4][4], -0.177, 0.02);
        assertEquals(pairs[4][5], 0.273, 0.02);
        // next three in degrees: tilt, roll, twist
        assertEquals(steps[4][0], -1.456, 0.1);
        assertEquals(steps[4][1], 6.583, 0.1);
        assertEquals(steps[4][2], 33.234, 0.5);
        // next three in Å, shift, slide, rise
        assertEquals(steps[4][3], -0.735, 0.02);
        assertEquals(steps[4][4], -0.978, 0.02);
        assertEquals(steps[4][5], 3.491, 0.02);


        structure = StructureIO.getStructure("1P71");

        bp = new MismatchedBasePairParameters(structure, false, false, false).analyze();
        assertEquals(17, bp.getLength());

        pairs = bp.getPairingParameters();
        steps = bp.getStepParameters();

        // this was tested against 3DNA as well.
        assertEquals(pairs[16][0], -11.822, 0.1);
        assertEquals(pairs[16][1], -11.405, 0.1);
        assertEquals(pairs[16][2], -9.669, 0.1);
        // next three in Å: shear, stretch, stagger
        assertEquals(pairs[16][3], 0.855, 0.02);
        assertEquals(pairs[16][4], -0.276, 0.02);
        assertEquals(pairs[16][5], -0.604, 0.02);
        // next three in degrees: tilt, roll, twist
        assertEquals(steps[16][0], 1.516, 0.1);
        assertEquals(steps[16][1], 9.291, 0.1);
        assertEquals(steps[16][2], 42.052, 1.0);
        // next three in Å, shift, slide, rise
        assertEquals(steps[16][3], -0.627, 0.02);
        assertEquals(steps[16][4], -0.858, 0.02);
        assertEquals(steps[16][5], 4.697, 0.02);

    }

}

