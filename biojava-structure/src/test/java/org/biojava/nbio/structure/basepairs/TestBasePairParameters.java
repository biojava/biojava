package org.biojava.nbio.structure.basepairs;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.basepairs.BasePairParameters;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * Created by luke czapla on 7/21/17.
 */
public class TestBasePairParameters {

    @Test
    public void testBasePair() {

        Structure structure;
        try {
            structure = StructureIO.getStructure("1KX5");
        } catch (IOException|StructureException e) {
            e.printStackTrace();
            structure = null;
            assertEquals(1, 2);
        }
        BasePairParameters bp = new BasePairParameters(structure);
        double[][] pairs = bp.getPairingParameters();
        double[][] steps = bp.getStepParameters();
        String sequence = bp.getPairSequence();

        assertEquals(sequence.trim().length(), 147);
        // below all this set of comparator data was from an external program, 3DNA.
        // next three in degrees: buckle, propeller, opening
        assertEquals(pairs[0][0], -3.796, 0.1);
        assertEquals(pairs[0][1], 4.482, 0.1);
        assertEquals(pairs[0][2], -0.730, 0.1);
        // next three in Å: shear, stretch, stagger
        assertEquals(pairs[0][3], -0.324, 0.01);
        assertEquals(pairs[0][4], -0.578, 0.01);
        assertEquals(pairs[0][5], -0.336, 0.01);
        // next three in degrees: tilt, roll, twist
        assertEquals(steps[1][0], 2.354, 0.1);
        assertEquals(steps[1][1], 0.785, 0.1);
        assertEquals(steps[1][2], 32.522, 1.0);
        // next three in Å, shift, slide, rise
        assertEquals(steps[1][3], -0.873, 0.01);
        assertEquals(steps[1][4], -0.607, 0.01);
        assertEquals(steps[1][5], 3.070, 0.01);

        try {
            structure = StructureIO.getStructure("3PHP");
        } catch (IOException|StructureException e) {
            e.printStackTrace();
            structure = null;
            assertEquals(1, 2);
        }
        bp = new TertiaryBasePairParameters(structure, true, false);
        assertEquals(9, bp.getPairingParameters().length);

    }

}

