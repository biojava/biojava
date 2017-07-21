import org.biojava.nbio.structure.secstruc.BasePairParameters;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Created by luke on 7/21/17.
 */
public class TestBasePairParameters {

    @Test
    public void testBasePair() {

        BasePairParameters bp = new BasePairParameters("1KX5");
        double[][] pairs = bp.getPairingParameters();
        double[][] steps = bp.getStepParameters();
        String sequence = bp.getPairSequence();

        assertEquals(sequence.trim().length(), 147);
        // below: next 3 in degrees, all this data was from an external program, 3DNA
        assertEquals(pairs[0][0], -3.796, 0.1);
        assertEquals(pairs[0][1], 4.482, 0.1);
        assertEquals(pairs[0][2], -0.730, 0.1);
        // in Å
        assertEquals(pairs[0][3], -0.324, 0.01);
        assertEquals(pairs[0][4], -0.578, 0.01);
        assertEquals(pairs[0][5], -0.336, 0.01);
        // in degrees
        assertEquals(steps[1][0], 2.354, 0.1);
        assertEquals(steps[1][1], 0.785, 0.1);
        assertEquals(steps[1][2], 32.522, 1.0);
        // in Å
        assertEquals(steps[1][3], -0.873, 0.01);
        assertEquals(steps[1][4], -0.607, 0.01);
        assertEquals(steps[1][5], 3.070, 0.01);

    }

}

