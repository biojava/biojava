import org.biojava.nbio.structure.secstruc.BasePairParameters;
import org.junit.Test;

/**
 * Created by luke on 7/21/17.
 */
public class TestBasePairParameters {
    @Test
    public void testBasePair() {
        BasePairParameters bp = new BasePairParameters("1KX5");
        double[][] pairs = bp.getPairingParameters();
        double[][] steps = bp.getStepParameters();
        int pos = 0;
        String sequence = bp.getPairSequence();
        System.out.println("buckle propeller opening shear stretch stagger tilt roll twist shift slide rise");
        for (int i = 0; i < pairs.length; i++) {
            while (sequence.charAt(pos) == ' ') pos++;
            System.out.print(sequence.charAt(pos) + ": ");
            for (int j = 0; j < 6; j++)
                System.out.print(pairs[i][j] + " ");
            System.out.print("  ");
            for (int j = 0; j < 6; j++)
                System.out.print(steps[i][j] + " ");
            System.out.println();
            pos++;
        }

    }

}

