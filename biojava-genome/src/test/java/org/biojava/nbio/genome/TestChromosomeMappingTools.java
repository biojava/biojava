package org.biojava.nbio.genome;

import org.biojava.nbio.genome.util.ChromosomeMappingTools;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Created by Yana Valasatava on 8/14/17.
 */
public class TestChromosomeMappingTools {

    @Test
    public void testGetCDSLengthForward() {

        List<Integer> exonStarts = new ArrayList<>(Arrays.asList(10, 30, 50, 70));
        List<Integer> exonEnds = new ArrayList<>(Arrays.asList(20, 40, 60, 80));
        int cdsStart = 35;
        int cdsEnd = 75;

        int cdsDesired = 23 - 3;
        ChromosomeMappingTools.setCoordinateSystem(0);
        int cdsTest = ChromosomeMappingTools.getCDSLengthForward(exonStarts, exonEnds, cdsStart, cdsEnd);

        assertEquals(cdsDesired, cdsTest);
    }

    @Test
    public void testGetCDSLengthReverseAsc() {

        List<Integer> exonStarts = new ArrayList<>(Arrays.asList(10, 50, 70));
        List<Integer> exonEnds = new ArrayList<>(Arrays.asList(20, 60, 80));
        int cdsStart = 55;
        int cdsEnd = 75;

        int cdsDesired = 12 - 3;
        ChromosomeMappingTools.setCoordinateSystem(0);
        int cdsTest = ChromosomeMappingTools.getCDSLengthReverse(exonStarts, exonEnds, cdsStart, cdsEnd);

        assertEquals(cdsDesired, cdsTest);
    }

    @Test
    public void testGetCDSLengthReverseDesc() {

        List<Integer> exonStarts = new ArrayList<>(Arrays.asList(70, 50, 10));
        List<Integer> exonEnds = new ArrayList<>(Arrays.asList(80, 60, 20));
        int cdsStart = 75;
        int cdsEnd = 50;

        int cdsDesired = 17 - 3;
        ChromosomeMappingTools.setCoordinateSystem(0);
        int cdsTest = ChromosomeMappingTools.getCDSLengthReverse(exonStarts, exonEnds, cdsStart, cdsEnd);

        assertEquals(cdsDesired, cdsTest);
    }
}
