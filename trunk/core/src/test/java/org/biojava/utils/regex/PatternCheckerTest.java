


package org.biojava.utils.regex;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;

public class PatternCheckerTest
    extends TestCase
{
    String [] patterns = 
        { "tcag",
          "tnag",
          "tc(ac|gt)ag",
          "tc[ac]gt",
          "tc?gt",
          "tc??gt" };

    String [] expected =
        { "tcag", 
          "t[acgt]ag",
          "tc(ac|gt)ag",
          "tc[ac]gt",
          "tc?gt",
          "tc??gt" };

    public void testPatternChecker()
        throws Exception
    {
        for (int i=0; i < patterns.length; i++) {
            String testString = patterns[i];
            PatternChecker checker = new PatternChecker(DNATools.getDNA());
            String result = checker.parse(testString);
//            System.out.println(testString + " " + result + " " + expected[i]);
            assertEquals(
                "failure on test " + i,
                result, 
                expected[i]);
        }
    }
}

