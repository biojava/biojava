

package org.biojava.utils.regex;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.SymbolList;

public class RegexTest
    extends TestCase
{
    private String symbols = "atatcgatttagatggatcccgcgta";
    private String patternString = "ggatcc";
    private SymbolList sl;

    protected void setUp()
        throws Exception
    {
        sl = DNATools.createDNA(symbols);
    }

    public void testRegex()
        throws Exception
    {
        // create a factory suitable for DNA
        PatternFactory fact = PatternFactory.makeFactory(DNATools.getDNA());

        assertNotNull("Could not create PatternFactory object.", fact);

        // make a pattern
        Pattern pattern = fact.compile(patternString);

        assertNotNull("Could not create Pattern object.", pattern);

        // make a matcher
        Matcher matcher = pattern.matcher(sl);
        assertNotNull("Could not create Matcher object.", matcher);

        // test matcher
        if (matcher.find()) {
            if ((matcher.start() != 15) && (matcher.end() != 21))
                fail("start and end yielded wrong result.");
        }
        else {
            fail("failed to find target at all!");
        }

        SymbolList subList = matcher.group();
        assertNotNull("failed to extract match group.", matcher);

        assertEquals("failed to find match group correctly.", subList.seqString(), patternString);
    }
}

