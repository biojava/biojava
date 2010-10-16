package org.biojavax.bio.seq.io;

import junit.framework.TestCase;
import org.biojava.bio.seq.io.ParseException;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichLocation;

/**
 * This class will test parsing of GenBank locations.
 * 
 * @author dsheoran
 */
public class GenbankLocationParserTest extends TestCase {

    /**
     * This method tests if GenbankLocationParser class can successfully parse
     * locations found in GenBank records.
     * @author dsheoran
     */
    public void testParseLocation() throws Exception {
        String[] location = new String[]{
            "467",
            "340..565",
            "<345..500",
            "<1..888",
            "(102.110)",
            "(23.45)..600",
            "(122.133)..(204.221)",
            "123^124",
            "145^177",
            "join(12..78,134..202)",
            "complement(1..23)",
            "complement(join(2691..4571,4918..5163)",
            "join(complement(4918..5163),complement(2691..4571))",
            "complement(34..(122.126))",
            "complement((122.126)..34)",
            "J00194:100..202",
            "(8298.8300)..10206",
            "join((8298.8300)..10206,1..855)"}; // Found in M32882
       
         for (String loc : location) {
            try {
                RichLocation parseLocation =
                        GenbankLocationParser.parseLocation(new SimpleNamespace("gb"), "dummy", loc);
                System.out.println(
                        "Location " + loc + ": " + GenbankLocationParser.writeLocation(parseLocation));
            } catch (ParseException pex) {
                fail("Location " + loc + "was not parsed properly");
            }
        }
    }
}

