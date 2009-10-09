package org.biojavax.bio.seq;

import java.util.ArrayList;
import java.util.Collection;
import junit.framework.TestCase;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;

/**
 *
 * @author George Waldon
 */
public class RichLocationToolsTest extends TestCase {

    public RichLocationToolsTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {

    }

    @Override
    protected void tearDown() throws Exception {

    }

    /** Test RichLocation.Tools.contruct //partial
     *
     */
    public void testconstruct() {
        //case members.size() == 1
        Collection<Location> members = new ArrayList<Location>();
        members.add(new RangeLocation(1,10));
        assertTrue(RichLocation.Tools.construct(members)instanceof RichLocation);
    }
}
