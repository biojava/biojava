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
 */

package org.biojava.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for the <code>SmallMap</code> class.
 * 
 * @author Len Trigg
 */
public class SmallMapTest extends TestCase {
    
    protected SmallMap mMap = null;

    public SmallMapTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        mMap = new SmallMap();
    }


    public void tearDown() throws Exception {
        mMap = null;
    }
    

    // A test that does some simple map contents manipulation
    public void testPutRemove() throws Exception {
        String annoTag = "anno1";
        String annoVal = "val1";
        String anno2Tag = "anno2";
        String anno2Val = "val2";
        String anno3Tag = "anno3";
        String anno3Val = "val3";

        // Test putting a single value in
        mMap.put(annoTag, annoVal);
        assertTrue(mMap.containsKey(annoTag));
        assertTrue(!mMap.containsKey(anno2Tag));

        // Test removing the only value
        mMap.remove(annoTag);
        assertTrue(!mMap.containsKey(annoTag));
        mMap.put(annoTag, annoVal);

        // Test adding a second value
        mMap.put(anno2Tag, anno2Val);
        assertTrue(mMap.containsKey(anno2Tag));

        // Test removing the first value
        // This triggered a bug in SmallMap
        mMap.remove(annoTag);
        assertTrue(!mMap.containsKey(annoTag));
        assertTrue(mMap.containsKey(anno2Tag));

        // Test removing the second value
        mMap.remove(anno2Tag);
        assertTrue(!mMap.containsKey(anno2Tag));

        assertEquals("Should be empty, but was: " + mMap.toString(), 0, mMap.size());
        mMap.put(annoTag, annoVal);
        mMap.put(anno2Tag, anno2Val);
        mMap.put(anno3Tag, anno3Val);
        assertEquals(3, mMap.size());
        assertTrue(mMap.containsKey(annoTag));
        assertTrue(mMap.containsKey(anno2Tag));
        assertTrue(mMap.containsKey(anno3Tag));
        assertEquals(annoVal, mMap.get(annoTag));
        assertEquals(anno2Val, mMap.get(anno2Tag));
        assertEquals(anno3Val, mMap.get(anno3Tag));

        mMap.remove(anno2Tag);
        assertEquals(2, mMap.size());
        assertTrue(mMap.containsKey(annoTag));
        assertTrue(!mMap.containsKey(anno2Tag));
        assertTrue(mMap.containsKey(anno3Tag));
        assertEquals(annoVal, mMap.get(annoTag));
        assertEquals(anno3Val, mMap.get(anno3Tag));

        mMap.remove(anno3Tag);
        assertEquals(1, mMap.size());
        assertTrue(mMap.containsKey(annoTag));
        assertTrue(!mMap.containsKey(anno2Tag));
        assertTrue(!mMap.containsKey(anno3Tag));
        assertEquals(annoVal, mMap.get(annoTag));

        mMap.remove(annoTag);
        assertEquals(0, mMap.size());
        assertTrue(!mMap.containsKey(annoTag));
        assertTrue(!mMap.containsKey(anno2Tag));
        assertTrue(!mMap.containsKey(anno3Tag));
    }


    public static Test suite() {
        return new TestSuite(SmallMapTest.class);
    }
    
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }    
}
