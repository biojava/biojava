/*
 * SimpleRankedCrossRefTest.java
 * JUnit based test
 *
 * Created on 12 November 2005, 15:11
 */

package org.biojavax;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ChangeListener.ChangeEventRecorder;


/**
 *
 * @author Mark Schreiber
 * @author gwaldon
 */
public class SimpleRankedCrossRefTest extends TestCase {
    SimpleCrossRef ref;
    RankedCrossRef xref;
    RankedCrossRef xref2;
    int rank = 1;
    ChangeEventRecorder cr;
    
    public SimpleRankedCrossRefTest(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
        ref = new SimpleCrossRef("dbname", "AC123456", 1);
        xref = new SimpleRankedCrossRef(ref, rank);
        cr = new ChangeEventRecorder();
        xref.addChangeListener(cr);
    }

    protected void tearDown() throws Exception {
        ref = null;
        xref.removeChangeListener(cr);
        xref = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleRankedCrossRefTest.class);
        
        return suite;
    }

    /**
     * Test of getCrossRef method, of class org.biojavax.SimpleRankedCrossRef.
     */
    public void testGetCrossRef() {
        System.out.println("testGetCrossRef");
        
        assertEquals(ref, xref.getCrossRef());
    }

    /**
     * Test of setRank method, of class org.biojavax.SimpleRankedCrossRef.
     */ 
    public void testSetRank() {
        System.out.println("testSetRank");
        try {
            xref.setRank(2);
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(RankedCrossRef.RANK, ev.getType());
            //old value should be Integer(1);
            assertEquals(new Integer(1), ev.getPrevious());
            //new value should be Integer(2);
             assertEquals(new Integer(2), ev.getChange());
             xref.setRank(1);
        } catch (ChangeVetoException cve) {
            fail("Unexpected exception: "+ cve);
        }
    }
    
    /**
     * Test of getRank method, of class org.biojavax.SimpleRankedCrossRef.
     */
    public void testGetRank() {
        System.out.println("testGetRank");
        
        assertEquals(rank, xref.getRank());
    }

    /**
     * Test of equals method, of class org.biojavax.SimpleRankedCrossRef.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(xref.equals(xref));
        assertFalse(xref.equals(ref));
        xref2 = new SimpleRankedCrossRef(ref, rank);
        assertTrue(xref.equals(xref2));
        assertTrue(xref2.equals(xref));
        xref2 = new SimpleRankedCrossRef(ref, 100);
        assertFalse(xref.equals(xref2));
        assertFalse(xref2.equals(xref));
        
        //although equal, not ==
        SimpleCrossRef ref2 = new SimpleCrossRef("dbname", "AC123456", 1);
        xref2 = new SimpleRankedCrossRef(ref2, rank);
        assertTrue(xref.equals(xref2));
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleRankedCrossRef.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(xref.compareTo(xref) == 0);
        
        xref2 = new SimpleRankedCrossRef(ref, rank);
        assertTrue(xref.compareTo(xref2) == 0);
        assertTrue(xref2.compareTo(xref) == 0);
        xref2 = new SimpleRankedCrossRef(ref, 100);
        assertTrue(xref.compareTo(xref2) < 0);
        assertTrue(xref2.compareTo(xref) > 0);
        
        SimpleCrossRef ref2 = new SimpleCrossRef("dbname", "AC123459", 1); //after
        xref2 = new SimpleRankedCrossRef(ref2, rank);
        assertTrue(xref.compareTo(xref2) < 0);
        assertTrue(xref2.compareTo(xref) > 0);
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleRankedCrossRef.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        xref2 = new SimpleRankedCrossRef(ref, rank);
        assertTrue(xref.hashCode() == xref2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleRankedCrossRef.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = "(#"+rank+") "+ref;
        assertEquals(expected, xref.toString());
    }
    
}
