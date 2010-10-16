/*
 * SimpleRankedDocRefTest.java
 * JUnit based test
 *
 * Created on 12 November 2005, 15:43
 */

package org.biojavax;

import java.util.Collections;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ChangeListener.ChangeEventRecorder;
import org.biojavax.bio.seq.Position;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.SimplePosition;
import org.biojavax.bio.seq.SimpleRichLocation;

/**
 *
 * @author Mark Schreiber
 * @author gwaldon
 */
public class SimpleRankedDocRefTest extends TestCase {
    DocRef dr;
    SimpleRankedDocRef ref;
    SimpleRankedDocRef ref2;
    int rank = 1;
    Integer start;
    Integer end;
    ChangeEventRecorder cr;
    
    public SimpleRankedDocRefTest(String testName) {
        super(testName);
        start = new Integer(1);
        end = new Integer(25);
        dr = new SimpleDocRef(Collections.singletonList(
                new SimpleDocRefAuthor("Hubert Hubertson", false, false)), "Journal of Voodoo Virology", "Viruses, what are they good for?");
    }

    protected void setUp() throws Exception {
        ref = new SimpleRankedDocRef(dr, start, end, rank);
        cr = new ChangeEventRecorder();
        ref.addChangeListener(cr);
    }

    protected void tearDown() throws Exception {
        ref.removeChangeListener(cr);
        ref = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleRankedDocRefTest.class);
        
        return suite;
    }
    
    /**
     * Test of setRank method, of class org.biojavax.SimpleRankedDocRef.
     */ 
    public void testSetRank() {
        System.out.println("testSetRank");
        try {
            ref.setRank(2);
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(RankedDocRef.RANK, ev.getType());
            //old value should be Integer(1);
            assertEquals(new Integer(1), ev.getPrevious());
            //new value should be Integer(2);
             assertEquals(new Integer(2), ev.getChange());
             ref.setRank(1);
        } catch (ChangeVetoException cve) {
            fail("Unexpected exception: "+ cve);
        }
    }

    /**
     * Test of getRank method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testGetRank() {
        System.out.println("testGetRank");
        
        assertEquals(rank, ref.getRank());
    }

    /**
     * Test of getDocumentReference method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testGetDocumentReference() {
        System.out.println("testGetDocumentReference");
        
        assertEquals(dr, ref.getDocumentReference());
    }

    /**
     * Test of getStart method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testGetStart() {
        System.out.println("testGetStart");
        
        assertEquals(start, ref.getStart());
    }

    /**
     * Test of getEnd method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testGetEnd() {
        System.out.println("testGetEnd");
        
        assertEquals(end, ref.getEnd());
    }
    
    /**
     * Test of setLocation method, of class org.biojavax.SimpleRankedDocRef.
     */ 
    public void testSetLocation() {
        System.out.println("testSetLocation");
        try {
            Position p1 = new SimplePosition(2);
            Position p2 = new SimplePosition(4);
            RichLocation loc = new SimpleRichLocation(p1,p2,0);
            ref.setLocation(loc);
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(RankedDocRef.LOCATION, ev.getType());
            //old value;
            Object o = ev.getPrevious();
            assertTrue(o instanceof RichLocation);
            RichLocation l = (RichLocation) o;
            assertTrue(l.getMin()==1);
            assertTrue(l.getMax()==25);
            //new value;
            o = ev.getChange();
            assertTrue(o instanceof RichLocation);
            l = (RichLocation) o;
            assertTrue(l.getMin()==2);
            assertTrue(l.getMax()==4);
             
            p1 = new SimplePosition(1);
            p2 = new SimplePosition(25);
            loc = new SimpleRichLocation(p1,p2,0);
            ref.setLocation(loc);
            
        } catch (ChangeVetoException cve) {
            fail("Unexpected exception: "+ cve);
        }
    }
    
    /**
     * Test of equals method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(ref.equals(ref));
        assertFalse(ref.equals(new Object()));
        assertFalse(ref.equals(null));
        //Two ranked document references are equal if they have the same rank 
        //and refer to the same document reference.
        ref2 = new SimpleRankedDocRef(dr, start, end, 1); //equal
        assertTrue(ref.equals(ref2));
        assertTrue(ref2.equals(ref));
        
        ref2 = new SimpleRankedDocRef(dr, new Integer(30), new Integer(60), 1); //not equal
        assertFalse(ref.equals(ref2));
        assertFalse(ref2.equals(ref));
        
        ref2 = new SimpleRankedDocRef(dr, start, end, 100); //not equal
        assertFalse(ref.equals(ref2));
        assertFalse(ref2.equals(ref));
        
        ref2 = new SimpleRankedDocRef(new SimpleDocRef(
                Collections.singletonList(new SimpleDocRefAuthor("Rev. Falliwell", false, false)), 
                "Kansas Journal of Creationism", "Un-intelligent design"), start, end, 1); //not equal
        assertFalse(ref.equals(ref2));
        assertFalse(ref2.equals(ref));
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(ref.compareTo(ref) == 0);

        //Two ranked document references are equal if they have the same rank and location 
        //and refer to the same document reference.
        ref2 = new SimpleRankedDocRef(dr, start, end, 1); //equal
        assertTrue(ref.compareTo(ref2) == 0);
        assertTrue(ref2.compareTo(ref) == 0);
        
        ref2 = new SimpleRankedDocRef(dr, new Integer(30), new Integer(60), 1); //not equal
        assertTrue(ref.compareTo(ref2) < 0);
        assertTrue(ref2.compareTo(ref) > 0);
        
        ref2 = new SimpleRankedDocRef(dr, start, end, 100); //not equal
        assertTrue(ref.compareTo(ref2) < 0);
        assertTrue(ref2.compareTo(ref) > 0);
        
        ref2 = new SimpleRankedDocRef(new SimpleDocRef(
                Collections.singletonList(new SimpleDocRefAuthor("Rev. Falliwell", false, false)), 
                "Kansas Journal of Creationism", "Evidence for the giant spaghetti monster"), start, end, 1); //not equal
        assertTrue(ref.compareTo(ref2) == ref.getDocumentReference().compareTo(ref2.getDocumentReference())); //everything else is the same
        assertTrue(ref2.compareTo(ref) == ref2.getDocumentReference().compareTo(ref.getDocumentReference())); //everything else is the same
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        ref2 = new SimpleRankedDocRef(dr, start, end, 1); //equal
        assertTrue(ref.hashCode() == ref2.hashCode());
        
        ref2 = new SimpleRankedDocRef(dr, new Integer(30), new Integer(60), 1); //not equal
        assertTrue(ref.hashCode() != ref2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleRankedDocRef.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = "(#"+rank+") "+dr;
        assertEquals(expected, ref.toString());
    }
}
