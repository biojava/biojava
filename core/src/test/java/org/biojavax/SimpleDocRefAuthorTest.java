/*
 * SimpleDocRefAuthorTest.java
 * JUnit based test
 *
 * Created on November 11, 2005, 4:43 PM
 */

package org.biojavax;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 *
 * @author Mark Schreiber
 */
public class SimpleDocRefAuthorTest extends TestCase {
    SimpleDocRefAuthor auth1;
    SimpleDocRefAuthor auth2;
    SimpleDocRefAuthor ed1;
    SimpleDocRefAuthor ed2;
    SimpleDocRefAuthor con1;
    SimpleDocRefAuthor con2;
    
    String name;
    
    public SimpleDocRefAuthorTest(String testName) {
        super(testName);
        name = "Hubert Hubertson";
    }

    protected void setUp() throws Exception {
        auth1 = new SimpleDocRefAuthor(name, false, false);
        auth2 = new SimpleDocRefAuthor(name);
        
        ed1 = new SimpleDocRefAuthor(name, false, true);
        ed2 = new SimpleDocRefAuthor(name+" (ed.)");
        
        con1 = new SimpleDocRefAuthor(name, true, false);
        con2 = new SimpleDocRefAuthor(name+" (consortium)");
    }

    protected void tearDown() throws Exception {
        auth1 = auth2 = ed1 = ed2 = con1 = con2 = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleDocRefAuthorTest.class);
        
        return suite;
    }

    /**
     * Test of getName method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        assertEquals(name, auth1.getName());
        assertEquals(name, auth2.getName());
        assertEquals(name, ed1.getName());
        assertEquals(name, ed2.getName());
        assertEquals(name, con1.getName());
        assertEquals(name, con2.getName());
    }

    /**
     * Test of getExtendedName method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testGetExtendedName() {
        System.out.println("testGetExtendedName");
        
        assertEquals(name, auth1.getExtendedName());
        assertEquals(ed1.getExtendedName(), ed2.getExtendedName());
        assertEquals(name+" (ed.)", ed1.getExtendedName());
        assertEquals(con1.getExtendedName(), con2.getExtendedName());
        assertEquals(name+" (consortium)", con1.getExtendedName());
    }

    /**
     * Test of isEditor method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testIsEditor() {
        assertTrue(ed1.isEditor());
        assertTrue(ed2.isEditor());
        assertFalse(auth1.isEditor());
        assertFalse(auth2.isEditor());
        assertFalse(con1.isEditor());
        assertFalse(con2.isEditor());
    }

    /**
     * Test of isConsortium method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testIsConsortium() {
        System.out.println("testIsConsortium");
        
        assertTrue(con1.isConsortium());
        assertTrue(con2.isConsortium());
        assertFalse(auth1.isConsortium());
        assertFalse(auth2.isConsortium());
        assertFalse(ed1.isConsortium());
        assertFalse(ed2.isConsortium());
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(con1.compareTo(con2) == 0);
        assertTrue(con2.compareTo(con1) == 0);
        assertTrue(ed1.compareTo(ed2) == 0);
        assertTrue(ed2.compareTo(ed1) == 0);
        assertTrue(auth1.compareTo(auth2) == 0);
        assertTrue(auth2.compareTo(auth1) == 0);

    }
    

    /**
     * Test of equals method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(ed1.equals(ed2));
        assertTrue(auth1.equals(auth2));
        assertTrue(con1.equals(con2));
        assertTrue(con2.equals(con1));
        assertTrue(ed2.equals(ed1));
        assertTrue(auth2.equals(auth1));
        assertFalse(ed1.equals(con1));
        assertFalse(ed1.equals(auth1));
        assertFalse(con1.equals(ed1));
        assertFalse(auth1.equals(ed1));
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        assertEquals(con1.hashCode(), con2.hashCode());
        assertEquals(auth1.hashCode(), auth2.hashCode());
        assertEquals(ed1.hashCode(), ed2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleDocRefAuthor.
     */
    public void testToString() {
        System.out.println("testToString");
        
        assertEquals(auth2.toString(), auth2.getExtendedName());
    }
    
}
