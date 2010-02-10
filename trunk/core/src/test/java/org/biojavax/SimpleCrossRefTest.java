/*
 * SimpleCrossRefTest.java
 * JUnit based test
 *
 * Created on November 10, 2005, 2:10 PM
 */

package org.biojavax;

import java.util.HashSet;
import java.util.Set;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeVetoException;



/**
 *
 * @author Mark Schreiber
 * @autor George Waldon
 */
public class SimpleCrossRefTest extends TestCase {
    SimpleCrossRef xref;
    String dbname = "test db";
    String accession = "1234TEST";
    int version = 99;
    
    public SimpleCrossRefTest(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
        xref = new SimpleCrossRef(dbname, accession, version);
    }

    protected void tearDown() throws Exception {
        xref = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleCrossRefTest.class);
        
        return suite;
    }

    /**
     * Test of getAnnotation  & getRichAnnotation methods of class org.biojavax.SimpleCrossRef.
     */
    public void testGetAnnotation() {
        System.out.println("testGetAnnotation");
        
        assertNotNull(xref.getAnnotation());
        assertNotNull(xref.getRichAnnotation());

        //should be an editable annotation
        try{
            xref.getAnnotation().setProperty("key", "value");
        }catch(Exception e){
            fail("Was expecting to be able to edit the annotation");
        }
        try{
            xref.getRichAnnotation().setProperty("key2", "value2");
        }catch(Exception e){
            fail("Was expecting to be able to edit the rich annotation");
        }
    }

    /**
     * Test of getNoteSet method, of class org.biojavax.SimpleCrossRef.
     */
    public void testGetNoteSet() {
        System.out.println("testGetNoteSet");
        
        assertNotNull(xref.getNoteSet());
    }

    /**
     * Test of setNoteSet method, of class org.biojavax.SimpleCrossRef.
     */
    public void testSetNoteSet() {
        System.out.println("testSetNoteSet");
        
        Set noteSet = new HashSet();
        try{
            xref.setNoteSet(noteSet);
            assertEquals(noteSet, xref.getNoteSet());
        }catch(ChangeVetoException ex){
            fail("Was expecting to be able to add a new note set, got ChangeVetoException");
        }
    }

    /**
     * Test of getAccession method, of class org.biojavax.SimpleCrossRef.
     */
    public void testGetAccession() {
        System.out.println("testGetAccession");
        
        assertEquals(accession, xref.getAccession());
    }

    /**
     * Test of getDbname method, of class org.biojavax.SimpleCrossRef.
     */
    public void testGetDbname() {
        System.out.println("testGetDbname");
        
        assertEquals(dbname, xref.getDbname());
    }

    /**
     * Test of getVersion method, of class org.biojavax.SimpleCrossRef.
     */
    public void testGetVersion() {
        System.out.println("testGetVersion");
        
        assertEquals(version, xref.getVersion());
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleCrossRef.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        SimpleCrossRef before = new SimpleCrossRef("a", accession, version);
        assertTrue(before.compareTo(xref) < 0);
        assertTrue(xref.compareTo(before) > 0);
        before = new SimpleCrossRef(dbname, "1111TEST", version);
        assertTrue(before.compareTo(xref) < 0);
        assertTrue(xref.compareTo(before) > 0);
        before = new SimpleCrossRef(dbname, accession, 0);
        assertTrue(before.compareTo(xref) < 0);
        assertTrue(xref.compareTo(before) > 0);
        
        SimpleCrossRef equal = new SimpleCrossRef(dbname, accession, version);
        assertTrue(xref.compareTo(equal) == 0);
        assertTrue(equal.compareTo(xref) == 0);
        
        SimpleCrossRef after = new SimpleCrossRef("z", accession, version);
        assertTrue(after.compareTo(xref) > 0);
        assertTrue(xref.compareTo(after) < 0);
        after = new SimpleCrossRef(dbname, "z", version);
        assertTrue(after.compareTo(xref) > 0);
        assertTrue(xref.compareTo(after) < 0);
        after = new SimpleCrossRef(dbname, accession, 9999);
        assertTrue(after.compareTo(xref) > 0);
        assertTrue(xref.compareTo(after) < 0);
        
        
    }

    /**
     * Test of equals method, of class org.biojavax.SimpleCrossRef.
     */
    public void testEquals() {
        System.out.println("testEquals");
        SimpleCrossRef xref2 = new SimpleCrossRef(dbname, accession, version);
        assertTrue(xref.equals(xref2));
        assertTrue(xref2.equals(xref));
        
        xref2 = new SimpleCrossRef("other", accession, version);
        assertFalse(xref.equals(xref2));
        assertFalse(xref2.equals(xref));
        
        xref2 = new SimpleCrossRef(dbname, "jjj", version);
        assertFalse(xref.equals(xref2));
        assertFalse(xref2.equals(xref));
        
        xref2 = new SimpleCrossRef(dbname, accession, 0);
        assertFalse(xref.equals(xref2));
        assertFalse(xref2.equals(xref));
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleCrossRef.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        SimpleCrossRef xref2 = new SimpleCrossRef(dbname, accession, version);
        assertEquals(xref.hashCode(), xref2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleCrossRef.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = dbname+":"+accession+"."+version;
        assertEquals(expected, xref.toString());
    }
    
}
