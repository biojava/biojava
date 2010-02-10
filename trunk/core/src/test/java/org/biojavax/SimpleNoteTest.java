/*
 * SimpleNoteTest.java
 * JUnit based test
 *
 * Created on 11 November 2005, 22:07
 */

package org.biojavax;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojavax.ontology.ComparableOntology;
import org.biojavax.ontology.ComparableTerm;


/**
 *
 * @author Mark Schreiber
 */
public class SimpleNoteTest extends TestCase {
    ComparableOntology ont;
    SimpleNote note;
    SimpleNote note2;
    SimpleNote note3;
    SimpleNote note4;
    ComparableTerm term;
    
    String value = "test_val";
    int rank = 1;
    
    public SimpleNoteTest(String testName) {
        super(testName);
        ont = RichObjectFactory.getDefaultOntology();
        term = ont.getOrCreateTerm("test_term");
    }

    protected void setUp() throws Exception {
        note = new SimpleNote(term, value, rank);
        note2 = new SimpleNote(term, value, rank);
        note3 = new SimpleNote(term, "another_val", rank);
        note4 = new SimpleNote(term, value, 22);
    }

    protected void tearDown() throws Exception {
        note = note2 = note3 = note4 = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleNoteTest.class);
        
        return suite;
    }

    /**
     * Test of getTerm method, of class org.biojavax.SimpleNote.
     */
    public void testGetTerm() {
        System.out.println("testGetTerm");
        
        assertEquals(term, note.getTerm());
    }

    /**
     * Test of setTerm method, of class org.biojavax.SimpleNote.
     */
    public void testSetTerm() {
        System.out.println("testSetTerm");
        
        try{
            note.setTerm(null);
            fail("Should have got IllegalArgumentException");
        }catch (Exception ex){}
        
        ComparableTerm term2 = ont.getOrCreateTerm("term2");
        try{
            note.setTerm(term2);
            assertEquals(term2, note.getTerm());
        }catch (Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of getValue method, of class org.biojavax.SimpleNote.
     */
    public void testGetValue() {
        System.out.println("testGetValue");
        
        assertEquals(value, note.getValue());
    }

    /**
     * Test of setValue method, of class org.biojavax.SimpleNote.
     */
    public void testSetValue() {
        System.out.println("testSetValue");
        
        String otherVal = "something else";
        try{
            note.setValue(otherVal);
            assertEquals(otherVal, note.getValue());
        }catch (Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of getRank method, of class org.biojavax.SimpleNote.
     */
    public void testGetRank() {
        System.out.println("testGetRank");
        
        assertEquals(rank, note.getRank());
    }

    /**
     * Test of setRank method, of class org.biojavax.SimpleNote.
     */
    public void testSetRank() {
        System.out.println("testSetRank");
        
        try{
            note.setRank(100);
            assertEquals(100, note.getRank());
        }catch (Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleNote.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(note.compareTo(note2) == 0);
        assertTrue(note.compareTo(note3) == 0);
        assertTrue(note.compareTo(note4) < 1);
        assertTrue(note4.compareTo(note) > 1);
    }

    /**
     * Test of equals method, of class org.biojavax.SimpleNote.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(note.equals(note2));
        assertTrue(note.equals(note3)); //notes are equal if they have the same term and rank
        assertFalse(note.equals(note4));
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleNote.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        assertEquals(note.hashCode(), note2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleNote.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = "(#"+rank+") "+term+": "+value;
        assertEquals(expected, note.toString());
    }
    
}
