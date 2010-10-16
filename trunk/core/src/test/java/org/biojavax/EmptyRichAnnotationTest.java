/*
 * EmptyRichAnnotationTest.java
 * JUnit based test
 *
 * Created on November 11, 2005, 4:05 PM
 */

package org.biojavax;

import java.util.Collections;
import java.util.HashSet;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeVetoException;


/**
 *
 * @author Mark Schreiber
 */
public class EmptyRichAnnotationTest extends TestCase {
    EmptyRichAnnotation ann;
    
    public EmptyRichAnnotationTest(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
        ann = new EmptyRichAnnotation();
    }

    protected void tearDown() throws Exception {
        ann = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(EmptyRichAnnotationTest.class);
        
        return suite;
    }

    /**
     * Test of getProperty method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testGetProperty() {
        System.out.println("testGetProperty");
        
        assertNull(ann.getProperty(new Object()));
    }

    /**
     * Test of getNote method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testGetNote() {
        System.out.println("testGetNote");
        
        assertNull(ann.getNote(null));
    }

    /**
     * Test of setProperty method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testSetProperty() {
        System.out.println("testSetProperty");
        
        try{
            ann.setProperty("key", "value");
            fail("should throw ChangeVetoException");
        }catch(ChangeVetoException ex){}
    }

    /**
     * Test of setNoteSet method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testSetNoteSet() {
        System.out.println("testSetNoteSet");
        
        try{
            ann.setNoteSet(new HashSet());
            fail("should throw ChangeVetoException");
        }catch(ChangeVetoException ex){}
    }

    /**
     * Test of addNote method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testAddNote() {
        System.out.println("testAddNote");
        
        try{
            ann.addNote(null);
            fail("should throw ChangeVetoException");
        }catch(ChangeVetoException ex){}
    }

    /**
     * Test of clear method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testClear() {
        System.out.println("testClear");
        
        //does nothing but shouldn't throw exception.
        try{
            ann.clear();
        }catch (ChangeVetoException ex){
            fail("wasn't expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of removeProperty method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testRemoveProperty() {
        System.out.println("testRemoveProperty");
        
        try{
            ann.removeProperty(new Object());
            fail("should throw ChangeVetoException");
        }catch(ChangeVetoException ex){}
    }

    /**
     * Test of removeNote method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testRemoveNote() {
        System.out.println("testRemoveNote");
        
        try{
            ann.removeNote(null);
            fail("should throw ChangeVetoException");
        }catch(ChangeVetoException ex){}
    }

    /**
     * Test of containsProperty method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testContainsProperty() {
        System.out.println("testContainsProperty");
        
        assertFalse(ann.containsProperty(new Object()));
    }

    /**
     * Test of contains method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testContains() {
        System.out.println("testContains");
        
        assertFalse(ann.contains(null));
    }

    /**
     * Test of keys method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testKeys() {
        System.out.println("testKeys");
        
        assertEquals(Collections.EMPTY_SET, ann.keys());
    }

    /**
     * Test of getNoteSet method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testGetNoteSet() {
        System.out.println("testGetNoteSet");
        
        assertEquals(Collections.EMPTY_SET, ann.getNoteSet());
    }

    /**
     * Test of asMap method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testAsMap() {
        System.out.println("testAsMap");
        
        assertNotNull(ann.asMap());
    }


    /**
     * Test of equals method, of class org.biojavax.EmptyRichAnnotation.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(ann.equals(new EmptyRichAnnotation()));
        assertFalse(ann.equals(new SimpleRichAnnotation()));
    }
    
}
