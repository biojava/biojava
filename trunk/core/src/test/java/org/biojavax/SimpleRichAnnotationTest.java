/*
 * SimpleRichAnnotationTest.java
 * JUnit based test
 *
 * Created on November 28, 2005, 4:01 PM
 */

package org.biojavax;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojavax.ontology.ComparableOntology;
import org.biojavax.ontology.ComparableTerm;

/**
 *
 * @author Mark Schreiber
 */
public class SimpleRichAnnotationTest extends TestCase {
    private SimpleRichAnnotation anno1;
    private Note note;
    
    public SimpleRichAnnotationTest(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
        anno1 = new SimpleRichAnnotation();
        anno1.setProperty("foo", "bar"); //should convert to a note
        note = (Note)anno1.getNoteSet().iterator().next();
    }

    protected void tearDown() throws Exception {
        anno1 = null;
        note = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleRichAnnotationTest.class);
        
        return suite;
    }

    /**
     * Test of clear method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testClear() {
        System.out.println("testClear");
        //should have one note after setUp()
        assertTrue(anno1.getNoteSet().size() == 1);
        try{
            anno1.clear();
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        //should have no notes now
        assertTrue(anno1.getNoteSet().size() == 0);
    }

    /**
     * Test of asMap method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testAsMap() {
        System.out.println("testAsMap");
        Map map = anno1.asMap();
        assertTrue(map.containsKey(note.getTerm()));
        assertTrue(map.size() == 1);
        assertEquals(map.get(note.getTerm()), "bar");
    }

    /**
     * Test of addNote method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testAddNote() {
        System.out.println("testAddNote");
        
        try{
            anno1.addNote(note);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        //should do nothing already added
        assertTrue(anno1.getNoteSet().size() == 1);
        
        ComparableOntology ont = RichObjectFactory.getDefaultOntology();
        ComparableTerm term = ont.getOrCreateTerm("foo2"); 
        Note note2 = new SimpleNote(term, "bar", 0);
        try{
            anno1.addNote(note2);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        assertTrue(anno1.contains(note));
        assertTrue(anno1.contains(note2));
        assertEquals(2, anno1.getNoteSet().size());
        assertTrue(anno1.keys().size() == 2);
        assertTrue(anno1.containsProperty("foo"));
        assertTrue(anno1.containsProperty(note.getTerm()));
        assertTrue(anno1.containsProperty(note2.getTerm()));
        assertTrue(anno1.containsProperty("foo2"));
    }

    /**
     * Test of contains method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testContains() {
        System.out.println("testContains");
        
        assertTrue(anno1.contains(note));
    }

    /**
     * Test of containsProperty method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testContainsProperty() {
        System.out.println("testContainsProperty");
        
        assertTrue(anno1.containsProperty("foo"));
        ComparableTerm term = RichObjectFactory.getDefaultOntology().getOrCreateTerm("foo");
        //this should also work.
        assertTrue(anno1.containsProperty(term));
    }

    /**
     * Test of getNote method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testGetNote() {
        System.out.println("testGetNote");
        
        assertNotNull(anno1.getNote(note));
        assertEquals(note, anno1.getNote(note));
    }

    /**
     * Test of getProperty method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testGetProperty() {
        System.out.println("testGetProperty");
        
        assertNotNull(anno1.getProperty("foo"));
        assertTrue(anno1.getProperty("foo") instanceof String);
        assertEquals("bar", anno1.getProperty("foo"));
        
        ComparableTerm term = 
                RichObjectFactory.getDefaultOntology().getOrCreateTerm("foo");
        //this should also work.
        assertNotNull(anno1.getProperty(term));
        assertTrue(anno1.getProperty(term) instanceof String);
        assertEquals("bar", anno1.getProperty(term));
    }

    /**
     * Test of getPropertys method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testGetPropertys() {
        System.out.println("testGetPropertys");
        
        //add a note with the same term but different rank
        Note note2 = new SimpleNote(note.getTerm(), "bar1", 1);
        
        try{
            anno1.addNote(note2);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        //should be two notes whichever way you do it
        assertEquals(2, anno1.getProperties(note2.getTerm()).length);
        assertEquals(2, anno1.getProperties(note.getTerm()).length);
        assertEquals(2, anno1.getProperties("foo").length);
        
        Note[] notes = anno1.getProperties(note.getTerm());
        assertEquals(note.getValue(), notes[0].getValue());
        assertEquals(note2.getValue(), notes[1].getValue());
        
        //should be an empty array not null;
        assertNotNull(anno1.getProperties("not_here"));
        assertEquals(0, anno1.getProperties("not_here").length);
    }
    
    /**
     * Test of keys method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testKeys() {
        System.out.println("testKeys");
        
        assertNotNull(anno1.keys());
        Set keys = anno1.keys();
        
        assertTrue(keys.contains(note.getTerm()));
        assertTrue(keys.size() == 1);
    }

    /**
     * Test of removeNote method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testRemoveNote() {
        System.out.println("testRemoveNote");
        
        try{
            anno1.removeNote(note);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        assertTrue(anno1.getNoteSet().size() == 0);
        assertTrue(anno1.keys().size() == 0);
        try{
            anno1.getNote(note);
            fail("Was expecting NoSuchElementException");
        }catch(NoSuchElementException ex){}
        try{
            anno1.getProperty(note.getTerm());
            fail("Was expecting NoSuchElementException");
        }catch(NoSuchElementException ex){}
    }
    
    public void testRemoveProperty() {
        System.out.println("testRemoveProperty");
        
        try{
            anno1.removeProperty("foo");
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        assertTrue(anno1.getNoteSet().size() == 0);
        assertTrue(anno1.keys().size() == 0);
        try{
            anno1.getNote(note);
            fail("Was expecting NoSuchElementException");
        }catch(NoSuchElementException ex){}
        try{
            anno1.getProperty("foo");
            fail("Was expecting NoSuchElementException");
        }catch(NoSuchElementException ex){}
    }
    
    /**
     * Test of removeProperty method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testRemoveProperty2() {
        System.out.println("testRemoveProperty2");
        
        try{
            anno1.removeProperty(note.getTerm());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        assertTrue(anno1.getNoteSet().size() == 0);
        assertTrue(anno1.keys().size() == 0);
        try{
            anno1.getNote(note);
            fail("Was expecting NoSuchElementException");
        }catch(NoSuchElementException ex){}
        try{
            anno1.getProperty(note.getTerm());
            fail("Was expecting NoSuchElementException");
        }catch(NoSuchElementException ex){}
    }    

    /**
     * Test of setProperty method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testSetProperty() {
        System.out.println("testSetProperty");
        
        try{
            anno1.setProperty("foo", "bar");
            //shouldn't add twice
            assertTrue(anno1.getNoteSet().size() == 1);
            assertTrue(anno1.keys().size() == 1);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            anno1.setProperty(null, null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        try{
            anno1.setProperty("foo2", "bar2");
            assertTrue(anno1.getNoteSet().size() == 2);
            assertTrue(anno1.keys().size() == 2);
            assertTrue(anno1.containsProperty("foo2"));
            assertNotNull(anno1.getProperty("foo2"));
            assertEquals("bar2", anno1.getProperty("foo2"));
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of getNoteSet method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testGetNoteSet() {
        System.out.println("testGetNoteSet");
        
        // already thoroughly tested but here we go..
        assertNotNull(anno1.getNoteSet());
        assertTrue(anno1.getNoteSet().size() == 1);
        //should be writable
        try{
            anno1.getNoteSet().remove(note);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        //should not be null at instantiation
        anno1 = new SimpleRichAnnotation();
        assertNotNull(anno1.getNoteSet());
        assertTrue(anno1.getNoteSet().size() == 0);
    }

    /**
     * Test of setNoteSet method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testSetNoteSet() {
        System.out.println("testSetNoteSet");
        
        // this method is really required by Hibernate. You should be allowed to
        // make this null.
        try{
            anno1.setNoteSet(null);
            assertNull(anno1.getNoteSet());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        try{ 
            Set s = new HashSet();
            anno1.setNoteSet(s);
            assertEquals(s, anno1.getNoteSet());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleRichAnnotation.
     */
    public void testToString() {
        System.out.println("testToString");
        
        StringBuffer sb = new StringBuffer();
        for (Iterator i = this.anno1.getNoteSet().iterator(); i.hasNext(); ) {
            sb.append("[");
            sb.append(i.next());
            sb.append("]");
            if (i.hasNext()) sb.append(",");
        }
        String expected = sb.toString();
        assertEquals(expected, anno1.toString());
    }
    
}
