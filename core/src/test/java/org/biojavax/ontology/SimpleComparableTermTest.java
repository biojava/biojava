/*
 * SimpleComparableTermTest.java
 * JUnit based test
 *
 * Created on 3 December 2005, 20:55
 */

package org.biojavax.ontology;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeListener.ChangeEventRecorder;
import org.biojavax.RankedCrossRef;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleCrossRef;
import org.biojavax.SimpleRankedCrossRef;

/**
 *
 * @author Mark Schreiber
 */
public class SimpleComparableTermTest extends TestCase {
    SimpleComparableTerm t1;
    SimpleComparableTerm t2;
    ComparableOntology ont;
    String name;
    Object[] synonyms;
    ChangeListener.ChangeEventRecorder cr;
    
    public SimpleComparableTermTest(String testName) {
        super(testName);
        ont = new SimpleComparableOntology("test_onto");
        name = "foo";
        synonyms = new String[]{"foo2", "foo3"};
    }

    protected void setUp() throws Exception {
        cr = new ChangeEventRecorder();
        t1 = new SimpleComparableTerm(ont, name, null);
        t1.addChangeListener(cr);
        t2 = new SimpleComparableTerm(ont, name, synonyms);
    }

    protected void tearDown() throws Exception {
        t1.removeChangeListener(cr);
        cr = null;
        t1 = null;
        t2 = null;
        RichObjectFactory.clearLRUCache(SimpleComparableTerm.class);
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleComparableTermTest.class);
        
        return suite;
    }

    /**
     * Test of hashCode method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        //cache is not considered
        assertEquals(t1.hashCode(), t2.hashCode());
    }

    /**
     * Test of equals method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(t1.equals(t1));
        assertTrue(t1.equals(t2));
        assertTrue(t2.equals(t1));
        
        //not equal
        t2 = new SimpleComparableTerm(ont, "bar", null);
        assertTrue(! t1.equals(t2));
        assertTrue(! t2.equals(t1));
        
        //not equal
        t2 = new SimpleComparableTerm(new SimpleComparableOntology("another_ont"), name, null);
        assertTrue(! t1.equals(t2));
        assertTrue(! t2.equals(t1));
    }

    /**
     * Test of compareTo method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(t1.compareTo(t1) == 0);
        assertTrue(t1.compareTo(t2) == 0);
        assertTrue(t2.compareTo(t1) == 0);
        
        //not equal
        t2 = new SimpleComparableTerm(ont, "bar", null);
        assertTrue( t1.compareTo(t2) > 0);
        assertTrue( t2.compareTo(t1) < 0);
        
        //not equal
        t2 = new SimpleComparableTerm(new SimpleComparableOntology("another_ont"), name, null);
        assertTrue( t1.compareTo(t2) > 0);
        assertTrue( t2.compareTo(t1) < 0);
    }

    /**
     * Test of addSynonym method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testAddSynonym() {
        System.out.println("testAddSynonym");
        String syn = "another_name";
        
        //synonyms don't generate change events
        t1.addSynonym(syn);
        assertTrue(Arrays.binarySearch(t1.getSynonyms(), syn) != -1);
        
        t2.addSynonym(syn);
        //lenght should be 3
        assertEquals(3, t2.getSynonyms().length);
    }

    /**
     * Test of removeSynonym method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testRemoveSynonym() {
        System.out.println("testRemoveSynonym");
        
        String synonym = "foo2";
        t2.removeSynonym(synonym);
        //length should be 1
        assertEquals(1, t2.getSynonyms().length);
        assertTrue(Arrays.binarySearch(t1.getSynonyms(), synonym) == -1);
        
    }

    /**
     * Test of getSynonyms method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetSynonyms() {
        System.out.println("testGetSynonyms");
        // mostly tested already
        // should not be null
        assertNotNull(t1.getSynonyms());
    }

    /**
     * Test of getRankedCrossRefs method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetRankedCrossRefs() {
        System.out.println("testGetRankedCrossRefs");
        
        // cannot be null even if empty
        assertNotNull(t1.getRankedCrossRefs());
        // must be writable for hibernate
        try{
            t1.getRankedCrossRefs().add(new SimpleRankedCrossRef(new SimpleCrossRef("db", "acc", 1),0));
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of setRankedCrossRefs method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testSetRankedCrossRefs() {
        System.out.println("testSetRankedCrossRefs");
        
        try{
            Set s = new TreeSet();
            s.add(new SimpleRankedCrossRef(new SimpleCrossRef("db", "acc", 1),0));
            t1.setRankedCrossRefs(s);
            assertEquals(s, t1.getRankedCrossRefs());
            
            //used by hibernate does not generate events
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
    }

    /**
     * Test of addRankedCrossRef method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testAddRankedCrossRef() {
        System.out.println("testAddRankedCrossRef");
        
        RankedCrossRef xref = new SimpleRankedCrossRef(new SimpleCrossRef("db", "acc", 1),0);
        try{
            t1.addRankedCrossRef(xref);
            assertTrue(t1.getRankedCrossRefs().contains(xref));
            
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(ComparableTerm.RANKEDCROSSREF, ev.getType());
            //old value should be null
            assertNull(ev.getPrevious());
            //new should be xref
            assertEquals(xref, ev.getChange());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of removeRankedCrossRef method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testRemoveRankedCrossRef() {
        System.out.println("testRemoveRankedCrossRef");
        
        RankedCrossRef xref = new SimpleRankedCrossRef(new SimpleCrossRef("db", "acc", 1),0);
        try{
            t1.addRankedCrossRef(xref); //add it
            assertTrue(t1.getRankedCrossRefs().contains(xref));
            t1.removeRankedCrossRef(xref); //now remove it
            
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(ComparableTerm.RANKEDCROSSREF, ev.getType());
            //old value should be xref
            assertEquals(xref, ev.getPrevious());
            //new should be null
            assertNull(ev.getChange());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of getName method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        assertEquals(name, t1.getName());
    }

    /**
     * Test of getDescription method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetDescription() {
        System.out.println("testGetDescription");
        
        //should be null if not set
        assertNull(t1.getDescription());
    }

    /**
     * Test of setDescription method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testSetDescription() {
        System.out.println("testSetDescription");
        
        //can be null
        try{t1.setDescription(null);}catch(Exception ex){ fail("Not expecting "+ex.getClass().getName());}
        
        String desc = "my_desc";
        try{
            t1.setDescription(desc);
            assertEquals(desc, t1.getDescription());
            
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(ComparableTerm.DESCRIPTION, ev.getType());
            //old value should be null
            assertNull(ev.getPrevious());
            //new should be desc
            assertEquals(desc, ev.getChange());
        }catch(Exception ex){fail("Not expecting "+ex.getClass().getName());}
    }

    /**
     * Test of getOntology method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetOntology() {
        System.out.println("testGetOntology");
        
        assertEquals(ont, t1.getOntology());
    }

    /**
     * Test of toString method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = ont+":"+name;
        assertEquals(expected, t1.toString());
        try{
            t1.setObsolete(new Boolean(true));
            expected = ont+":"+name+" [obsolete]";
            assertEquals(expected, t1.toString());
        }catch(Exception ex){fail("Not expecting "+ex.getClass().getName());}
    }

    /**
     * Test of getAnnotation method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetAnnotation() {
        System.out.println("testGetAnnotation");
        
        // always Empty
        assertEquals(RichAnnotation.EMPTY_ANNOTATION, t1.getAnnotation());
    }

    /**
     * Test of getIdentifier method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetIdentifier() {
        System.out.println("testGetIdentifier");
        
        //null at start
        assertNull(t1.getIdentifier());
    }

    /**
     * Test of setIdentifier method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testSetIdentifier() {
        System.out.println("testSetIdentifier");
        
        //can be null
        try{t1.setDescription(null);}catch(Exception ex){ fail("Not expecting "+ex.getClass().getName());}
        
        String ident = "my_ident";
        try{
            t1.setIdentifier(ident);
            assertEquals(ident, t1.getIdentifier());
            
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(ComparableTerm.IDENTIFIER, ev.getType());
            //old value should be null
            assertNull(ev.getPrevious());
            //new should be desc
            assertEquals(ident, ev.getChange());
        }catch(Exception ex){fail("Not expecting "+ex.getClass().getName());}
    }

    /**
     * Test of getObsolete method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testGetObsolete() {
        System.out.println("testGetObsolete");
        
        // false at start
        assertFalse(t1.getObsolete().booleanValue());
    }

    /**
     * Test of setObsolete method, of class org.biojavax.ontology.SimpleComparableTerm.
     */
    public void testSetObsolete() {
        System.out.println("testSetObsolete");
        
        Boolean obs = Boolean.TRUE;
        try{
            t1.setObsolete(obs);
            assertEquals(obs, t1.getObsolete());
            
            //should generate an event
            ChangeEvent ev = cr.getEvent();
            assertNotNull(ev);
            //of the correct type
            assertEquals(ComparableTerm.OBSOLETE, ev.getType());
            //old value should be null
            assertEquals(Boolean.FALSE, ev.getPrevious());
            //new should be true
            assertEquals(obs, ev.getChange());
        }catch(Exception ex){fail("Not expecting "+ex.getClass().getName());}
    }
    
}
