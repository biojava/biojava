/*
 * SimpleComparableTripleTest.java
 * JUnit based test
 *
 * Created on 4 December 2005, 14:57
 */

package org.biojavax.ontology;
import java.util.Set;
import java.util.TreeSet;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.ontology.Triple;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener.ChangeEventRecorder;
import org.biojavax.RichAnnotation;


/**
 *
 * @author Mark Schreiber
 */
public class SimpleComparableTripleTest extends TestCase {
    
    ChangeEventRecorder cr;
    SimpleComparableTriple trip;
    ComparableTerm subj;
    ComparableTerm obj;
    ComparableTerm pred;
    ComparableOntology ont;
    ComparableTerm subj2;
    ComparableTerm obj2;
    ComparableTerm pred2;
    ComparableOntology ont2;
    
    public SimpleComparableTripleTest(String testName) {
        super(testName);
        ont = new SimpleComparableOntology("test");
        subj = new SimpleComparableTerm(ont, "subj", null);
        obj = new SimpleComparableTerm(ont, "obj", null);
        pred = new SimpleComparableTerm(ont, "pred", null);
        ont2 = new SimpleComparableOntology("test2");
        subj2 = new SimpleComparableTerm(ont, "subj2", null);
        obj2 = new SimpleComparableTerm(ont, "obj2", null);
        pred2 = new SimpleComparableTerm(ont, "pred2", null);
    }

    protected void setUp() throws Exception {
        cr = new ChangeEventRecorder();
        trip = new SimpleComparableTriple(ont, subj, obj, pred);
        trip.addChangeListener(cr);
    }

    protected void tearDown() throws Exception {
        trip.removeChangeListener(cr);
        trip = null;
        cr = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleComparableTripleTest.class);
        
        return suite;
    }

    /**
     * Test of compareTo method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        /*
         * Triples are equal only if they are from the same ontology and share the
         * same subject, object and predicate
         */
        
        //canonical
        assertTrue(trip.compareTo(trip) == 0);
        //equal
        SimpleComparableTriple trip2 = new SimpleComparableTriple(ont, subj, obj, pred);
        assertTrue(trip.compareTo(trip2) ==0);
        assertTrue(trip2.compareTo(trip) ==0);
        
        //not equal
        trip2 = new SimpleComparableTriple(ont2, subj, obj, pred);
        assertTrue(trip.compareTo(trip2) < 0);
        assertTrue(trip2.compareTo(trip) > 0);
        
        //not equal
        trip2 = new SimpleComparableTriple(ont, subj2, obj, pred);
        assertTrue(trip.compareTo(trip2) < 0);
        assertTrue(trip2.compareTo(trip) > 0);
        
        //not equal
        trip2 = new SimpleComparableTriple(ont, subj, obj2, pred);
        assertTrue(trip.compareTo(trip2) < 0);
        assertTrue(trip2.compareTo(trip) > 0);
        
        //not equal
        trip2 = new SimpleComparableTriple(ont, subj, obj, pred2);
        assertTrue(trip.compareTo(trip2) < 0);
        assertTrue(trip2.compareTo(trip) > 0);
    }

    /**
     * Test of equals method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        /*
         * Triples are equal only if they are from the same ontology and share the
         * same subject, object and predicate
         */
        
        //canonical
        assertTrue(trip.equals(trip));
        //equal
        Triple trip2 = new SimpleComparableTriple(ont, subj, obj, pred);
        assertTrue(trip.equals(trip2));
        assertTrue(trip2.equals(trip));
        
        //not equal
        trip2 = new SimpleComparableTriple(ont2, subj, obj, pred);
        assertTrue(! trip.equals(trip2));
        assertTrue(! trip2.equals(trip));
        
        //not equal
        trip2 = new SimpleComparableTriple(ont, subj2, obj, pred);
        assertTrue(! trip.equals(trip2));
        assertTrue(! trip2.equals(trip));
        
        //not equal
        trip2 = new SimpleComparableTriple(ont, subj, obj2, pred);
        assertTrue(! trip.equals(trip2));
        assertTrue(! trip2.equals(trip));
        
        //not equal
        trip2 = new SimpleComparableTriple(ont, subj, obj, pred2);
        assertTrue(! trip.equals(trip2));
        assertTrue(! trip2.equals(trip));
    }

    /**
     * Test of hashCode method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        SimpleComparableTriple trip2 = new SimpleComparableTriple(ont, subj, obj, pred);
        //canonical
        assertTrue(trip.hashCode() == trip2.hashCode());
        //equal
        assertTrue(trip.hashCode() == trip2.hashCode());
        
    }

    /**
     * Test of getName method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        //same as toString
        assertEquals(trip.toString(), trip.getName());
    }

    /**
     * Test of getSubject method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetSubject() {
        System.out.println("testGetSubject");
        
        assertEquals(subj, trip.getSubject());
    }

    /**
     * Test of getObject method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetObject() {
        System.out.println("testGetObject");
        
        assertEquals(obj, trip.getObject());
    }

    /**
     * Test of getPredicate method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetPredicate() {
        System.out.println("testGetPredicate");
        
        assertEquals(pred, trip.getPredicate());
    }

    /**
     * Test of addDescriptor method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testAddDescriptor() {
        System.out.println("testAddDescriptor");
        
        ComparableTerm term = new SimpleComparableTerm(ont, "foo", null);
        // cannot be null
        try{
            trip.addDescriptor(null);
            fail("Cannot be null");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            trip.addDescriptor(term);
            //should generate an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the correct type
            assertEquals(ComparableTriple.DESCRIPTOR, ce.getType());
            //old was
            assertNull(ce.getPrevious());
            //new is
            assertEquals(term, ce.getChange());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of removeDescriptor method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testRemoveDescriptor() {
        System.out.println("testRemoveDescriptor");
        
        ComparableTerm term = new SimpleComparableTerm(ont, "foo", null);
        
        try{
            trip.addDescriptor(term);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        // cannot be null
        try{
            trip.removeDescriptor(null);
            fail("Cannot be null");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            trip.removeDescriptor(term);
            //should generate an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the correct type
            assertEquals(ComparableTriple.DESCRIPTOR, ce.getType());
            //old was
            assertEquals(term, ce.getPrevious());
            //new is
            assertNull(ce.getChange());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of getDescriptors method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetDescriptors() {
        System.out.println("testGetDescriptors");
        
        // empty set at start
        assertNotNull(trip.getDescriptors());
        assertTrue(trip.getDescriptors().size() == 0);
    }

    /**
     * Test of setDescriptors method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testSetDescriptors() {
        System.out.println("testSetDescriptors");
        
        Set desc= new TreeSet();
        // can be null
        try{
            trip.setDescriptors(null);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            trip.setDescriptors(desc);
            //no event for this hibernate method
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of removeSynonym method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testRemoveSynonym() {
        System.out.println("testRemoveSynonym");
        
        //throws an exception.
        try{
            trip.removeSynonym("Synonym");
            fail("Should throw UnsupportedOperationException");
        }catch(UnsupportedOperationException ex){}
    }

    /**
     * Test of addSynonym method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testAddSynonym() {
        System.out.println("testAddSynonym");
        
        // throws an exception.
        try{
            trip.addSynonym("Synonym");
            fail("Should throw UnsupportedOperationException");
        }catch(UnsupportedOperationException ex){}
        
    }

    /**
     * Test of getSynonyms method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetSynonyms() {
        System.out.println("testGetSynonyms");
        
        //not null
        assertNotNull(trip.getSynonyms());
        // always an empty list.
        assertTrue(trip.getSynonyms().length == 0);
    }

    /**
     * Test of getOntology method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetOntology() {
        System.out.println("testGetOntology");
        
        assertEquals(ont, trip.getOntology());
    }

    /**
     * Test of getDescription method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetDescription() {
        System.out.println("testGetDescription");
        
        // always the empty string
        assertEquals("", trip.getDescription());
    }

    /**
     * Test of getAnnotation method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testGetAnnotation() {
        System.out.println("testGetAnnotation");
        
        assertEquals(RichAnnotation.EMPTY_ANNOTATION, trip.getAnnotation());
    }

    /**
     * Test of toString method, of class org.biojavax.ontology.SimpleComparableTriple.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = ont+":"+pred+"("+subj+","+obj+")";
        assertEquals(expected, trip.toString());
    }
    
}
