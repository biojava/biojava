/*
 * SimpleComparableOntologyTest.java
 * JUnit based test
 *
 * Created on December 9, 2005, 9:59 AM
 */

package org.biojavax.ontology;

import java.util.NoSuchElementException;
import java.util.Set;
import java.util.TreeSet;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.ontology.AlreadyExistsException;
import org.biojava.ontology.Term;
import org.biojava.ontology.Triple;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.RichObjectFactory;


/**
 *
 * @author Mark Schreiber
 */
public class SimpleComparableOntologyTest extends TestCase {
    SimpleComparableOntology ont;
    String name;
    ComparableTerm term;
    Triple t3;
    
    ChangeListener.ChangeEventRecorder cr;
            
    public SimpleComparableOntologyTest(String testName) {
        super(testName);
        name = "a_test";
    }

    protected void setUp() throws Exception {
        cr = new ChangeListener.ChangeEventRecorder();
        ont = new SimpleComparableOntology(name);
        term = ont.getOrCreateTerm("foo");
        t3 = ont.createTriple(term, term, term, "t3", "");
        
        ont.addChangeListener(cr);
        
    }

    protected void tearDown() throws Exception {
        ont.removeChangeListener(cr);
        term = null;
        ont = null;
        cr = null;
        RichObjectFactory.clearLRUCache(SimpleComparableTerm.class);
        RichObjectFactory.clearLRUCache(SimpleComparableTriple.class);
        RichObjectFactory.clearLRUCache(SimpleComparableOntology.class);
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleComparableOntologyTest.class);
        
        return suite;
    }

    /**
     * Test of compareTo method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        // comparison is only by name
        assertTrue(ont.compareTo(ont) == 0);
        
        //equal
        SimpleComparableOntology ont2 = new SimpleComparableOntology(name);
        assertTrue(ont.compareTo(ont2) == 0);
        assertTrue(ont2.compareTo(ont) == 0);
        
        //not equal
        ont2 = new SimpleComparableOntology("xyz");
        assertTrue(ont.compareTo(ont2) < 1);
        assertTrue(ont2.compareTo(ont) > 1);
    }

    /**
     * Test of equals method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        //basic tests
        assertFalse(ont.equals(null));
        assertFalse(ont.equals(new Object()));
        
        // comparison is only by name
        assertTrue(ont.equals(ont));
        
        //equal
        SimpleComparableOntology ont2 = new SimpleComparableOntology(name);
        assertTrue(ont.equals(ont2));
        assertTrue(ont2.equals(ont));
        
        //not equal
        ont2 = new SimpleComparableOntology("xyz");
        assertTrue(! ont.equals(ont2));
        assertTrue(! ont2.equals(ont));
    }

    /**
     * Test of hashCode method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        SimpleComparableOntology ont2 = new SimpleComparableOntology(name);
        assertTrue(ont.hashCode() == ont2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testToString() {
        System.out.println("testToString");
        
        assertEquals(name, ont.toString());
    }

    /**
     * Test of containsTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testContainsTerm() {
        System.out.println("testContainsTerm");
        
        assertTrue(ont.containsTerm(term.getName()));
        assertFalse(ont.containsTerm("bar"));
        
    }

    /**
     * Test of getTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetTerm() {
        System.out.println("testGetTerm");
        
        assertEquals(term, ont.getTerm(term.getName()));
        try{
            ont.getTerm("foobar");
            fail("Expected NoSuchElementException");
        }catch(NoSuchElementException ex){}
    }

    /**
     * Test of getOrCreateTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetOrCreateTerm() {
        System.out.println("testGetOrCreateTerm");
        
        int size = ont.getTerms().size();
        assertEquals(term, ont.getOrCreateTerm(term.getName()));
        //size should still be the same as it was already there
        assertEquals(size, ont.getTerms().size());
        
        Term foo = ont.getOrCreateTerm("foobar");
        assertNotNull(foo);
        assertEquals(size+1, ont.getTerms().size());
        
    }

    /**
     * Test of getOrImportTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetOrImportTerm() {
        System.out.println("testGetOrImportTerm");
        
        int size = ont.getTerms().size();
        Term t = ont.getOrImportTerm(term);
        assertEquals(term, t);
        //size should still be the same as it was already there
        assertEquals(size, ont.getTerms().size());
        
        Term t2 = new SimpleComparableTerm(
                new SimpleComparableOntology("test2"),
                "hkhjj", null);
        
        t2 = ont.getOrImportTerm(t2);
        assertTrue(ont.containsTerm(t2.getName()));
        assertEquals(size+1, ont.getTerms().size());
    }

    /**
     * Test of createTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testCreateTerm() {
        System.out.println("testCreateTerm");
        
        try{
            Term term = ont.createTerm("bar","");
            assertEquals(term, ont.getTerm("bar"));
            
            assertNotNull(cr.getEvent());
            ChangeEvent ce = cr.getEvent();
            assertEquals(ComparableOntology.TERM, ce.getType());
            assertEquals(term, ce.getChange());
            assertNull(ce.getPrevious());
            
            try{
                ont.createTerm("bar", "");
                fail("Expected AlreadyExistsException");
            }catch(AlreadyExistsException ex){}
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of importTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testImportTerm() {
        System.out.println("testImportTerm");
        
        int size = ont.getTerms().size();
        try{
            ont.importTerm(term, term.getName());
        }catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
        }
        
        //size should still be the same as it was already there
        assertEquals(size, ont.getTerms().size());
        
        Term t2 = new SimpleComparableTerm(
                new SimpleComparableOntology("test2"),
                "hkhjj", null);
        
        try{
            ont.importTerm(t2, t2.getName());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        assertTrue(ont.containsTerm(t2.getName()));
        assertEquals(size+1, ont.getTerms().size());
    }

    /**
     * Test of createTriple method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testCreateTriple() {
        System.out.println("testCreateTriple");
        
        int size = ont.getTriples(null, null, null).size();
        Term obj = ont.getOrCreateTerm("foo");
        Term subj = ont.getOrCreateTerm("bar");
        Term pred = ont.getOrCreateTerm("is_a");
        
        try{
            Triple trip = ont.createTriple(subj, obj, pred, "Triple", "");
            assertTrue(ont.containsTriple(subj, obj, pred));
            assertEquals(size+1, ont.getTriples(null, null, null).size());
            
            
            assertNotNull(cr.getEvent());
            ChangeEvent ce = cr.getEvent();
            assertEquals(ComparableOntology.TRIPLE, ce.getType());
            assertEquals(trip, ce.getChange());
            assertNull(ce.getPrevious());
            
            
            //can't add it twice
            try{
                trip = ont.createTriple(subj, obj, pred, "Triple", "");
                fail("Expected AlreadyExistsException");
            }catch(AlreadyExistsException ex){}
            
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of deleteTerm method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testDeleteTerm() {
        System.out.println("testDeleteTerm");
        
        int termSize = ont.getTerms().size();
        int tripleSize = ont.getTriples(null, null, null).size();
        
        try{
            ont.deleteTerm(term);
            assertFalse(ont.containsTerm(term.getName()));
            assertEquals(termSize -1, ont.getTerms().size());
            //should have lost the triple too.
            assertEquals(tripleSize -1, ont.getTriples(null, null, null).size());
            
            assertNotNull(cr.getEvent());
            ChangeEvent ce = cr.getEvent();
            assertEquals(ComparableOntology.TERM, ce.getType());
            assertEquals(term, ce.getPrevious());
            assertNull(ce.getChange());
            
        }catch(Exception ex){
           fail("Not expecting "+ex.getClass().getName()); 
        }
    }

    /**
     * Test of getTriples method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetTriples() {
        System.out.println("testGetTriples");
        
        assertEquals(1, ont.getTriples(null, null, null).size());
        assertEquals(1, ont.getTriples(term, null, null).size());
        assertEquals(1, ont.getTriples(null, term, null).size());
        assertEquals(1, ont.getTriples(null, null, term).size());
        assertEquals(1, ont.getTriples(term, null, term).size());
        assertEquals(1, ont.getTriples(null, term, term).size());
        assertEquals(1, ont.getTriples(term, term, null).size());
        assertEquals(1, ont.getTriples(term, term, term).size());
        
        assertTrue(ont.getTriples(null, null, null).contains(t3));
    }

    /**
     * Test of setTripleSet method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testSetTripleSet() {
        System.out.println("testSetTripleSet");
        
        Set s = new TreeSet();
        try{
            ont.setTripleSet(s);
        }catch(ChangeVetoException ex){
            fail("Not expecting "+ex.getClass().getName()); 
        }
    }

    /**
     * Test of getTripleSet method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetTripleSet() {
        System.out.println("testGetTripleSet");
        
        assertNotNull(ont.getTripleSet());
        assertTrue(ont.getTripleSet().contains(t3));
        //should be writable
        try{
            Term term2 = ont.createTerm("foo2",null);
            ont.getTripleSet().add(new SimpleComparableTriple(
                    ont, term, term, 
                    (ComparableTerm)term2));
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName()); 
        }
    }

    /**
     * Test of getTerms method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetTerms() {
        System.out.println("testGetTerms");
        
        assertNotNull(ont.getTerms());
        assertEquals(1, ont.getTerms().size());
    }

    /**
     * Test of setTermSet method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testSetTermSet() {
        System.out.println("testSetTermSet");
        
        Set s = new TreeSet();
        try{
            ont.setTermSet(s);
        }catch(ChangeVetoException ex){
            fail("Not expecting "+ex.getClass().getName()); 
        }
    }

    /**
     * Test of getTermSet method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetTermSet() {
        System.out.println("testGetTermSet");
        
        assertNotNull(ont.getTermSet());
        assertTrue(ont.getTermSet().contains(term));
        //should be writable
        try{
            Term term2 = ont.createTerm("foo2",null);
            ont.getTripleSet().add(term2);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName()); 
        }
    }

    /**
     * Test of containsTriple method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testContainsTriple() {
        System.out.println("testContainsTriple");
        
        assertTrue(ont.containsTriple(t3.getSubject(), t3.getObject(), t3.getPredicate()));
        assertFalse(ont.containsTriple(t3.getSubject(), t3.getObject(), ont.getOrCreateTerm("another_term")));
    }

    /**
     * Test of createVariable method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testCreateVariable() {
        System.out.println("testCreateVariable");
        
        try{
            ont.createVariable("name", "desc");
            fail("Should throw an Exception");
        }catch(Exception ex){}
    }

    /**
     * Test of getDescription method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetDescription() {
        System.out.println("testGetDescription");
        
        // null until set
        assertNull(ont.getDescription());
    }

    /**
     * Test of setDescription method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testSetDescription() {
        System.out.println("testSetDescription");
        String desc = "desc";
        try{
            ont.setDescription(desc);
            assertEquals(desc, ont.getDescription());
            
            assertNotNull(cr.getEvent());
            ChangeEvent ce = cr.getEvent();
            assertEquals(ComparableOntology.DESCRIPTION, ce.getType());
            assertEquals(desc, ce.getChange());
            assertNull(ce.getPrevious());
            
        }catch(ChangeVetoException ex){
            fail("Not expecting "+ex.getClass().getName()); 
        }
    }

    /**
     * Test of getName method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        assertEquals(name, ont.getName());
    }

    /**
     * Test of getOps method, of class org.biojavax.ontology.SimpleComparableOntology.
     */
    public void testGetOps() {
        System.out.println("testGetOps");
        
        // should be not null
        assertNotNull(ont.getOps());
    }
    
}
