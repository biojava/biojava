/*
 * SimpleBioEntryTest.java
 * JUnit based test
 *
 * Created on November 30, 2005, 3:09 PM
 */

package org.biojavax.bio;

import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener.ChangeEventRecorder;
import org.biojavax.Comment;
import org.biojavax.DocRefAuthor;
import org.biojavax.Namespace;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleComment;
import org.biojavax.SimpleCrossRef;
import org.biojavax.SimpleDocRef;
import org.biojavax.SimpleDocRefAuthor;
import org.biojavax.SimpleNamespace;
import org.biojavax.SimpleNote;
import org.biojavax.SimpleRankedCrossRef;
import org.biojavax.SimpleRankedDocRef;
import org.biojavax.bio.taxa.NCBITaxon;
import org.biojavax.bio.taxa.SimpleNCBITaxon;
import org.biojavax.ontology.ComparableTerm;

/**
 *
 * @author Mark Schreiber
 * @autor George Waldon
 */
public class SimpleBioEntryTest extends TestCase {
    SimpleBioEntry be;
    String name;
    Namespace ns;
    String acc;
    int version;
    ChangeEventRecorder cr;
    
    public SimpleBioEntryTest(String testName) {
        super(testName);
        name = "aaa";
        ns = RichObjectFactory.getDefaultNamespace();
        acc = "test_acc";
        version = 0;
        cr = new ChangeEventRecorder();
    }
    
    protected void setUp() throws Exception {
        be = new SimpleBioEntry(ns, name, acc, version);
        cr = new ChangeEventRecorder();
        be.addChangeListener(cr);
    }
    
    protected void tearDown() throws Exception {
        be.removeChangeListener(cr);
        be = null;
        cr = null;
    }
    
    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleBioEntryTest.class);
        
        return suite;
    }
    
    /**
     * Test of getRankedCrossRefs method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetRankedCrossRefs() {
        System.out.println("testGetRankedCrossRefs");
        
        //at the start should be not-null and empty;
        assertNotNull(be.getRankedCrossRefs());
        assertEquals(0, be.getRankedCrossRefs().size());
        //hibernate needs it to be writable
        try{
            be.getRankedCrossRefs().add(new Object());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of setTaxon method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testSetTaxon() {
        System.out.println("testSetTaxon");
        
        NCBITaxon tax = new SimpleNCBITaxon(1621);
        
        try{
            be.addChangeListener(cr);
            be.setTaxon(tax);
            
            //should have fired an event;
            assertNotNull(cr.getEvent());
            //type should be TAXON
            assertEquals(BioEntry.TAXON, cr.getEvent().getType());
            //old should be null
            assertNull(cr.getEvent().getPrevious());
            //new should be tax
            assertEquals(tax, cr.getEvent().getChange());
            
            //tax and get taxon should be equal
            assertEquals(tax, be.getTaxon());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of getAnnotation & getRichAnnotation methods of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetAnnotation() {
        System.out.println("testGetAnnotation");
        
        //should be not null!
        assertNotNull(be.getAnnotation());
        assertNotNull(be.getRichAnnotation());
    }
    
    /**
     * Test of getNoteSet method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetNoteSet() {
        System.out.println("testGetNoteSet");
        
        //should not be null;
        assertNotNull(be.getNoteSet());
        //hibernate needs it to be writable
        try{
            be.getNoteSet().add(new Object());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of setNoteSet method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testSetNoteSet() {
        System.out.println("testSetNoteSet");
        Set notes = new TreeSet();
        notes.add(new SimpleNote(
                RichObjectFactory.getDefaultOntology().getOrCreateTerm("foo"),
                "bar", 0));
        try{
            be.setNoteSet(notes);
            
            //doesn't generate a change event, should it??
            
            //should get notes back.
            assertEquals(notes, be.getNoteSet());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of getComments method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetComments() {
        System.out.println("testGetComments");
        
        //should be not-null and empty at start
        Set s = be.getComments();
        assertNotNull(s);
        assertEquals(0, s.size());
    }
    
    /**
     * Test of getRankedDocRefs method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetRankedDocRefs() {
        System.out.println("testGetRankedDocRefs");
        
        //should be not-null and empty at the start
        Set s = be.getRankedDocRefs();
        assertNotNull(s);
        assertEquals(0, s.size());
    }
    
    /**
     * Test of getRelationships method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetRelationships() {
        System.out.println("testGetRelationships");
        
        //should be not-null and empty at the start
        Set s = be.getRelationships();
        assertNotNull(s);
        assertEquals(0, s.size());
    }
    
    /**
     * Test of setIdentifier method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testSetIdentifier() {
        System.out.println("testSetIdentifier");
        
        String id = "new id";
        try{
            be.setIdentifier(id);
            //should get back id
            assertEquals(id, be.getIdentifier());
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.IDENTIFIER, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now id
            assertEquals(id, ce.getChange());
            
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of setDivision method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testSetDivision() {
        System.out.println("testSetDivision");
        
        String div = "new div";
        try{
            be.setDivision(div);
            //should get back div
            assertEquals(div, be.getDivision());
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.DIVISION, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now div
            assertEquals(div, ce.getChange());
            
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of setDescription method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testSetDescription() {
        System.out.println("testSetDescription");
        
        String d = "new desc";
        try{
            be.setDescription(d);
            //should get back d
            assertEquals(d, be.getDescription());
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.DESCRIPTION, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now d
            assertEquals(d, ce.getChange());
            
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
    /**
     * Test of getAccession method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetAccession() {
        System.out.println("testGetAccession");
        
        // should be acc
        assertEquals(this.acc, be.getAccession());
    }
    
    /**
     * Test of getDescription method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetDescription() {
        System.out.println("testGetDescription");
        
        // should be null at start
        assertNull(be.getDescription());
    }
    
    /**
     * Test of getDivision method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetDivision() {
        System.out.println("testGetDivision");
        
        // should be null at start
        assertNull(be.getDivision());
    }
    
    /**
     * Test of getIdentifier method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetIdentifier() {
        System.out.println("testGetIdentifier");
        
        // should be null at start
        assertNull(be.getIdentifier());
    }
    
    /**
     * Test of getName method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        // should be name
        assertEquals(this.name, be.getName());
    }
    
    /**
     * Test of getNamespace method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetNamespace() {
        System.out.println("testGetNamespace");
        
        // should be ns
        assertEquals(this.ns, be.getNamespace());
    }
    
    /**
     * Test of getTaxon method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetTaxon() {
        System.out.println("testGetTaxon");
        
        // should be null at start
        assertNull(be.getTaxon());
    }
    
    /**
     * Test of getVersion method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testGetVersion() {
        System.out.println("testGetVersion");
        
        // should be version
        assertEquals(this.version, be.getVersion());
    }
    
    /**
     * Test of equals method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        //Two bioentries are equal if they share the same namespace, name,
        //accession and version.
        assertFalse(be.equals(new Object()));
        assertFalse(be.equals(null));
        assertTrue(be.equals(be));
        
        BioEntry be2 = new SimpleBioEntry(be.getNamespace(), be.getName(),
                be.getAccession(), be.getVersion());
        assertTrue(be.equals(be2));
        assertTrue(be2.equals(be));
        
        //should still be equal
        try{
            be2.setDescription("test");
            assertTrue(be.equals(be2));
            assertTrue(be2.equals(be));
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        //should not be equal
        
        be2 = new SimpleBioEntry(
                new SimpleNamespace("new"), be.getName(),
                be.getAccession(), be.getVersion());
        assertFalse(be.equals(be2));
        assertFalse(be2.equals(be));
        
        
        //should not be equal
        
        be2 = new SimpleBioEntry(
                be.getNamespace(), "different",
                be.getAccession(), be.getVersion());
        assertFalse(be.equals(be2));
        assertFalse(be2.equals(be));
        
        
        //should not be equal
        
        be2 = new SimpleBioEntry(
                be.getNamespace(), be.getName(),
                "different", be.getVersion());
        assertFalse(be.equals(be2));
        assertFalse(be2.equals(be));
        
        
        //should not be equal
        be2 = new SimpleBioEntry(
                be.getNamespace(), be.getName(),
                be.getAccession(), be.getVersion()+1);
        assertFalse(be.equals(be2));
        assertFalse(be2.equals(be));
        
    }
    
    /**
     * Test of compareTo method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        /*
         * Bioentries are ordered first by namespace, then name, accession, and
         * finally version.
         */
        
        assertTrue(be.compareTo(be) == 0);
        
        BioEntry be2 = new SimpleBioEntry(be.getNamespace(), be.getName(),
                be.getAccession(), be.getVersion());
        assertTrue(be.compareTo(be2) == 0);
        assertTrue(be2.compareTo(be) == 0);
        
        //should still be equal
        try{
            be2.setDescription("test");
            assertTrue(be.compareTo(be2) == 0);
            assertTrue(be2.compareTo(be) == 0);
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        
        //should not be equal
        be2 = new SimpleBioEntry(
                new SimpleNamespace("new"), be.getName(),
                be.getAccession(), be.getVersion());
        assertTrue(be.compareTo(be2) < 0);
        assertTrue(be2.compareTo(be) > 0);
        
        
        //should not be equal
        be2 = new SimpleBioEntry(
                be.getNamespace(), "different",
                be.getAccession(), be.getVersion());
        assertTrue(be.compareTo(be2) < 0);
        assertTrue(be2.compareTo(be) > 0);
        
        
        //should not be equal
        be2 = new SimpleBioEntry(
                be.getNamespace(), be.getName(),
                "different", be.getVersion());
        assertTrue(be.compareTo(be2) > 0);
        assertTrue(be2.compareTo(be) < 0);
        
        
        //should not be equal
        be2 = new SimpleBioEntry(
                be.getNamespace(), be.getName(),
                be.getAccession(), be.getVersion()+1);
        assertTrue(be.compareTo(be2) < 0);
        assertTrue(be2.compareTo(be) > 0);
    }
    
    /**
     * Test of hashCode method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        BioEntry be2 = new SimpleBioEntry(
                be.getNamespace(), be.getName(),
                be.getAccession(), be.getVersion());
        
        assertTrue(be.hashCode() == be2.hashCode());
        
        //should still be equal
        try{
            be2.setDescription("test");
            assertTrue(be.hashCode() == be2.hashCode());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }       
    }
    
    /**
     * Test of toString method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = be.getNamespace()+":"+be.getName()+"/"+be.getAccession()+"."+be.getVersion(); 
        assertEquals(expected, be.toString());
    }
    
    /**
     * Test of addRankedCrossRef method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testAddRankedCrossRef() {
        System.out.println("testAddRankedCrossRef");
        
        //should not be able to add null
        try{
            be.addRankedCrossRef(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        RankedCrossRef xref = new SimpleRankedCrossRef(
                new SimpleCrossRef("dbname", "AC123456", 1), 0);
        try{
            be.addRankedCrossRef(xref);
            assertTrue(be.getRankedCrossRefs().contains(xref));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.RANKEDCROSSREF, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now xref
            assertEquals(xref, ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of removeRankedCrossRef method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testRemoveRankedCrossRef() {
        System.out.println("testRemoveRankedCrossRef");
        
        //first add one
        RankedCrossRef xref = new SimpleRankedCrossRef(
                new SimpleCrossRef("dbname", "AC123456", 1), 0);
        try{
            be.addRankedCrossRef(xref);
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        //cannot remove null
        try{
            be.removeRankedCrossRef(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        try{
            be.removeRankedCrossRef(xref);
            assertFalse(be.getRankedCrossRefs().contains(xref));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.RANKEDCROSSREF, ce.getType());
            //was xref
            assertEquals(xref, ce.getPrevious());
            //now null
            assertNull(ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of addRankedDocRef method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testAddRankedDocRef() {
        System.out.println("testAddRankedDocRef");
        
        //should not be able to add null
        try{
            be.addRankedDocRef(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        DocRefAuthor author = new SimpleDocRefAuthor("Hemmingway");
        RankedDocRef ref = new SimpleRankedDocRef(
                new SimpleDocRef(Collections.singletonList(author), "a book", "a title"),
                new Integer(1), new Integer(10), 0);
        try{
            be.addRankedDocRef(ref);
            assertTrue(be.getRankedDocRefs().contains(ref));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.RANKEDDOCREF, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now xref
            assertEquals(ref, ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of removeRankedDocRef method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testRemoveRankedDocRef() {
        System.out.println("testRemoveRankedDocRef");
        
        //should not be able to remove null
        try{
            be.removeRankedDocRef(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        DocRefAuthor author = new SimpleDocRefAuthor("Hemmingway");
        RankedDocRef ref = new SimpleRankedDocRef(
                new SimpleDocRef(Collections.singletonList(author), "a book", "a title"),
                new Integer(1), new Integer(10), 0);
        
        
        //first add one
        try{
            be.addRankedDocRef(ref);
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            be.removeRankedDocRef(ref);
            assertFalse(be.getRankedCrossRefs().contains(ref));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.RANKEDDOCREF, ce.getType());
            //was ref
            assertEquals(ref, ce.getPrevious());
            //now null
            assertNull(ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of addComment method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testAddComment() {
        System.out.println("testAddComment");
        
        //should not be able to add null
        try{
            be.addComment(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        Comment com = new SimpleComment("comment", 0);
        try{
            be.addComment(com);
            assertTrue(be.getComments().contains(com));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.COMMENT, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now com
            assertEquals(com, ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of removeComment method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testRemoveComment() {
        System.out.println("testRemoveComment");
        
        //should not be able to remove null
        try{
            be.removeComment(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        Comment com = new SimpleComment("comment", 0);
        
        
        //first add one
        try{
            be.addComment(com);
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            be.removeComment(com);
            assertFalse(be.getComments().contains(com));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.COMMENT, ce.getType());
            //was com
            assertEquals(com, ce.getPrevious());
            //now null
            assertNull(ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of addRelationship method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testAddRelationship() {
        System.out.println("testAddRelationship");
        
        //should not be able to add null
        try{
            be.addRelationship(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        BioEntry be2 = new SimpleBioEntry(
                be.getNamespace(), "different",
                be.getAccession(), be.getVersion());
        ComparableTerm term =
                RichObjectFactory.getDefaultOntology().getOrCreateTerm("foo");
        BioEntryRelationship rel = new SimpleBioEntryRelationship(
                be,be2, term, new Integer(0));
        try{
            be.addRelationship(rel);
            assertTrue(be.getRelationships().contains(rel));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.RELATIONS, ce.getType());
            //was null
            assertNull(ce.getPrevious());
            //now rel
            assertEquals(rel, ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of removeRelationship method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testRemoveRelationship() {
        System.out.println("testRemoveRelationship");
        
        //should not be able to remove null
        try{
            be.removeRelationship(null);
            fail("Expected IllegalArgumentException");
        }catch(IllegalArgumentException ex){}
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
        
        BioEntry be2 = new SimpleBioEntry(
                be.getNamespace(), "different",
                be.getAccession(), be.getVersion());
        ComparableTerm term =
                RichObjectFactory.getDefaultOntology().getOrCreateTerm("foo");
        BioEntryRelationship rel = new SimpleBioEntryRelationship(
                be,be2, term, new Integer(0));
        
        //first add one
        try{
            be.addRelationship(rel);
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
        }
        
        try{
            be.removeRelationship(rel);
            assertFalse(be.getRelationships().contains(rel));
            
            //should have generated an event
            ChangeEvent ce = cr.getEvent();
            assertNotNull(ce);
            //of the right type
            assertEquals(BioEntry.RELATIONS, ce.getType());
            //was rel
            assertEquals(rel, ce.getPrevious());
            //now null
            assertNull(ce.getChange());
        }
         catch(Exception ex){
             fail("Not expecting "+ex.getClass().getName());
         }
    }
    
    /**
     * Test of setRankedCrossRefs method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testSetRankedCrossRefs() {
        System.out.println("testSetRankedCrossRefs");
        
        Set s = new TreeSet();
        s.add(new SimpleRankedCrossRef(new SimpleCrossRef("dbname", "AC123456", 1), 0));
        try{
            be.setRankedCrossRefs(s);
            assertEquals(s, be.getRankedCrossRefs());
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
    }
    
}
