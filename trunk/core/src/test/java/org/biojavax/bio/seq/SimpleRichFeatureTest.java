package org.biojavax.bio.seq;

import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import junit.framework.TestCase;

import org.biojava.bio.symbol.SymbolList;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNamespace;
import org.biojavax.SimpleNote;
import org.biojavax.SimpleRichAnnotation;
import org.biojavax.ontology.ComparableTerm;

/**
 * Unit tests for SimpleRichFaeture.
 * @author Bubba Puryear
 */
public class SimpleRichFeatureTest extends TestCase {
    private SimpleRichFeature richFeature;

    /**
     * A quickie utility function to streamline some of the test code.
     * @param term the name of a term to get or create
     * @return said Term in all it's OO glory
     */
    private ComparableTerm makeTerm(String term) {
		return RichObjectFactory.getDefaultOntology().getOrCreateTerm(term);
	}

    protected void setUp() throws Exception {
        SimpleRichSequence parent = new SimpleRichSequence(new SimpleNamespace("testns"), "foo", "test_acc", 1, SymbolList.EMPTY_LIST, null);
        RichFeature.Template templ = new RichFeature.Template();
        templ.location = new SimpleRichLocation(new SimplePosition(1), 1);
        templ.sourceTerm = makeTerm("unknown");
        templ.typeTerm = makeTerm("CDS");
        templ.annotation = new SimpleRichAnnotation();
        this.richFeature = new SimpleRichFeature(parent, templ);
        this.richFeature.setRank(0);
    }
   
    /**
     * Test of equals method, of class org.biojavax.bio.seq.SimpleRichFeature.
     */
    public void testEquals() {
        //Two richFeatures are equal if they share the same parent, type term, source term and rank.
        assertFalse(this.richFeature.equals(new Object()));
        assertFalse(this.richFeature.equals(null));
        assertTrue(this.richFeature.equals(this.richFeature));

        SimpleRichFeature richFeature2 = null;
        try {
			richFeature2 = new SimpleRichFeature(this.richFeature.getParent(), this.richFeature.makeTemplate());
	        richFeature2.setRank(0);
		} catch (Exception e) {
			fail("Unexpected exception: "+e);
		}
        assertTrue(this.richFeature.equals(richFeature2));
        assertTrue(richFeature2.equals(this.richFeature));

        //should still be equal
        try{
            Set notes = new HashSet();
            notes.add(new SimpleNote(makeTerm("something"), "a value", 1));
            richFeature2.setNoteSet(notes);
            assertTrue(this.richFeature.equals(richFeature2));
            assertTrue(richFeature2.equals(this.richFeature));
        }catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }

        //should not be equal
        SimpleRichSequence otherParent = new SimpleRichSequence(new SimpleNamespace("anotherns"), "foo", "another_acc", 2, SymbolList.EMPTY_LIST, null);
        try {
			richFeature2 = new SimpleRichFeature(otherParent, this.richFeature.makeTemplate());
	        richFeature2.setRank(0);
		} catch (Exception e) {
            fail("Not expecting "+e.getClass().getName());
		}
        assertFalse(this.richFeature.equals(richFeature2));
        assertFalse(richFeature2.equals(this.richFeature));


        //should not be equal
        RichFeature.Template otherTempl = new RichFeature.Template();
        otherTempl.location = new SimpleRichLocation(new SimplePosition(1), 1);
		otherTempl.sourceTerm = makeTerm("unknown");
        otherTempl.typeTerm = makeTerm("promoter");
        otherTempl.annotation = new SimpleRichAnnotation();
        
        try {
			richFeature2 = new SimpleRichFeature(this.richFeature.getParent(), otherTempl);
	        richFeature2.setRank(0);
		} catch (Exception e) {
			fail("Unexpected exception: "+e);
		}
        assertFalse(this.richFeature.equals(richFeature2));
        assertFalse(richFeature2.equals(this.richFeature));
    }


    /**
     * Test of compareTo method, of class org.biojavax.bio.seq.SimpleRichFeature.
     */
    public void testCompareTo() {
        /*
         * RichFeatures are ordered by parent, type term, source term and rank.
         */

        assertTrue(this.richFeature.compareTo(this.richFeature) == 0);

        SimpleRichFeature richFeature2 = null;
        try {
        	richFeature2 = new SimpleRichFeature(this.richFeature.getParent(), this.richFeature.makeTemplate());
	        richFeature2.setRank(0);
        } catch (Exception e) {
			fail("Unexpected exception: "+e);
        }
        assertTrue(this.richFeature.compareTo(richFeature2) == 0);
        assertTrue(richFeature2.compareTo(this.richFeature) == 0);

        //should still be equal
        try{
            Set notes = new HashSet();
            notes.add(new SimpleNote(makeTerm("something"), "a value", 1));
            richFeature2.setNoteSet(notes);
            assertTrue(this.richFeature.compareTo(richFeature2) == 0);
            assertTrue(richFeature2.compareTo(this.richFeature) == 0);
        } catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }

        //should not be equal
        SimpleRichSequence otherParent = new SimpleRichSequence(new SimpleNamespace("anotherns"), "foo", "another_acc", 2, SymbolList.EMPTY_LIST, null);
        try {
        	richFeature2 = new SimpleRichFeature(otherParent, this.richFeature.makeTemplate());
	        richFeature2.setRank(0);
        } catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        assertTrue(this.richFeature.compareTo(richFeature2) > 0);
        assertTrue(richFeature2.compareTo(this.richFeature) < 0);


        //should not be equal
        RichFeature.Template otherTempl = new RichFeature.Template();
        otherTempl.location = new SimpleRichLocation(new SimplePosition(1), 1);
		otherTempl.sourceTerm = makeTerm("unknown");
        otherTempl.typeTerm = makeTerm("promoter");
        otherTempl.annotation = new SimpleRichAnnotation();
        try {
        	richFeature2 = new SimpleRichFeature(this.richFeature.getParent(), otherTempl);
	        richFeature2.setRank(0);
        } catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        assertTrue(this.richFeature.compareTo(richFeature2) < 0);
        assertTrue(richFeature2.compareTo(this.richFeature) > 0);


        //should not be equal
        otherTempl.typeTerm = makeTerm("CDS");
        otherTempl.sourceTerm = makeTerm("known");
        try {
        	richFeature2 = new SimpleRichFeature(this.richFeature.getParent(), otherTempl);
	        richFeature2.setRank(0);
        } catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        assertTrue(this.richFeature.compareTo(richFeature2) > 0);
        assertTrue(richFeature2.compareTo(this.richFeature) < 0);
    }

    /**
     * Test of hashCode method, of class org.biojavax.bio.SimpleBioEntry.
     */
    public void testHashCode() {
        SimpleRichFeature richFeature2 = null;
        try {
        	richFeature2 = new SimpleRichFeature(this.richFeature.getParent(), this.richFeature.makeTemplate());
	        richFeature2.setRank(0);
        } catch (Exception e) {
			fail("Unexpected exception: "+e);
        }
        assertTrue(this.richFeature.hashCode() == richFeature2.hashCode());

        //should still be equal
        try{
            Set notes = new HashSet();
            notes.add(new SimpleNote(makeTerm("something"), "a value", 1));
            richFeature2.setNoteSet(notes);
        } catch(Exception ex){
            fail("Not expecting "+ex.getClass().getName());
        }
        assertTrue(this.richFeature.hashCode() == richFeature2.hashCode());        
    }

    /**
     * A test that SimpleRichFeature behaves well with the TreeSet class.
     */
    public void testFeatureSet() {
    	Set set = new TreeSet();
    	set.add(this.richFeature);
    	assertTrue(set.contains(this.richFeature));
    	assertEquals(1, set.size());
    	set.remove(this.richFeature);
    	assertFalse(set.contains(this.richFeature));
    	assertEquals(0, set.size());
    }
 
}
