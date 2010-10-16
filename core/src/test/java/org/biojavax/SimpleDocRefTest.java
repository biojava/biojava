/*
 * SimpleDocRefTest.java
 * JUnit based test
 *
 * Created on November 11, 2005, 1:37 PM
 */

package org.biojavax;

import java.util.Collections;
import java.util.List;
import java.util.zip.Checksum;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojavax.utils.CRC64Checksum;


/**
 *
 * @author Mark Schreiber
 */
public class SimpleDocRefTest extends TestCase {
    SimpleDocRef ref;
    List authors;
    String location;
    String title;
    
    public SimpleDocRefTest(String testName) {
        super(testName);
        authors = Collections.singletonList(
                new SimpleDocRefAuthor("Hubert Hubertson", false, false));
        location = "Journal of Voodo Virology (7) 222-265";
        title = "ADE, myth or lie?";
    }

    protected void setUp() throws Exception {
        ref = new SimpleDocRef(authors, location, title);
    }

    protected void tearDown() throws Exception {
        ref = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleDocRefTest.class);
        
        return suite;
    }

    /**
     * Test of setRemark method, of class org.biojavax.SimpleDocRef.
     */
    public void testSetRemark() {
        System.out.println("testSetRemark");
        
        //should be able to do this
        try{
            ref.setRemark(null);
        }catch(Exception e){
            fail("Should be able to set remark to null without"+
                    e.getClass().getName());
        }
        try{
            String remark = "Remarkable!";
            ref.setRemark(remark);
            assertEquals(remark, ref.getRemark());
        }catch(Exception e){
            fail("Should be able to set remark without"+
                    e.getClass().getName());
        }
    }

    /**
     * Test of setTitle method, of class org.biojavax.SimpleDocRef.
     */
    
    //this method is now private. don't need to test it
//    public void testSetTitle() {
//        System.out.println("testSetTitle");
//        
//        //should be able to do this
//        try{
//            ref.setTitle(null);
//        }catch(Exception e){
//            fail("Should be able to set title to null without"+
//                    e.getClass().getName());
//        }
//        try{
//            String title = "Title";
//            ref.setTitle(title);
//            assertEquals(title, ref.getTitle());
//        }catch(Exception e){
//            fail("Should be able to set title without"+
//                    e.getClass().getName());
//        }
//    }

    /**
     * Test of setCrossref method, of class org.biojavax.SimpleDocRef.
     */
    public void testSetCrossref() {
        System.out.println("testSetCrossref");
        
        try{
            ref.setCrossref(null);
        }catch(Exception e){
            fail("Should be able to set crossref to null without"+
                    e.getClass().getName());
        }
        
        try{
            SimpleCrossRef xref = new SimpleCrossRef("another DB","AC123456", 1);
            ref.setCrossref(xref);
            assertEquals(xref, ref.getCrossref());
        }catch(Exception e){
            fail("Should be able to set crossref without"+
                    e.getClass().getName());
        }
    }

    /**
     * Test of getAuthors method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetAuthors() {
        System.out.println("testGetAuthors");
        
        assertNotNull(ref.getAuthors());
    }

    /**
     * Test of getAuthorList method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetAuthorList() {
        System.out.println("testGetAuthorList");
        
        assertEquals(DocRefAuthor.Tools.generateAuthorString(authors,true),
                ref.getAuthors());
    }

    /**
     * Test of getCRC method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetCRC() {
        System.out.println("testGetCRC");
        
        SimpleDocRef ref2 = new SimpleDocRef(authors, location, title);
        assertNotNull(ref.getCRC());
        assertEquals(ref.getCRC(), ref2.getCRC());
        
        StringBuffer sb = new StringBuffer();
        sb.append(ref.getAuthors());
        sb.append((ref.getTitle() ==null || ref.getTitle().equals(""))?"<undef>":ref.getTitle());
        sb.append((ref.getLocation() ==null || ref.getLocation().equals(""))?"<undef>":ref.getLocation());
        Checksum cs = new CRC64Checksum();
        cs.update(sb.toString().getBytes(), 0, sb.length());
        
        assertEquals(cs.toString(), ref.getCRC());
    }

    /**
     * Test of getRemark method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetRemark() {
        System.out.println("testGetRemark");
        //should be null until intitialized
        assertNull(ref.getRemark());
    }

    /**
     * Test of getCrossref method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetCrossref() {
        System.out.println("testGetCrossref");
        //should be null until intitialized
        assertNull(ref.getCrossref());
    }

    /**
     * Test of getLocation method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetLocation() {
        System.out.println("testGetLocation");
        
        assertNotNull(ref.getLocation());
        assertEquals(ref.getLocation(), location);
    }

    /**
     * Test of getTitle method, of class org.biojavax.SimpleDocRef.
     */
    public void testGetTitle() {
        System.out.println("testGetTitle");
        
        //should be equal to the set value
        assertEquals(this.title ,ref.getTitle());
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleDocRef.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        DocRef before = new SimpleDocRef(authors, "A", title);
        assertTrue(before.compareTo(ref) < 0);
        assertTrue(ref.compareTo(before) > 0);
        before = new SimpleDocRef(Collections.singletonList(
                new SimpleDocRefAuthor("A", false, false)), location, title);
        assertTrue(before.compareTo(ref) < 0);
        assertTrue(ref.compareTo(before) > 0);
        before = new SimpleDocRef(authors, location, "AAA");
        assertTrue(before.compareTo(ref) < 0);
        assertTrue(ref.compareTo(before) > 0);
        before = new SimpleDocRef(authors, location, null);
        assertTrue(before.compareTo(ref) < 0);
        assertTrue(ref.compareTo(before) > 0);
        
        DocRef equal = new SimpleDocRef(authors, location, title);
        assertTrue(ref.compareTo(ref) ==0);
        assertTrue(ref.compareTo(equal) ==0);
        assertTrue(equal.compareTo(ref) ==0);
    }

    /**
     * Test of equals method, of class org.biojavax.SimpleDocRef.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(ref.equals(ref));
        assertTrue(ref.equals(new SimpleDocRef(authors, location, title)));
        assertTrue(new SimpleDocRef(authors, location, title).equals(ref));
        
        assertFalse(new SimpleDocRef(authors, "A", title).equals(ref));
        assertFalse(ref.equals(new SimpleDocRef(authors, "A", title)));
        assertFalse(ref.equals(new SimpleDocRef(authors, location, "The long awaited biojava book!")));
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleDocRef.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        assertEquals(ref.hashCode(), new SimpleDocRef(authors, location, title).hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleDocRef.
     */
    public void testToString() {
        System.out.println("testToString");
        
        assertEquals("Hubert Hubertson; "+location, ref.toString());
    }
    
}
