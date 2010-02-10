/*
 * SimpleNamespaceTest.java
 * JUnit based test
 *
 * Created on 11 November 2005, 21:28
 */

package org.biojavax;

import java.net.URI;
import java.net.URISyntaxException;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.utils.ChangeVetoException;

/**
 *
 * @author Mark Schreiber
 */
public class SimpleNamespaceTest extends TestCase {
    SimpleNamespace ns;
    SimpleNamespace ns2;
    
    String name = "test space";
    String acronym = "by any other name";
    String authority = "respect mah authoritaaah!!!";
    String desc = "description";
    String uriString = "file:///foo.bar";
    
    public SimpleNamespaceTest(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
        ns = (SimpleNamespace)RichObjectFactory.getObject(SimpleNamespace.class, new Object[]{name});
        ns2 = (SimpleNamespace)RichObjectFactory.getObject(SimpleNamespace.class, new Object[]{name});
    }

    protected void tearDown() throws Exception {
        ns = null;
        ns2 = null;
        //clear the buffer
        RichObjectFactory.clearLRUCache(SimpleNamespace.class);
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleNamespaceTest.class);
        
        return suite;
    }

    /**
     * Test of setAcronym method, of class org.biojavax.SimpleNamespace.
     */
    public void testSetAcronym() {
        System.out.println("testSetAcronym");
        
        try{
            ns.setAcronym(acronym);
            assertEquals(acronym, ns.getAcronym());
        }catch(ChangeVetoException ex){
            fail("Was not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of setAuthority method, of class org.biojavax.SimpleNamespace.
     */
    public void testSetAuthority() {
        System.out.println("testSetAuthority");
        
        
        try{
            ns.setAuthority(authority);
            assertEquals(authority, ns.getAuthority());
        }catch(ChangeVetoException ex){
            fail("Was not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of setDescription method, of class org.biojavax.SimpleNamespace.
     */
    public void testSetDescription() {
        System.out.println("testSetDescription");
        
        
        try{
            ns.setDescription(desc);
            assertEquals(desc, ns.getDescription());
        }catch(ChangeVetoException ex){
            fail("Was not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of setURI method, of class org.biojavax.SimpleNamespace.
     */
    public void testSetURI() throws URISyntaxException{
        System.out.println("testSetURI");
        
        URI uri = new URI(uriString);
        try{
            ns.setURI(uri);
            assertEquals(uri, ns.getURI());
        }catch(ChangeVetoException ex){
            fail("Was not expecting "+ex.getClass().getName());
        }
    }

    /**
     * Test of getAcronym method, of class org.biojavax.SimpleNamespace.
     */
    public void testGetAcronym() {
        System.out.println("testGetAcronym");
        //System.out.println(ns.getAcronym());
        //should be set because of retreival from the LRU cache
        assertEquals(acronym, ns.getAcronym());
    }

    /**
     * Test of getAuthority method, of class org.biojavax.SimpleNamespace.
     */
    public void testGetAuthority() {
        System.out.println("testGetAuthority");
        
        assertEquals(authority,ns.getAuthority());
    }

    /**
     * Test of getDescription method, of class org.biojavax.SimpleNamespace.
     */
    public void testGetDescription() {
        System.out.println("testGetDescription");
        
        //should be null at startup
        //assertNull(ns.getDescription());
        try{
            ns.setDescription(desc);
        }catch(Exception ex){
            fail("Was not expecting "+ex.getClass().getName());
        }
        assertEquals(desc, ns.getDescription());
    }

    /**
     * Test of getName method, of class org.biojavax.SimpleNamespace.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        assertEquals(name, ns.getName());
    }

    /**
     * Test of getURI method, of class org.biojavax.SimpleNamespace.
     */
    public void testGetURI() {
        System.out.println("testGetURI");
        
       assertEquals(uriString, ns.getURI().toString());
    }

    /**
     * Test of compareTo method, of class org.biojavax.SimpleNamespace.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(ns.compareTo(RichObjectFactory.getDefaultNamespace()) >1);
        assertTrue(RichObjectFactory.getDefaultNamespace().compareTo(ns) <1);
        assertTrue(ns.compareTo(ns2) == 0);
    }

    /**
     * Test of equals method, of class org.biojavax.SimpleNamespace.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(ns == ns2);
        assertTrue(ns.equals(ns2));
        assertFalse(ns.equals(RichObjectFactory.getDefaultNamespace()));
    }

    /**
     * Test of hashCode method, of class org.biojavax.SimpleNamespace.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        assertTrue(ns.hashCode() == ns2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.SimpleNamespace.
     */
    public void testToString() {
        System.out.println("testToString");
        
        assertEquals(name, ns.toString());
    }
    
}
