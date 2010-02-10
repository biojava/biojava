/*
 * SimpleNCBITaxonNameTest.java
 * JUnit based test
 *
 * Created on December 12, 2005, 3:49 PM
 */

package org.biojavax.bio.taxa;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;


/**
 *
 * @author Mark Schreiber
 */
public class SimpleNCBITaxonNameTest extends TestCase {
    
    SimpleNCBITaxonName tax;
    
    String nameClass;
    String name;
    
    public SimpleNCBITaxonNameTest(String testName) {
        super(testName);
        nameClass = "abc";
        name = "aaa";
    }

    protected void setUp() throws Exception {
        
        tax = new SimpleNCBITaxonName(nameClass, name);
       
    }

    protected void tearDown() throws Exception {
        tax = null;
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(SimpleNCBITaxonNameTest.class);
        
        return suite;
    }

    /**
     * Test of setNameClass method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testSetNameClass() {
        System.out.println("testSetNameClass");
        
        String newNameClass = "newClass";
        tax.setNameClass(newNameClass);
        assertEquals(newNameClass, tax.getNameClass());
    }

    /**
     * Test of getNameClass method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testGetNameClass() {
        System.out.println("testGetNameClass");
        
        assertEquals(nameClass, tax.getNameClass());
    }

    /**
     * Test of setName method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testSetName() {
        System.out.println("testSetName");
        
        String newName = "newName";
        tax.setName(newName);
        assertEquals(newName, tax.getName());
    }

    /**
     * Test of getName method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testGetName() {
        System.out.println("testGetName");
        
        assertEquals(name, tax.getName());
    }

    /**
     * Test of equals method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testEquals() {
        System.out.println("testEquals");
        
        assertTrue(tax.equals(tax));
        assertFalse(tax.equals(new Object()));
        
        SimpleNCBITaxonName tax2 = new SimpleNCBITaxonName(nameClass,name);
        
        //equal
        assertTrue(tax.equals(tax2));
        assertTrue(tax2.equals(tax));
        
        //not equal
        tax2 = new SimpleNCBITaxonName("foo",name);
        assertTrue(! tax.equals(tax2));
        assertTrue(! tax2.equals(tax));
        
        //not equal
        tax2 = new SimpleNCBITaxonName(nameClass,"foo");
        assertTrue(! tax.equals(tax2));
        assertTrue(! tax2.equals(tax));
        
        //not equal
        tax2 = new SimpleNCBITaxonName("foo","foo");
        assertTrue(! tax.equals(tax2));
        assertTrue(! tax2.equals(tax));
    }

    /**
     * Test of compareTo method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testCompareTo() {
        System.out.println("testCompareTo");
        
        assertTrue(tax.compareTo(tax) ==0);
        
        
        SimpleNCBITaxonName tax2 = new SimpleNCBITaxonName(nameClass,name);
        
        //equal
        assertTrue(tax.compareTo(tax2) == 0);
        assertTrue(tax2.compareTo(tax) == 0);
        
        //not equal
        tax2 = new SimpleNCBITaxonName("foo",name);
        assertTrue(tax.compareTo(tax2) < 1);
        assertTrue(tax2.compareTo(tax) > 1);
        
        //not equal
        tax2 = new SimpleNCBITaxonName(nameClass,"foo");
        assertTrue( tax.compareTo(tax2) < 1);
        assertTrue( tax2.compareTo(tax) > 1);
        
        //not equal
        tax2 = new SimpleNCBITaxonName("foo","foo");
        assertTrue( tax.compareTo(tax2) < 1);
        assertTrue( tax2.compareTo(tax) > 1);
    }

    /**
     * Test of hashCode method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testHashCode() {
        System.out.println("testHashCode");
        
        SimpleNCBITaxonName tax2 = new SimpleNCBITaxonName(nameClass,name);
        assertEquals(tax.hashCode(), tax2.hashCode());
    }

    /**
     * Test of toString method, of class org.biojavax.bio.taxa.SimpleNCBITaxonName.
     */
    public void testToString() {
        System.out.println("testToString");
        
        String expected = nameClass+":"+name;
        assertEquals(expected, tax.toString());
    }
    
}
