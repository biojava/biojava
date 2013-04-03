/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;


/**
 *
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 */
public class HetatomImplTest extends TestCase{

    int bigTestNumber = 60000;

    public HetatomImplTest() {

    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

//    /**
//     * Test of has3D method, of class HetatomImpl.
//     */
//    @Test
//    public void testHas3D() {
//        System.out.println("has3D");
//        HetatomImpl instance = new HetatomImpl();
//        boolean expResult = false;
//        boolean result = instance.has3D();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setPDBFlag method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetPDBFlag() {
//        System.out.println("setPDBFlag");
//        boolean flag = false;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setPDBFlag(flag);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPDBCode method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetPDBCode() {
//        System.out.println("getPDBCode");
//        HetatomImpl instance = new HetatomImpl();
//        String expResult = "";
//        String result = instance.getPDBCode();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setPDBCode method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetPDBCode() {
//        System.out.println("setPDBCode");
//        String pdb = "";
//        HetatomImpl instance = new HetatomImpl();
//        instance.setPDBCode(pdb);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setPDBName method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetPDBName() throws Exception {
//        System.out.println("setPDBName");
//        String s = "";
//        HetatomImpl instance = new HetatomImpl();
//        instance.setPDBName(s);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPDBName method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetPDBName() {
//        System.out.println("getPDBName");
//        HetatomImpl instance = new HetatomImpl();
//        String expResult = "";
//        String result = instance.getPDBName();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of addAtom method, of class HetatomImpl.
//     */
//    @Test
//    public void testAddAtom() {
//        System.out.println("addAtom");
//        Atom atom = null;
//        HetatomImpl instance = new HetatomImpl();
//        instance.addAtom(atom);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of clearAtoms method, of class HetatomImpl.
//     */
//    @Test
//    public void testClearAtoms() {
//        System.out.println("clearAtoms");
//        HetatomImpl instance = new HetatomImpl();
//        instance.clearAtoms();
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of size method, of class HetatomImpl.
//     */
//    @Test
//    public void testSize() {
//        System.out.println("size");
//        HetatomImpl instance = new HetatomImpl();
//        int expResult = 0;
//        int result = instance.size();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAtoms method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetAtoms() {
//        System.out.println("getAtoms");
//        HetatomImpl instance = new HetatomImpl();
//        List expResult = null;
//        List result = instance.getAtoms();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setAtoms method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetAtoms() {
//        System.out.println("setAtoms");
//        List<Atom> atoms = null;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setAtoms(atoms);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAtom method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetAtom_String() throws Exception {
//        System.out.println("getAtom");
//        String name = "";
//        HetatomImpl instance = new HetatomImpl();
//        Atom expResult = null;
//        Atom result = instance.getAtom(name);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAtomByPDBname method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetAtomByPDBname() throws Exception {
//        System.out.println("getAtomByPDBname");
//        String name = "";
//        HetatomImpl instance = new HetatomImpl();
//        Atom expResult = null;
//        Atom result = instance.getAtomByPDBname(name);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAtom method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetAtom_int() throws Exception {
//        System.out.println("getAtom");
//        int position = 0;
//        HetatomImpl instance = new HetatomImpl();
//        Atom expResult = null;
//        Atom result = instance.getAtom(position);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of hasAtom method, of class HetatomImpl.
//     */
//    @Test
//    public void testHasAtom() {
//        System.out.println("hasAtom");
//        String fullName = "";
//        HetatomImpl instance = new HetatomImpl();
//        boolean expResult = false;
//        boolean result = instance.hasAtom(fullName);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getType method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetType() {
//        System.out.println("getType");
//        HetatomImpl instance = new HetatomImpl();
//        String expResult = "";
//        String result = instance.getType();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of toString method, of class HetatomImpl.
//     */
//    @Test
//    public void testToString() {
//        System.out.println("toString");
//        HetatomImpl instance = new HetatomImpl();
//        String expResult = "";
//        String result = instance.toString();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of hasAminoAtoms method, of class HetatomImpl.
//     */
//    @Test
//    public void testHasAminoAtoms() {
//        System.out.println("hasAminoAtoms");
//        HetatomImpl instance = new HetatomImpl();
//        boolean expResult = false;
//        boolean result = instance.hasAminoAtoms();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setProperties method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetProperties() {
//        System.out.println("setProperties");
//        Map<String, Object> props = null;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setProperties(props);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getProperties method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetProperties() {
//        System.out.println("getProperties");
//        HetatomImpl instance = new HetatomImpl();
//        Map expResult = null;
//        Map result = instance.getProperties();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setProperty method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetProperty() {
//        System.out.println("setProperty");
//        String key = "";
//        Object value = null;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setProperty(key, value);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getProperty method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetProperty() {
//        System.out.println("getProperty");
//        String key = "";
//        HetatomImpl instance = new HetatomImpl();
//        Object expResult = null;
//        Object result = instance.getProperty(key);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of iterator method, of class HetatomImpl.
//     */
//    @Test
//    public void testIterator() {
//        System.out.println("iterator");
//        HetatomImpl instance = new HetatomImpl();
//        Iterator expResult = null;
//        Iterator result = instance.iterator();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of clone method, of class HetatomImpl.
//     */
//    @Test
//    public void testClone() {
//        System.out.println("clone");
//        HetatomImpl instance = new HetatomImpl();
//        Object expResult = null;
//        Object result = instance.clone();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setParent method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetParent() {
//        System.out.println("setParent");
//        Chain parent = null;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setParent(parent);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getParent method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetParent() {
//        System.out.println("getParent");
//        HetatomImpl instance = new HetatomImpl();
//        Chain expResult = null;
//        Chain result = instance.getParent();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getId method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetId() {
//        System.out.println("getId");
//        HetatomImpl instance = new HetatomImpl();
//        long expResult = 0L;
//        long result = instance.getId();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setId method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetId() {
//        System.out.println("setId");
//        long id = 0L;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setId(id);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getChemComp method, of class HetatomImpl.
//     */
//    @Test
//    public void testGetChemComp() {
//        System.out.println("getChemComp");
//        HetatomImpl instance = new HetatomImpl();
//        ChemComp expResult = null;
//        ChemComp result = instance.getChemComp();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setChemComp method, of class HetatomImpl.
//     */
//    @Test
//    public void testSetChemComp() {
//        System.out.println("setChemComp");
//        ChemComp cc = null;
//        HetatomImpl instance = new HetatomImpl();
//        instance.setChemComp(cc);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
    /**
     * Test of setChain method, of class HetatomImpl.
     */
    @Test
    public void testSetGetChain() {
//        System.out.println("setGetChain");
        Chain chain = new ChainImpl();
        chain.setChainID("A");
        HetatomImpl instance = new HetatomImpl();
        instance.setChain(chain);
        Chain expResult = chain;
        Chain result = instance.getChain();
        assertEquals(expResult, result);
    }

    /**
     * Test of getChainId method, of class HetatomImpl.
     */
    @Test
    public void testGetChainId() {
//        System.out.println("getChainId");
        Chain chain = new ChainImpl();
        chain.setChainID("A");
        HetatomImpl instance = new HetatomImpl();
        instance.setChain(chain);
        String expResult = "A";
        String result = instance.getChainId();
        assertEquals(expResult, result);
    }

    /**
     * Test of getSeqNum method, of class HetatomImpl.
     */
    @Test
    public void testSetGetResidueNumber() {
//        System.out.println("setGetResidueNumber");
        ResidueNumber residueNumber = new ResidueNumber("A", 42, ' ');
        HetatomImpl instance = new HetatomImpl();
        instance.setResidueNumber(residueNumber);
        ResidueNumber expResult = residueNumber;
        ResidueNumber result = instance.getResidueNumber();
        assertEquals(expResult, result);

    }

    @Test
    public void testGetResidueNumberUsage() {
//        System.out.println("testGetResidueNumberUsage");
        List<Group> resNumgroups = new ArrayList<Group>();

        for (int i = 0; i < bigTestNumber; i++) {
            ResidueNumber resNum = new ResidueNumber("A", i, ' ');
            HetatomImpl hetAtom = new HetatomImpl();
            hetAtom.setResidueNumber(resNum);
            resNumgroups.add(hetAtom);
        }

        List<Integer> integers = new ArrayList<Integer>();

        for (Group group : resNumgroups) {
            ResidueNumber resnum = group.getResidueNumber();
            integers.add(resnum.getSeqNum());
        }
        assertEquals(bigTestNumber, integers.size());
    }

    @Test
    public void testSetResidueNumberUsage() {
     
        List<Group> resNumgroups = new ArrayList<Group>();

        for (int i = 0; i < bigTestNumber; i++) {
            ResidueNumber resNum = new ResidueNumber("A", i, ' ');
            HetatomImpl hetAtom = new HetatomImpl();
            hetAtom.setResidueNumber(resNum);
            resNumgroups.add(hetAtom);
        }
        int groupsSize = resNumgroups.size();
        assertEquals(bigTestNumber, groupsSize);
    }
    
  
}