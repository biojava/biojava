/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.structure;

import junit.framework.TestCase;
import org.junit.*;

import java.util.HashSet;
import java.util.Set;


/**
 *
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 */
public class ResidueNumberTest extends TestCase {

    public ResidueNumberTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Override
	@Before
    public void setUp() {
    }

    @Override
	@After
    public void tearDown() {
    }

    /**
     * Test of getChainId method, of class ResidueNumber.
     */
    @Test
    public void testGetSetChainId() {
//        System.out.println("getChainId");
        ResidueNumber instance = new ResidueNumber("A", 42, ' ');
        String expResult = "A";
        String result = instance.getChainId();
        assertEquals(expResult, result);
    }

    /**
     * Test of getInsCode method, of class ResidueNumber.
     */
    @Test
    public void testGetSetInsCode() {
//        System.out.println("getInsCode");
        ResidueNumber instance = new ResidueNumber("A", 42, ' ');
        Character expResult = ' ';
        Character result = instance.getInsCode();
        assertEquals(expResult, result);
    }

    /**
     * Test of getSeqNum method, of class ResidueNumber.
     */
    @Test
    public void testGetSetResidueNumber() {
//        System.out.println("getSeqNum");
        ResidueNumber instance = new ResidueNumber("A", 42, ' ');
        Integer expResult = 42;
        Integer result = instance.getSeqNum();
        assertEquals(expResult, result);

    }


    /**
     * Test of equals method, of class ResidueNumber.
     */
    @Test
    public void testEquals() {
//        System.out.println("equals");
        ResidueNumber number1 = new ResidueNumber("A", 42, ' ');
        ResidueNumber number2 = new ResidueNumber("A", 42, ' ');
        boolean expResult = true;
        boolean result = number2.equals(number1);
        assertEquals(expResult, result);

        Set<ResidueNumber> numberSet= new HashSet<ResidueNumber>();
        numberSet.add(number1);
        numberSet.add(number2);
        assertEquals(1, numberSet.size());

    }

    /**
     * Test of hashCode method, of class ResidueNumber.
     */
    @Test
    public void testHashCode() {
//        System.out.println("hashCode");
        ResidueNumber instance = new ResidueNumber("A", 42, ' ');
        int expResult = 93290;
        int result = instance.hashCode();
        assertEquals(expResult, result);
    }

    /**
     * Test of toString method, of class ResidueNumber.
     */
    @Test
    public void testToString() {
//        System.out.println("toString");
        ResidueNumber instance = new ResidueNumber("A", 42, ' ');
        String expResult = "42";
        String result = instance.toString();
        assertEquals(expResult, result);

    }

    /**
     * Test of toPDB method, of class ResidueNumber.
     */
    @Test
    public void testToPDB() {

        ResidueNumber instance  = new ResidueNumber("A", 42, ' ');
        ResidueNumber instance2 = new ResidueNumber("A", 42, null);

        String expResult = "A  42  ";
        String result1 = instance.toPDB();
        assertEquals(expResult, result1);

        String result2 = instance2.toPDB();
        assertEquals(expResult, result2);
    }



}