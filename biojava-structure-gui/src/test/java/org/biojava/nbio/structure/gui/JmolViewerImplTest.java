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

package org.biojava.nbio.structure.gui;


import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author Jules
 */
public class JmolViewerImplTest {



	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testMe(){
	   Assert.assertTrue(true);
	}

//    /**
//     * Test of setStructure method, of class JmolViewerImpl.
//     */
//    public void testSetStructure() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        try {
//            System.out.println("setStructure");
//            PDBFileReader pdbr = new PDBFileReader();
//            pdbr.setAutoFetch(true);
//            //        pdbr.setPath("/Users/andreas/WORK/PDB/");
//            String pdbCode = "5pti";
//            Structure structure = pdbr.getStructureById(pdbCode);
//            JmolViewerImpl instance = new JmolViewerImpl();
//            instance.setStructure(structure);
//            //try {
//                //Thread.sleep(10000);
//            //} catch (InterruptedException ex) {
//             //   Logger.getLogger(JmolViewerImplTest.class.getName()).log(Level.SEVERE, null, ex);
//            //}
//            // TODO review the generated test code and remove the default call to fail.
//            //fail("The test case is a prototype.");
//        } catch (IOException ex) {
//            Logger.getLogger(JmolViewerImplTest.class.getName()).log(Level.SEVERE, null, ex);
//        }
//    }
//
//    /**
//     * Test of clear method, of class JmolViewerImpl.
//     */
//    public void testClear() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("clear");
//        JmolViewerImpl instance = new JmolViewerImpl();
//        instance.clear();
//        // TODO review the generated test code and remove the default call to fail.
//        //fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getColor method, of class JmolViewerImpl.
//     */
//    public void testGetColor() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("getColor");
//        JmolViewerImpl instance = new JmolViewerImpl();
//        Color expResult = null;
//        Color result = instance.getColor();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        //fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSelection method, of class JmolViewerImpl.
//     */
//    public void testGetSelection() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("getSelection");
//        JmolViewerImpl instance = new JmolViewerImpl();
//        Selection expResult = null;
//        Selection result = instance.getSelection();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        //fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of repaint method, of class JmolViewerImpl.
//     */
//    public void testRepaint() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("repaint");
//        JmolViewerImpl instance = new JmolViewerImpl();
//        instance.repaint();
//        // TODO review the generated test code and remove the default call to fail.
//        //fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setColor method, of class JmolViewerImpl.
//     */
//    public void testSetColor() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("setColor");
//        Color red = null;
//        JmolViewerImpl instance = new JmolViewerImpl();
//        instance.setColor(red);
//        // TODO review the generated test code and remove the default call to fail.
//        //fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setSelection method, of class JmolViewerImpl.
//     */
//    public void testSetSelection() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("setSelection");
//        Selection selection = null;
//        JmolViewerImpl instance = new JmolViewerImpl();
//        instance.setSelection(selection);
//        // TODO review the generated test code and remove the default call to fail.
//      //  fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setStyle method, of class JmolViewerImpl.
//     */
//    public void testSetStyle() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("setStyle");
//        RenderStyle wireframe = null;
//        JmolViewerImpl instance = new JmolViewerImpl();
//        instance.setStyle(wireframe);
//        // TODO review the generated test code and remove the default call to fail.
//        //fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setZoom method, of class JmolViewerImpl.
//     */
//    public void testSetZoom() {
//    	if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//        System.out.println("setZoom");
//        int i = 0;
//        JmolViewerImpl instance = new JmolViewerImpl();
//        instance.setZoom(i);
//        // TODO review the generated test code and remove the default call to fail.
//       // fail("The test case is a prototype.");
//    }

}
