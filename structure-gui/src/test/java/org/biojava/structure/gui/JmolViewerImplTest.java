/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.structure.gui;

import java.awt.Color;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava3.structure.gui.JmolViewerImpl;
import org.biojava3.structure.gui.RenderStyle;
import org.biojava3.structure.gui.Selection;

/**
 *
 * @author Jules
 */
public class JmolViewerImplTest extends TestCase {
    
    public JmolViewerImplTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of setStructure method, of class JmolViewerImpl.
     */
    public void testSetStructure() {
        try {
            System.out.println("setStructure");
            PDBFileReader pdbr = new PDBFileReader();
            pdbr.setAutoFetch(true);
            //        pdbr.setPath("/Users/andreas/WORK/PDB/");
            String pdbCode = "5pti";
            Structure structure = pdbr.getStructureById(pdbCode);
            JmolViewerImpl instance = new JmolViewerImpl();
            instance.setStructure(structure);
            //try {
                //Thread.sleep(10000);
            //} catch (InterruptedException ex) {
             //   Logger.getLogger(JmolViewerImplTest.class.getName()).log(Level.SEVERE, null, ex);
            //}
            // TODO review the generated test code and remove the default call to fail.
            //fail("The test case is a prototype.");
        } catch (IOException ex) {
            Logger.getLogger(JmolViewerImplTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Test of clear method, of class JmolViewerImpl.
     */
    public void testClear() {
        System.out.println("clear");
        JmolViewerImpl instance = new JmolViewerImpl();
        instance.clear();
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of getColor method, of class JmolViewerImpl.
     */
    public void testGetColor() {
        System.out.println("getColor");
        JmolViewerImpl instance = new JmolViewerImpl();
        Color expResult = null;
        Color result = instance.getColor();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of getSelection method, of class JmolViewerImpl.
     */
    public void testGetSelection() {
        System.out.println("getSelection");
        JmolViewerImpl instance = new JmolViewerImpl();
        Selection expResult = null;
        Selection result = instance.getSelection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of repaint method, of class JmolViewerImpl.
     */
    public void testRepaint() {
        System.out.println("repaint");
        JmolViewerImpl instance = new JmolViewerImpl();
        instance.repaint();
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of setColor method, of class JmolViewerImpl.
     */
    public void testSetColor() {
        System.out.println("setColor");
        Color red = null;
        JmolViewerImpl instance = new JmolViewerImpl();
        instance.setColor(red);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of setSelection method, of class JmolViewerImpl.
     */
    public void testSetSelection() {
        System.out.println("setSelection");
        Selection selection = null;
        JmolViewerImpl instance = new JmolViewerImpl();
        instance.setSelection(selection);
        // TODO review the generated test code and remove the default call to fail.
      //  fail("The test case is a prototype.");
    }

    /**
     * Test of setStyle method, of class JmolViewerImpl.
     */
    public void testSetStyle() {
        System.out.println("setStyle");
        RenderStyle wireframe = null;
        JmolViewerImpl instance = new JmolViewerImpl();
        instance.setStyle(wireframe);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of setZoom method, of class JmolViewerImpl.
     */
    public void testSetZoom() {
        System.out.println("setZoom");
        int i = 0;
        JmolViewerImpl instance = new JmolViewerImpl();
        instance.setZoom(i);
        // TODO review the generated test code and remove the default call to fail.
       // fail("The test case is a prototype.");
    }

}
