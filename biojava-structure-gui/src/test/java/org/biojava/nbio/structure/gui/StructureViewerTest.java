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

package org.biojava.nbio.structure.gui;

import junit.framework.TestCase;
import org.biojava.nbio.structure.Structure;

import java.awt.*;

/**
 *
 * @author Jules
 */
public class StructureViewerTest extends TestCase {
    
    public StructureViewerTest(String testName) {
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
     * Test of setStructure method, of class StructureViewer.
     */
    public void testSetStructure() {
        
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        Structure structure = null;
        StructureViewer instance = new StructureViewerImpl();
        instance.setStructure(structure);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of repaint method, of class StructureViewer.
     */
    public void testRepaint() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        
        StructureViewer instance = new StructureViewerImpl();
        instance.repaint();
        // TODO review the generated test code and remove the default call to fail.
       // fail("The test case is a prototype.");
    }

    /**
     * Test of setSelection method, of class StructureViewer.
     */
    public void testSetSelection() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        
        Selection selection = null;
        StructureViewer instance = new StructureViewerImpl();
        instance.setSelection(selection);
        // TODO review the generated test code and remove the default call to fail.
       // fail("The test case is a prototype.");
    }

    /**
     * Test of getSelection method, of class StructureViewer.
     */
    public void testGetSelection() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        
        StructureViewer instance = new StructureViewerImpl();
        Selection expResult = null;
        Selection result = instance.getSelection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
      //  fail("The test case is a prototype.");
    }

    /**
     * Test of setColor method, of class StructureViewer.
     */
    public void testSetColor() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
    	
        
        Color red = null;
        StructureViewer instance = new StructureViewerImpl();
        instance.setColor(red);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of getColor method, of class StructureViewer.
     */
    public void testGetColor() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        
        StructureViewer instance = new StructureViewerImpl();
        Color expResult = null;
        Color result = instance.getColor();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
       // fail("The test case is a prototype.");
    }

    /**
     * Test of setStyle method, of class StructureViewer.
     */
    public void testSetStyle() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        
        RenderStyle wireframe = null;
        StructureViewer instance = new StructureViewerImpl();
        instance.setStyle(wireframe);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of clear method, of class StructureViewer.
     */
    public void testClear() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        StructureViewer instance = new StructureViewerImpl();
        instance.clear();
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }

    /**
     * Test of setZoom method, of class StructureViewer.
     */
    public void testSetZoom() {
    	if (  java.awt.GraphicsEnvironment.isHeadless())
    		return;
        int i = 0;
        StructureViewer instance = new StructureViewerImpl();
        instance.setZoom(i);
        // TODO review the generated test code and remove the default call to fail.
       // fail("The test case is a prototype.");
    }

    public class StructureViewerImpl implements StructureViewer {

        @Override
		public void setStructure(Structure structure) {
        }

        @Override
		public void repaint() {
        }

        @Override
		public void setSelection(Selection selection) {
        }

        @Override
		public Selection getSelection() {
            return null;
        }

        @Override
		public void setColor(Color red) {
        }

        @Override
		public Color getColor() {
            return null;
        }

        @Override
		public void setStyle(RenderStyle wireframe) {
        }

        @Override
		public void clear() {
        }

        @Override
		public void setZoom(int i) {
        }
    }

}
