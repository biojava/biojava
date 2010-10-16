///*
// *                  BioJava development code
// *
// * This code may be freely distributed and modified under the
// * terms of the GNU Lesser General Public Licence.  This should
// * be distributed with the code.  If you do not have a copy,
// * see:
// *
// *      http://www.gnu.org/copyleft/lesser.html
// *
// * Copyright for this code is held jointly by the individual
// * authors.  These should be listed in @author doc comments.
// *
// * For more information on the BioJava project and its aims,
// * or to join the biojava-l mailing list, visit the home page
// * at:
// *
// *      http://www.biojava.org/
// *
// * Created on Jan 19, 2010
// *
// */
//package org.biojava.structure.gui;
//
//import java.awt.Color;
//import java.io.File;
//import junit.framework.TestCase;
//import org.biojava.bio.structure.Structure;
//import org.biojava.bio.structure.io.PDBFileParser;
//import org.biojava.bio.structure.io.PDBFileReader;
//import org.biojava3.structure.gui.OpenAstexViewerImpl;
//import org.biojava3.structure.gui.RenderStyle;
//import org.biojava3.structure.gui.Selection;
//import org.biojava3.structure.gui.StructureViewer;
//
//
///**
// *
// * @author Jules Jacobsen
// */
//public class OpenAstexViewerImplTest extends TestCase {
//
//    Structure structure;
//    StructureViewer instance;
//
//    public OpenAstexViewerImplTest(String testName) {
//        super(testName);
//    }
//
//    @Override
//    protected void setUp() throws Exception {
//        instance = new OpenAstexViewerImpl();
//        PDBFileReader pdbr = new PDBFileReader();
//        pdbr.setAutoFetch(true);
////        pdbr.setPath("/Users/andreas/WORK/PDB/");
//
//        String pdbCode = "5pti";
//
//        structure = pdbr.getStructureById(pdbCode);
//    }
//
//    @Override
//    protected void tearDown() throws Exception {
//        instance.clear();
//    }
//
//    /**
//     * Test of setStructure method, of class OpenAstexViewerImpl.
//     */
//    public void testSetStructure() {
//        System.out.println("setStructure");
//        instance.setStructure(structure);
//        //assertTrue(instance.);
//    }
//
//    /**
//     * Test of repaint method, of class OpenAstexViewerImpl.
//     */
//    public void testRepaint() {
//        System.out.println("repaint");
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        instance.repaint();
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setSelection method, of class OpenAstexViewerImpl.
//     */
//    public void testSetSelection() {
//        System.out.println("setSelection");
//        Selection selection = null;
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        instance.setSelection(selection);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSelection method, of class OpenAstexViewerImpl.
//     */
//    public void testGetSelection() {
//        System.out.println("getSelection");
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        Selection expResult = null;
//        Selection result = instance.getSelection();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setColor method, of class OpenAstexViewerImpl.
//     */
//    public void testSetColor() {
//        System.out.println("setColor");
//        Color red = null;
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        instance.setColor(red);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getColor method, of class OpenAstexViewerImpl.
//     */
//    public void testGetColor() {
//        System.out.println("getColor");
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        Color expResult = null;
//        Color result = instance.getColor();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setStyle method, of class OpenAstexViewerImpl.
//     */
//    public void testSetStyle() {
//        System.out.println("setStyle");
//        RenderStyle wireframe = null;
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        instance.setStyle(wireframe);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of clear method, of class OpenAstexViewerImpl.
//     */
//    public void testClear() {
//        System.out.println("clear");
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        instance.clear();
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setZoom method, of class OpenAstexViewerImpl.
//     */
//    public void testSetZoom() {
//        System.out.println("setZoom");
//        int i = 0;
//        OpenAstexViewerImpl instance = new OpenAstexViewerImpl();
//        instance.setZoom(i);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//}
