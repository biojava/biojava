package org.biojava.nbio.structure.gui;
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
//package org.biojava.nbio.structure.gui;
//
//import astex.MoleculeRenderer;
//import astex.MoleculeViewer;
//import astex.splitter.SplitterBar;
//import astex.splitter.SplitterLayout;
//import java.awt.BorderLayout;
//import java.awt.Color;
//import java.awt.FlowLayout;
//import java.awt.Frame;
//import java.awt.Label;
//import java.awt.Panel;
//import java.awt.TextField;
//import org.biojava.nbio.structure.Structure;
//
///**
// *
// * @author Jules
// */
//public class OpenAstexViewerImpl implements StructureViewer {
//
//    private MoleculeViewer moleculeViewer;
//    private MoleculeRenderer renderer;
//
//    @SuppressWarnings("static-access")
//    public OpenAstexViewerImpl() {
//
//        Frame frame = new Frame();
//        frame.setLayout(new BorderLayout());
//
//        moleculeViewer = new MoleculeViewer();
//        frame.setLayout(new BorderLayout());
//
//        moleculeViewer = new MoleculeViewer();
//
//        moleculeViewer.setApplication(true);
//
//        frame.addWindowListener(moleculeViewer);
//
//        if (false) {
//            moleculeViewer.setUsePopupMenu(true);
//            moleculeViewer.createMenuBar();
//        } else {
//            frame.setMenuBar(moleculeViewer.createMenuBar());
//        }
//
//
//        frame.setLayout(new SplitterLayout(SplitterLayout.HORIZONTAL));
//
//        frame.add("3", moleculeViewer);
//        SplitterBar sb = new SplitterBar();
//
//        frame.add(sb);
//
//        frame.add("1", new Panel());
//
//        frame.pack();
//        frame.setVisible(true);
//
//    }
//
//    public void setStructure(Structure structure) {
//        moleculeViewer.loadMolecule(structure.toPDB());
//        moleculeViewer.repaint();
//
//    }
//
//    public void repaint() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public void setSelection(Selection selection) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public Selection getSelection() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public void setColor(Color red) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public Color getColor() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public void setStyle(RenderStyle wireframe) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public void clear() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    public void setZoom(int i) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//}
