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


import org.biojava.nbio.structure.Structure;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author Jules
 */
public class ViewerTest {
	StructureViewer viewer;
	Structure structure;

	@Before
	public void setUp(){
//		if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//		//viewer = new OpenAstexViewer();
//		viewer = new JmolViewerImpl();
//		//viewer = new RCSBViewer();
//		structure = new StructureImpl();
	}

	/**
	 * First we want to get a viewer object
	 */


	/**
	 * then load a PDB file.
	 */
	@Test
	public void testStructureLoad(){

//		if (  java.awt.GraphicsEnvironment.isHeadless())
//    		return;
//
//		PDBFileReader parser = new PDBFileReader();
//		parser.setAutoFetch(true);
//		try {
//			structure = parser.getStructureById("4hhb");
//
//			viewer.setStructure(structure);
//
//			// manipulate the coodriantes
//			//
//			//Calc.rotate(structure,Matrix m);
//
//			viewer.repaint();
//
//			Selection selection = new SelectionImpl();
//
//			//selection can be a whole structure, mol_id, chain, residue, atom or SCOP, Pfam, UniProt features
//
//			viewer.setSelection(selection);
//
//			viewer.setColor(Color.RED);
//
//			viewer.setStyle(RenderStyle.WIREFRAME);
//
//			viewer.clear();
//
//			viewer.setZoom(50);
//		} catch (Exception e){
//			fail(e.getMessage());
//		}


	}

}
