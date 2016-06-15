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
package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.io.PDBFileReader;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

public class MyOpenPdbFileListener
implements ActionListener {
	@Override
	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();
		if ( cmd.equals("Open PDB file")){
			final JFileChooser fc = new JFileChooser();

			//					In response to a button click:
			int returnVal = fc.showOpenDialog(null);
			if ( returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();

				PDBFileReader reader = new PDBFileReader();
				try {
					Structure s = reader.getStructure(file);
					StructureAlignmentJmol jmol = new StructureAlignmentJmol(null,null,null);
					jmol.setStructure(s);

					jmol.evalString("set antialiasDisplay on; select all;spacefill off; wireframe off; backbone off; cartoon;color cartoon chain; select ligand;wireframe 0.16;spacefill 0.5; select all; color cartoon structure;");
					jmol.evalString("save STATE state_1");
				} catch (Exception ex){
					ex.printStackTrace();
				}


			}
		}
	}
}

