package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JFileChooser;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.io.PDBFileReader;

public class MyOpenPdbFileListener 
implements ActionListener {
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

