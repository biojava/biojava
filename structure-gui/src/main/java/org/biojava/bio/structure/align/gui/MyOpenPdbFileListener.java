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
					
					jmol.evalString("select * ; color chain;");
					jmol.evalString("select *; spacefill off; wireframe off; backbone 0.5;  model 0;  ");

				} catch (Exception ex){
					ex.printStackTrace();
				}


			}
		}				
	}		
}

