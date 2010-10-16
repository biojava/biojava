package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;

public class MySaveFileListener implements ActionListener {

	AFPChain afpChain ;
	private boolean printFatCat;
	public MySaveFileListener (AFPChain afpChain){
		this.afpChain = afpChain;
		printFatCat = false;
	}

	public void actionPerformed(ActionEvent evt) {

		if ( afpChain == null) {
			JOptionPane.showMessageDialog(null,"Could not save alignment, no alignment being displayed.");
			return;
			
		}
		JFileChooser fc = new JFileChooser();
		int returnVal = fc.showSaveDialog(null);

		if (returnVal != JFileChooser.APPROVE_OPTION) {
			System.err.println("User canceled file save.");
			return;
		}
		File selFile = fc.getSelectedFile();

		if ( selFile == null)
			return;

		System.out.println("Saving alignment to file: " + selFile.getName());
		//if ( ! selFile.canWrite()) {
		//	JOptionPane.showMessageDialog(null,"Don't have permission to write to file " + selFile);
		//}

		// get the XML serialization of the alignment
		try {
			UserConfiguration config = WebStartMain.getWebStartConfig();
			AtomCache cache = new AtomCache(config);
			
			//TODO use the right ca atoms, since this will fail for a custom file!
			Atom[] ca1 =cache.getAtoms(afpChain.getName1());
			Atom[] ca2 =cache.getAtoms(afpChain.getName2());
			
			String output = "";
			if ( ! printFatCat) {			
				output = AFPChainXMLConverter.toXML(afpChain, ca1,ca2);
			} else {
				output = afpChain.toFatcat(ca1, ca2);
			}

			// write to the file
			BufferedWriter out = new BufferedWriter(new FileWriter(selFile));
			out.write(output);
			out.close();


		} catch (Exception e){
			e.printStackTrace();
			JOptionPane.showMessageDialog(null,"Could not save file. Exception: " + e.getMessage());

		}

	}

	public void setFatCatOutput(boolean b) {
		printFatCat = b;

	}

}
