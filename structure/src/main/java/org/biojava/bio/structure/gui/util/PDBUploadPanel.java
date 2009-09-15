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
 * created at May 25, 2008
 */
package org.biojava.bio.structure.gui.util;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBFileReader;

/** A JPanel to upload 2 PDB files.
 * 
 * @author Andreas Prlic
 * @since 1.7
 */
public class PDBUploadPanel 
extends JPanel 
implements StructurePairSelector {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	JTextField filePath1;
	JTextField filePath2;

    public static Logger logger =  Logger.getLogger("org.biojava");
	
	public PDBUploadPanel(){

		Box vBox = Box.createVerticalBox();

		filePath1 = new JTextField(20);
		filePath2 = new JTextField(20);
		JPanel p1 = getLocalFilePanel(1,filePath1);
		JPanel p2 = getLocalFilePanel(2,filePath2);
		vBox.add(p1);
		vBox.add(p2);

		this.add(vBox);
	}

	
	
	public Structure getStructure1() throws StructureException{
		
		return getStructure(filePath1);
	}
	
public Structure getStructure2() throws StructureException{
		
		return getStructure(filePath2);
	}
	
	private Structure getStructure(JTextField filePath) throws StructureException{
		PDBFileReader reader = new PDBFileReader();
		String path = filePath.getText();
		File f = new File(path);
		Structure s = null;
		try {
			s = reader.getStructure(f);
		} catch (IOException  e){
			logger.warning(e.getMessage());
			//e.printStackTrace();
			throw new StructureException(e);
		}
		return s;

	}


	private JPanel getLocalFilePanel(int pos ,JTextField filePath){

		 JPanel panel = new JPanel();
		 panel.setBorder(BorderFactory.createLineBorder(Color.black));
	        
		
		JLabel l01 = new JLabel("File "+pos+":");
		panel.add(l01);
	
		panel.add(filePath);
		Action action3 = new ChooseAction(filePath);
		JButton chooser = new JButton(action3);
		panel.add(chooser);
		return panel;

	}
}

class ChooseAction extends AbstractAction{
	
	JTextField textField;
	public ChooseAction (JTextField textField){
		super("Choose");
		this.textField = textField;
	}
	public static final long serialVersionUID = 0l;
	// This method is called when the button is pressed
	public void actionPerformed(ActionEvent evt) {
		// Perform action...
		final JFileChooser fc = new JFileChooser();

//		In response to a button click:
		int returnVal = fc.showOpenDialog(null);
		if ( returnVal == JFileChooser.APPROVE_OPTION) {
			File file = fc.getSelectedFile();
			textField.setText(file.getAbsolutePath());
			textField.repaint();
			

		}

	}
}
