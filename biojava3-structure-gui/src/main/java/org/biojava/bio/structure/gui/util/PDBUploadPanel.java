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


import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.logging.Logger;

import javax.swing.AbstractAction;
import javax.swing.Action;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.StructureIOFile;

/** A JPanel to upload 2 custom PDB files.
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



	public static Logger logger =  Logger.getLogger("org.biojava");

	
	
	private JComboBox fileType ;
	
	JTextField filePath1;
	JTextField filePath2;
	JTextField chain1;
	JTextField chain2;
	
	public static JComboBox getFileFormatSelect(){
		JComboBox fileType = new JComboBox();
			fileType = new JComboBox(new String[] {UserConfiguration.PDB_FORMAT,UserConfiguration.MMCIF_FORMAT});
			fileType.setSelectedIndex(0);
			fileType.setMaximumSize(new Dimension(10,50));
			
		return fileType;
	}
	public PDBUploadPanel(){
		this(true);
	}
	public PDBUploadPanel(boolean show2boxes){
		Box vBox = Box.createVerticalBox();

		filePath1 = new JTextField(20);
		filePath2 = new JTextField(20);
		chain1 = new JTextField(1);
		chain2 = new JTextField(1);
		
		JPanel p1 = getLocalFilePanel(1,filePath1,chain1);
		JPanel p2 = getLocalFilePanel(2,filePath2,chain2);

		vBox.add(p1);
		if ( show2boxes)
			vBox.add(p2);
		
		JLabel ftype = new JLabel("File format:");
		Box hBox = Box.createHorizontalBox();
		hBox.add(Box.createGlue());
		hBox.add(ftype);
		fileType = getFileFormatSelect();
		hBox.add(fileType);
		hBox.add(Box.createGlue());
		
		vBox.add(hBox);

		this.add(vBox);
	}


	public String getFilePath1(){
		return filePath1.getText();
	}
	
	public String getChain1(){
		return chain1.getText();
	}

	public Structure getStructure1() throws StructureException{

		return getStructure(filePath1,chain1);
	}

	public Structure getStructure2() throws StructureException{

		return getStructure(filePath2,chain2);
	}

	private Structure getStructure(JTextField filePath,JTextField chainId) throws StructureException{
		//PDBFileReader reader = new PDBFileReader();
		
		StructureIOFile reader = null;
		String fileFormat = (String)fileType.getSelectedItem();
		if ( fileFormat.equals(UserConfiguration.PDB_FORMAT)){
			reader = new PDBFileReader();
		} else if ( fileFormat.equals(UserConfiguration.MMCIF_FORMAT)){
			reader = new MMCIFFileReader();
		} else {
			throw new StructureException("Unkown file format " + fileFormat);
		}
		
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
		
		Structure reduced = StructureTools.getReducedStructure(s, chainId.getText());
		
		String fileURL = "";
		try {
						
			URL u ;
			
			if ( chainId.getText() == null || chainId.getText().equals("")){
			
				u = f.toURI().toURL();
			} else {
				u = new URL(f.toURI().toURL().toString() + "?chainId=" + chainId.getText());
			}
			fileURL = u.toString() ; 
					
		} catch (Exception e){
			e.printStackTrace();
		}
		
		reduced.setPDBCode(fileURL);
		reduced.setName(fileURL);
		return reduced;

	}

	

	private JPanel getLocalFilePanel(int pos ,JTextField filePath, JTextField  chainId){

		JPanel panel = new JPanel();
		//panel.setBorder(BorderFactory.createLineBorder(Color.black));

		JLabel l01 = new JLabel("File "+pos+":");
		panel.add(l01);

		panel.add(filePath);
		Action action3 = new ChooseAction(filePath);
		JButton chooser = new JButton(action3);
		panel.add(chooser);
		
		JLabel chainLabel = new JLabel("Chain "+pos+": (optional)");
		panel.add(chainLabel);
		panel.add(chainId);
		
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
