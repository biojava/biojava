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
 * Created on Sep 28, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.JFileChooser;
import javax.swing.JTextField;

import org.biojava.bio.structure.align.webstart.PersistentConfig;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;


/** Ask the user to provide a directory containting PDB files. 
 * Sets the idr in the provided textField.
 * @author Andreas Prlic
 *
 */
public class ChooseDirAction extends AbstractAction{

	JTextField textField;
	UserConfiguration config;
	public ChooseDirAction (JTextField textField, UserConfiguration config){
		super("Choose");
		this.config = config;
		this.textField = textField;
	}
	public static final long serialVersionUID = 0l;
	// This method is called when the button is pressed
	public void actionPerformed(ActionEvent evt) {
		// Perform action...
		JFileChooser chooser = new JFileChooser();
		String txt = textField.getText();
		
		if ( config == null) {
			System.out.println("config == null, calling getWebStartConfig...");
			config = WebStartMain.getWebStartConfig();
		}
		if ( txt != null){
			chooser.setCurrentDirectory(new java.io.File(txt));
			config.setPdbFilePath(txt);
			try {
				PersistentConfig webstartConfig = new PersistentConfig();

				webstartConfig.save(config);

			} catch (Exception e){
				e.printStackTrace();
			}
		} 
		chooser.setDialogTitle("Choose directory that contains your PDB files");
		chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		//
		// disable the "All files" option.
		//
		chooser.setAcceptAllFileFilterUsed(false);
		//    


//		In response to a button click:
		int returnVal = chooser.showOpenDialog(null);
		if ( returnVal == JFileChooser.APPROVE_OPTION) {
			File file = chooser.getSelectedFile();
			textField.setText(file.getAbsolutePath());
			textField.repaint();
		}

	}
}
