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
 * Created on Nov 3, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;


import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;


import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;

import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;


import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.gui.util.PDBUploadPanel;

public class DBSearchGUI extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5657960663049062301L;
	

	StructureAlignment algorithm;
	SelectPDBPanel tab1;
	JTabbedPane tabPane;

	PDBUploadPanel tab2;
	JTabbedPane listPane;
	JButton abortB;
	AlignmentCalcDB alicalc;
	JProgressBar progress;
	ProgressThreadDrawer drawer;
	JTextField outFileLocation;
	
	static final ResourceManager resourceManager = ResourceManager.getResourceManager("ce");
	
	
	public DBSearchGUI(){
	
		
		tab1 = new SelectPDBPanel(false);
		
		tab2 = new PDBUploadPanel(false);
		
		tabPane = new JTabbedPane();
		tabPane.addTab("Select PDB ID", null, tab1,"Select PDB ID to align");

		
		tabPane.addTab("Custom files",null, tab2,"Align your own files.");
		
		listPane = createListPane();
		
		// build up UO

		Box vBox = Box.createVerticalBox();
				
		vBox.add(tabPane);
	
		vBox.add(listPane);
					
		//vBox.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));
		vBox.add(Box.createGlue());
				
		this.add(vBox);
		
		this.setVisible(true);

	}
	
	public JTabbedPane getTabPane()
   {
      return tabPane;
   }

   public void setTabPane(JTabbedPane tabPane)
   {
      this.tabPane = tabPane;
   }



   public SelectPDBPanel getSelectPDBPanel(){
	   return tab1;
	}
	public PDBUploadPanel getPDBUploadPanel(){
	   return tab2;
	}
	public String getOutFileLocation(){
	   return outFileLocation.getText();
	}


	private JTabbedPane createListPane() {
		JTabbedPane tabP = new JTabbedPane();
		
		JPanel dir = new JPanel();

		
		outFileLocation = new JTextField(20);
		JButton chB = new JButton("Select");
		
		Box hBox = Box.createHorizontalBox();
		hBox.add(outFileLocation);
		hBox.add(Box.createGlue());
		hBox.add(chB);
		
		dir.add(hBox);
		
		chB.addActionListener(new ActionListener() {
			
			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setDialogTitle("Select Output Directory");
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				//
				// disable the "All files" option.
				//
				chooser.setAcceptAllFileFilterUsed(false);
				//    


//				In response to a button click:
				int returnVal = chooser.showSaveDialog(null);
				if ( returnVal == JFileChooser.APPROVE_OPTION) {
					File file = chooser.getSelectedFile();
					outFileLocation.setText(file.getAbsolutePath());
					outFileLocation.repaint();
				}
				
			}
		});
		
		tabP.addTab("Select Output Directory", null, dir,
		"Configure the folder that will contain the results.");
		
		return tabP;
	}


	
	public void notifyCalcFinished(){
		if ( drawer != null)
			drawer.interrupt();
		abortB.setEnabled(false);		
		progress.setIndeterminate(false);
		
	}

	

	public StructureAlignment getStructureAlignment() {

		return algorithm;
	}
}

