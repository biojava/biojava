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

import javax.swing.Box;
import javax.swing.JButton;

import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;


import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.gui.util.PDBUploadPanel;
import org.biojava.bio.structure.gui.util.ScopSelectPanel;

public class DBSearchGUI extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5657960663049062301L;


	StructureAlignment algorithm;
	SelectPDBPanel tab1;
	JTabbedPane tabPane;

	PDBUploadPanel tab2;
	ScopSelectPanel tab3;

	JPanel listPane;
	JButton abortB;
	AlignmentCalcDB alicalc;
	JProgressBar progress;
	ProgressThreadDrawer drawer;
	JTextField outFileLocation;
	
	Boolean useDomainSplit = true;
	static final ResourceManager resourceManager = ResourceManager.getResourceManager("ce");


	public DBSearchGUI(){


		tab1 = new SelectPDBPanel(false);

		tab2 = new PDBUploadPanel(false);
		tab3 = new ScopSelectPanel(false);

		tabPane = new JTabbedPane();
		tabPane.addTab("Select PDB ID", null, tab1,"Select PDB ID to align");

		tabPane.addTab("Domains",null, tab3,"Domains");

		tabPane.addTab("Custom files",null, tab2,"Align your own files.");

		listPane = createListPane();

		// build up UO

		Box vBox = Box.createVerticalBox();

		vBox.add(tabPane);

		vBox.add(listPane);

		//domainSelectPane = createDomainSelectPane();

		//vBox.add(domainSelectPane);

		//vBox.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));
		vBox.add(Box.createGlue());

		this.add(vBox);

		this.setVisible(true);

	}

	public boolean isDomainSplit(){
		return useDomainSplit;
	}
	
	public JTabbedPane getTabPane()
	{
		return tabPane;
	}

	public void setTabPane(JTabbedPane tabPane)
	{
		this.tabPane = tabPane;
	}

	public ScopSelectPanel getScopSelectPanel(){
		return tab3;
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


	private JPanel createListPane() {
		//JTabbedPane tabP = new JTabbedPane();

		
		JLabel lable = new JLabel("Select Output Directory");
		JPanel dir = new JPanel();


		outFileLocation = new JTextField(20);
		JButton chB = new JButton("Select");

		Box fileSelectBox = Box.createHorizontalBox();
		fileSelectBox.add(lable);
		fileSelectBox.add(outFileLocation);
		fileSelectBox.add(chB);
		fileSelectBox.add(Box.createGlue());
		
		
		Box hBox = Box.createVerticalBox();		
		hBox.add(fileSelectBox);

		Box panel =createDomainSelectPane();
		hBox.add(panel);
		
		dir.add(hBox);

		chB.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setMultiSelectionEnabled(false);
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
					outFileLocation.setText(file.getPath());
					outFileLocation.repaint();
				}

			}
		});

		//tabP.addTab("Select Output Directory", null, dir,
		//		"Configure the folder that will contain the results.");


		return dir;
	}


	private Box createDomainSelectPane() {
		

		

		useDomainSplit  = true;

		String[] petStrings = { "Split proteins in Domains", "Use whole chains" };

		//Create the combo box, select item at index 4.
		//Indices start at 0, so 4 specifies the pig.
		JComboBox domainList = new JComboBox(petStrings);
		domainList.setSelectedIndex(0);
		domainList.setToolTipText("Either align whole chains or SCOP domains and domains assigned with PDP, where no SCOP available.");
		domainList.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				JComboBox box = (JComboBox)arg0.getSource();
				int index = box.getSelectedIndex();
				if ( index == 0)
					useDomainSplit = true;
				else 
					useDomainSplit = false;

			}
		});

		JLabel label= new JLabel("Domains:");
		
		Box domainBox = Box.createHorizontalBox();
		domainBox.add(label);
		
		domainBox.add(domainList);
		domainBox.add(Box.createGlue());
		//Box hBox = Box.createHorizontalBox();
		
		//hBox.add(Box.createGlue());

		

		return domainBox;
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

