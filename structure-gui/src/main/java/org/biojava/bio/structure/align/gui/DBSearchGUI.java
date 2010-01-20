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

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.gui.util.PDBDirPanel;
import org.biojava.bio.structure.gui.util.PDBUploadPanel;

public class DBSearchGUI extends JFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5657960663049062301L;
	

	StructureAlignment algorithm;
	SelectPDBPanel tab1;
	JTabbedPane tabPane;
	JTextField pdbDir;
	PDBUploadPanel tab2;
	JTabbedPane listPane;
	JButton abortB;
	AlignmentCalcDB alicalc;
	JProgressBar progress;
	ProgressThreadDrawer drawer;
	JTextField outFileLocation;
	
	static final ResourceManager resourceManager = ResourceManager.getResourceManager("ce");
	private static final String DB_SEARCH_TITLE = "DB Search Structure Alignment - Main - V." + resourceManager.getString("ce.version");
	
	private static final  DBSearchGUI me = new DBSearchGUI();
	
	public static void main(String[] args){
		
		System.setProperty(PDBDirPanel.PDB_DIR,"/Users/andreas/WORK/PDB/");
		DBSearchGUI.getInstance();
		
		
	}
	
	private DBSearchGUI(){
		JMenuBar menu = MenuCreator.initAlignmentGUIMenu(this);

		this.setJMenuBar(menu);

		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setTitle(DB_SEARCH_TITLE);


		Action action2 = new AbstractAction("Exit") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				abortCalc();
				dispose();
				System.exit(0);
			}
		};
		JButton closeB = new JButton(action2);

		Action action1 = new AbstractAction("Align") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				//System.out.println("calc structure alignment");
				calcDBSearch();

			}

		};

		JButton submitB = new JButton(action1);

		JLabel algoLabel = new JLabel("Select alignment algorithm: ");
		String[] algorithms = StructureAlignmentFactory.getAllAlgorithmNames();
		try {
			algorithm = StructureAlignmentFactory.getAlgorithm(algorithms[0]);
		} catch (StructureException e){
			e.printStackTrace();
		}
		JComboBox algorithmList = new JComboBox(algorithms);
		algorithmList.setSelectedIndex(0);

		Action actionAlgorithm = new AbstractAction("Algorithm") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				JComboBox cb = (JComboBox)evt.getSource();
				String algorithmName = (String)cb.getSelectedItem();
				// Perform action...
				//System.out.println("calc structure alignment");
				updateAlgorithm(algorithmName);

			}
		};

		progress =new JProgressBar();
		progress.setIndeterminate(false);
		progress.setMaximumSize(new Dimension(10,100));

		
		Action action3 = new AbstractAction("Abort") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				abortCalc();
			}
		};
		abortB = new JButton(action3);

		abortB.setEnabled(false);

		algorithmList.addActionListener(actionAlgorithm);

		tab1 = new SelectPDBPanel(false);
		pdbDir = tab1.getPDBDirField();
		
		tab2 = new PDBUploadPanel(false);
		
		tabPane = new JTabbedPane();
		tabPane.addTab("Select PDB ID", null, tab1,"Select PDB ID to align");

		
		tabPane.addTab("Custom files",null, tab2,"Align your own files.");
		
	
		JPanel dir = tab1.getPDBDirPanel(pdbDir);

		JTabbedPane configPane = new JTabbedPane();

		configPane.addTab("Local PDB install", null, dir,
		"Configure your local PDB setup.");

	
		listPane = createListPane();
		
		// build up UO


		Box vBox = Box.createVerticalBox();

		Box hBoxAlgo = Box.createHorizontalBox();
		hBoxAlgo.add(Box.createGlue());
		hBoxAlgo.add(algoLabel);		
		hBoxAlgo.add(algorithmList);
		hBoxAlgo.add(Box.createGlue());
		vBox.add(hBoxAlgo);
		
		vBox.add(tabPane);
		vBox.add(configPane);
		vBox.add(listPane);
					
		Box hBox = Box.createHorizontalBox();
	
		hBox.add(closeB);

		
		hBox.add(Box.createGlue());
		hBox.add(progress);
		//hBox.add(Box.createGlue());
		hBox.add(abortB);
		//hBox.add(Box.createGlue());
		hBox.add(submitB);
		
		vBox.add(hBox);

		vBox.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));

		this.getContentPane().add(vBox);
		this.pack();
		this.setVisible(true);

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
				chooser.setDialogTitle("Select output file.");
				chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
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
		
		tabP.addTab("Select output file", null, dir,
		"Configure the file that will contain the results.");
		
		return tabP;
	}


	private void updateAlgorithm(String algorithmName) {

		//String algorithmName = (String)algorithmList.getSelectedItem();
		try {
			algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
		} catch (StructureException ex){
			ex.printStackTrace();
		}

	}

	private void calcDBSearch() {
		// TODO Auto-generated method stub
		System.out.println("run DB search");
		
		String name1 = tab1.getName1();
		tab1.persistCurrentConfig();
		Structure s = null;
		try {
			s = tab1.getStructure1();
		} catch (Exception e){
			e.printStackTrace();
			JOptionPane.showMessageDialog(null,"Could not perform DB search. Exception: " + e.getMessage());
			return;
		}
		if ( s == null) {
			JOptionPane.showMessageDialog(null,"Could not perform DB search. Structure: "  + name1 + " not found!");
			return;
		}
		
		String file = outFileLocation.getText();
		if ( file == null || file.equals("")){
			JOptionPane.showMessageDialog(null,"Plrease select a file to contain the DB search results.");
			return;
		}
		
		UserConfiguration config = tab1.getConfiguration();
		alicalc = new AlignmentCalcDB(this, s,  name1,config,file);
		
		
		
		abortB.setEnabled(true);
		progress.setIndeterminate(true);
		drawer = new ProgressThreadDrawer(progress);
		drawer.start();
		
		Thread t = new Thread(alicalc);
		t.start();
	}

	private void abortCalc(){
		System.err.println("Interrupting alignment ...");
		if ( alicalc != null )
			alicalc.interrupt();
		notifyCalcFinished();
		

	}
	
	public void notifyCalcFinished(){
		if ( drawer != null)
			drawer.interrupt();
		abortB.setEnabled(false);		
		progress.setIndeterminate(false);
		
	}

	public static final DBSearchGUI getInstance(){
				
		if (!  me.isVisible())
			me.setVisible(true);
		
		if ( ! me.isActive())
			me.requestFocus();
		
		return me;
	}


	public StructureAlignment getStructureAlignment() {

		return algorithm;
	}
}

