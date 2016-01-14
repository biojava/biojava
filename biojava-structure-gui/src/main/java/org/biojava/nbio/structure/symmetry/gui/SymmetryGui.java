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
package org.biojava.nbio.structure.symmetry.gui;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.gui.AlignmentCalculationRunnable;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.ParameterGUI;
import org.biojava.nbio.structure.align.gui.SelectPDBPanel;
import org.biojava.nbio.structure.align.util.ResourceManager;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.gui.util.PDBUploadPanel;
import org.biojava.nbio.structure.gui.util.ScopSelectPanel;
import org.biojava.nbio.structure.gui.util.StructurePairSelector;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;

/** 
 * A JFrame that allows to trigger a symmetry analysis, either from files
 * in a directory or after manual upload
 * Adapted from the AlignmentGui class in biojava.
 * Only one structure is inputted and only the CeSymm algorihm can be chosen.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryGui extends JFrame {

	private final static long serialVersionUID = 0l;

	private CeSymm ceSymm = new CeSymm();
	private JButton abortB;

	private SelectPDBPanel  tab1 ;
	private PDBUploadPanel  tab2;
	private ScopSelectPanel tab3;

	private Thread thread;
	private AlignmentCalculationRunnable alicalc;
	private JTabbedPane masterPane;
	private JTabbedPane tabPane;
	private JProgressBar progress;

	public static void main(String[] args){
		SymmetryGui.getInstance();
	}

	static final ResourceManager resourceManager = 
			ResourceManager.getResourceManager("ce");

	private static final String MAIN_TITLE = 
			"Symmetry Analysis Tool: CE-Symm - V.1.0";

	private static final SymmetryGui me = new SymmetryGui();

	public static SymmetryGui getInstance(){

		AbstractUserArgumentProcessor.printAboutMe();

		AligUIManager.setLookAndFeel();

		if (!me.isVisible()) me.setVisible(true);
		if (! me.isActive()) me.requestFocus();

		return me;
	}

	public static SymmetryGui getInstanceNoVisibilityChange(){
		return me;
	}

	private SymmetryGui() {
		super();

		thread = null;

		JMenuBar menu = MenuCreator.initAlignmentGUIMenu(this);
		this.setJMenuBar(menu);

		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setTitle(MAIN_TITLE);

		tab1 = new SelectPDBPanel(false);
		tab2 = new PDBUploadPanel(false);
		tab3 = new ScopSelectPanel(false);

		//setup tabPane
		tabPane = new JTabbedPane();

		tabPane.addTab("Select PDB ID", null, tab1, "Select PDB ID to analyze");
		tabPane.addTab("Domain",null, tab3,"Select domain to analyze.");
		tabPane.addTab("Custom file",null, tab2,"Analyze your own file.");

		Box hBoxAlgo = setupAlgorithm();
		Box vBox = Box.createVerticalBox();

		vBox.add(tabPane);
		vBox.add(Box.createGlue());

		masterPane = new JTabbedPane();
		masterPane.addTab("Symmetry Analysis", vBox);

		Box vBoxMain = Box.createVerticalBox();
		vBoxMain.add(hBoxAlgo);

		vBoxMain.add(masterPane);
		vBoxMain.add(initButtons());

		this.getContentPane().add(vBoxMain);
		this.pack();
		this.setVisible(true);
	}

	private Box setupAlgorithm() {

		String[] algorithms = {"JCE-symmetry"};
		JLabel algoLabel = new JLabel("Symmetry algorithm: ");

		JComboBox algorithmList = new JComboBox(algorithms);
		algorithmList.setSelectedIndex(0);

		Action paramAction = new AbstractAction("Parameters") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			@Override
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				configureParameters();
			}
		};

		JButton parameterButton = new JButton(paramAction);

		Box hBoxAlgo = Box.createHorizontalBox();
		hBoxAlgo.add(Box.createGlue());
		hBoxAlgo.add(algoLabel);      
		hBoxAlgo.add(algorithmList);
		hBoxAlgo.add(Box.createGlue());
		hBoxAlgo.add(parameterButton);
		hBoxAlgo.add(Box.createGlue());
		return hBoxAlgo;
	}

	private Box initButtons(){

		progress =new JProgressBar();
		progress.setIndeterminate(false);
		progress.setMaximumSize(new Dimension(10,100));
		progress.setVisible(false);

		Action action1 = new AbstractAction("Analyze") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			@Override
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				//System.out.println("calc structure alignment");
				int selectedIndex = masterPane.getSelectedIndex();
				if (selectedIndex == 0)
					calcAlignment();
				else {
					System.err.println("Unknown TAB: " + selectedIndex);
				}
			}
		};

		JButton submitB = new JButton(action1);

		Action action3 = new AbstractAction("Abort") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			@Override
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				abortCalc();
			}
		};

		abortB = new JButton(action3);

		abortB.setEnabled(false);

		Action action2 = new AbstractAction("Exit") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			@Override
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				abortCalc();
				dispose();
				System.exit(0);
			}
		};

		JButton closeB = new JButton(action2);
		Box hBox = Box.createHorizontalBox();
		hBox.add(closeB);
		hBox.add(Box.createGlue());
		hBox.add(progress);
		hBox.add(abortB);
		//hBox.add(Box.createGlue());
		hBox.add(submitB);

		return hBox;
	}

	protected void configureParameters() {
		CeSymm algorithm = getSymmetryAlgorithm();
		System.out.println("configure parameters for " + 
				algorithm.getAlgorithmName());

		// show a new config GUI
		ConfigStrucAligParams params = algorithm.getParameters();
		new ParameterGUI(params, algorithm.getAlgorithmName());
	}

	public void cleanUp() {

		if ( alicalc != null) {
			alicalc.cleanup();
		}
	}

	private void calcAlignment() {

		int pos = tabPane.getSelectedIndex();
		StructurePairSelector tab = null;

		if (pos == 0 ){
			tab = tab1;         

		} else if (pos == 1){
			tab = tab3;

		} else if (pos == 2){
			tab = tab2;
		}


		try {
			Structure s = tab.getStructure1();

			if ( s == null) {
				System.err.println("Please select structure");
				return ;
			}

			String name = "custom";

			if  ( pos == 0){
				name = tab1.getName1();
			} else {
				name = s.getName();
			}

			System.out.println("Analyzing: " + name);


			alicalc = new SymmetryCalc(this,s,name);


			thread = new Thread(alicalc);
			thread.start();
			abortB.setEnabled(true);
			progress.setIndeterminate(true);
			ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
			drawer.start();
		} catch (StructureException e){
			JOptionPane.showMessageDialog(null,
					"Could not align structures. Exception: " + e.getMessage());
		}
	}

	public void notifyCalcFinished(){
		abortB.setEnabled(false);
		thread = null;
		progress.setIndeterminate(false);
		this.repaint();
	}

	private void abortCalc(){
		System.err.println("Interrupting alignment ...");
		if ( alicalc != null )
			alicalc.interrupt();
		notifyCalcFinished();
	}


	public CeSymm getSymmetryAlgorithm() {
		return ceSymm;
	}

}

class ProgressThreadDrawer extends Thread {

	JProgressBar progress;
	static int interval = 300;

	public ProgressThreadDrawer(JProgressBar progress) {
		this.progress = progress;
	}


	@Override
	public void run() {
		progress.setVisible(true);
		boolean finished = false;
		while ( ! finished) {
			try {
				progress.repaint();
				if ( ! progress.isIndeterminate() ){
					finished =false;
					break;
				}

				sleep(interval);
			} catch (InterruptedException e){
			}
			progress.repaint();
		}
		progress.setVisible(false);
		progress = null;
	}
}
