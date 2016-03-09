/*
 *                  BioJava development code
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
 * Created on Jul 16, 2006
 *
 */
package org.biojava.nbio.structure.align.gui;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.MultipleStructureAligner;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcMain;
import org.biojava.nbio.structure.align.webstart.AligUIManager;
import org.biojava.nbio.structure.gui.util.SelectMultiplePanel;

/**
 * A JFrame that allows to trigger a multiple structure alignment,
 * either from files in a directory or after manual upload.
 * <p>
 * The current version allows to select the parameters of
 * the pairwise alignment algorithm and the parameters of
 * the multiple alignment algorithm.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class MultipleAlignmentGUI extends JFrame {

	private final static long serialVersionUID =0l;
	private final static String version = "1.0";

	private MultipleStructureAligner multiple;
	private StructureAlignment pairwise;

	private SelectMultiplePanel tab;
	private JTabbedPane tabPane;

	private Thread thread;
	private AlignmentCalculationRunnable alicalc;
	private JProgressBar progress;
	private JButton abortB;

	private static final String MAIN_TITLE =
			"Multiple Structure Alignment - Main - V." + version;

	private static final MultipleAlignmentGUI me =
			new MultipleAlignmentGUI();

	public static void main(String[] args){
		MultipleAlignmentGUI.getInstance();
	}

	public static MultipleAlignmentGUI getInstance(){

		//TODO change about me
		AbstractUserArgumentProcessor.printAboutMe();
		AligUIManager.setLookAndFeel();

		if (!me.isVisible()) me.setVisible(true);
		if (!me.isActive()) me.requestFocus();

		return me;
	}

	public static MultipleAlignmentGUI getInstanceNoVisibilityChange(){
		return me;
	}

	protected MultipleAlignmentGUI() {
		super();

		thread = null;
		JMenuBar menu = MenuCreator.initAlignmentGUIMenu(this);

		this.setJMenuBar(menu);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setTitle(MAIN_TITLE);

		tab = new SelectMultiplePanel();

		// setup tabPane
		tabPane = new JTabbedPane();
		tabPane.addTab("Select Structures", null,
				tab, "Input Structure identifiers");

		Box hBoxPair = setupPairwiseAlgorithm();
		Box hBoxMult = setupMultipleAlgorithm();
		Box vBox = Box.createVerticalBox();

		vBox.add(tabPane);
		vBox.add(Box.createGlue());

		Box vBoxMain = Box.createVerticalBox();
		vBoxMain.add(hBoxPair);
		vBoxMain.add(hBoxMult);
		vBoxMain.add(tabPane);

		vBoxMain.add(initButtons());
		this.getContentPane().add(vBoxMain);

		this.pack();
		this.setVisible(true);
	}

	private Box setupPairwiseAlgorithm() {

		String[] pairAlgo = StructureAlignmentFactory.getAllAlgorithmNames();
		try {
			pairwise = StructureAlignmentFactory.getAlgorithm(pairAlgo[0]);
		} catch (StructureException e){
			e.printStackTrace();
		}
		JLabel algoLabel = new JLabel("Select pairwise aligner: ");

		JComboBox algorithmList = new JComboBox(pairAlgo);
		algorithmList.setSelectedIndex(0);

		Action actionAlgorithm = new AbstractAction("Algorithm") {
			public static final long serialVersionUID = 0l;
			@Override
			public void actionPerformed(ActionEvent evt) {
				JComboBox cb = (JComboBox)evt.getSource();
				String algorithmName = (String)cb.getSelectedItem();
				updatePairwiseAlgorithm(algorithmName);
			}
		};
		algorithmList.addActionListener(actionAlgorithm);

		Action paramAction = new AbstractAction("Parameters") {
			public static final long serialVersionUID = 0l;
			@Override
			public void actionPerformed(ActionEvent evt) {
				StructureAlignment p = getPairwiseStructureAligner();
				ConfigStrucAligParams params = p.getParameters();
				new ParameterGUI(params, p.getAlgorithmName());
			}
		};
		JButton parameterButton = new JButton(paramAction);

		Box hBoxAlgoPair = Box.createHorizontalBox();
		hBoxAlgoPair.add(Box.createGlue());
		hBoxAlgoPair.add(algoLabel);
		hBoxAlgoPair.add(algorithmList);
		hBoxAlgoPair.add(Box.createGlue());
		hBoxAlgoPair.add(parameterButton);
		hBoxAlgoPair.add(Box.createGlue());

		return hBoxAlgoPair;
	}

	private Box setupMultipleAlgorithm() {

		//TODO change in the future when more multiple algorithms are added
		String[] multAlgo = {MultipleMcMain.algorithmName};
		multiple = new MultipleMcMain(pairwise);

		JLabel multLabel = new JLabel("Select multiple aligner: ");
		JComboBox multList = new JComboBox(multAlgo);
		multList.setSelectedIndex(0);

		Action actionMultiple = new AbstractAction("Algorithm") {
			public static final long serialVersionUID = 0l;
			@Override
			public void actionPerformed(ActionEvent evt) {
				updateMultipleAlgorithm();
			}
		};
		multList.addActionListener(actionMultiple);

		Action paramAction = new AbstractAction("Parameters") {
			public static final long serialVersionUID = 0l;
			@Override
			public void actionPerformed(ActionEvent evt) {
				MultipleStructureAligner m = getMultipleStructureAligner();
				ConfigStrucAligParams params = m.getParameters();
				new ParameterGUI(params, m.getAlgorithmName());
			}
		};
		JButton parameterButton = new JButton(paramAction);

		Box hBoxAlgo = Box.createHorizontalBox();
		hBoxAlgo.add(Box.createGlue());
		hBoxAlgo.add(multLabel);
		hBoxAlgo.add(multList);
		hBoxAlgo.add(Box.createGlue());
		hBoxAlgo.add(parameterButton);
		hBoxAlgo.add(Box.createGlue());

		return hBoxAlgo;
	}

	private Box initButtons(){

		//Progress Bar
		progress = new JProgressBar();
		progress.setIndeterminate(false);
		progress.setMaximumSize(new Dimension(10,100));
		progress.setVisible(false);

		Action action1 = new AbstractAction("Align") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			@Override
			public void actionPerformed(ActionEvent evt) {
				calcAlignment();
			}
		};
		JButton submitB = new JButton(action1);

		Action action3 = new AbstractAction("Abort") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			@Override
			public void actionPerformed(ActionEvent evt) {
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
		hBox.add(submitB);

		return hBox;
	}

	public void cleanUp() {
		if (alicalc != null) alicalc.cleanup();
	}

	private void calcAlignment() {

		try {
			List<Structure> structures = tab.getStructures();

			if ( structures.size() < 2) {
				System.err.println("please input more than 1 structure");
				return;
			}

			List<StructureIdentifier> names = tab.getNames();

			String message = "aligning: ";
			for (StructureIdentifier name:names){
				message += name.getIdentifier() + " ";
			}
			System.out.println(message);

			alicalc = new MultipleAlignmentCalc(this, structures, names);

			thread = new Thread(alicalc);
			thread.start();
			abortB.setEnabled(true);
			progress.setIndeterminate(true);
			ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
			drawer.start();

		} catch (StructureException e){
			JOptionPane.showMessageDialog(null,"Could not align structures. "
					+ "Exception: " + e.getMessage());
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
		if (alicalc != null) alicalc.interrupt();
		notifyCalcFinished();
	}


	public MultipleStructureAligner getMultipleStructureAligner() {
		return multiple;
	}

	public StructureAlignment getPairwiseStructureAligner() {
		return pairwise;
	}

	private void updatePairwiseAlgorithm(String algorithmName) {
		try {
			pairwise = StructureAlignmentFactory.getAlgorithm(algorithmName);
			//Update also the multiple structure algorithm
			ConfigStrucAligParams params = multiple.getParameters();
			updateMultipleAlgorithm();
			multiple.setParameters(params);

		} catch (StructureException ex){
			ex.printStackTrace();
		}
	}

	private void updateMultipleAlgorithm() {
		//TODO a factory would be needed to select the MultipleAligner
		multiple = new MultipleMcMain(pairwise);
	}
}
