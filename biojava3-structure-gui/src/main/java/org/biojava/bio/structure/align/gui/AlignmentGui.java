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
package org.biojava.bio.structure.align.gui;


import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.util.logging.Logger;
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

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.ProgressThreadDrawer;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.AligUIManager;

import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.gui.util.PDBUploadPanel;
import org.biojava.bio.structure.gui.util.ScopSelectPanel;
import org.biojava.bio.structure.gui.util.StructurePairSelector;


/** A JFrame that allows to trigger a pairwise structure alignment,
 * either from files in a directory,
 * or after manual upload.
 *
 * @author Andreas Prlic
 *
 * @since 1.7
 *
 *
 *
 */
public class AlignmentGui extends JFrame{

	private final static long serialVersionUID =0l;

	public static Logger logger =  Logger.getLogger("org.biojava.spice");

	StructureAlignment algorithm;

	JButton abortB;

	SelectPDBPanel  tab1 ;
	PDBUploadPanel  tab2;
	ScopSelectPanel tab3;

	Thread thread;
	AlignmentCalculationRunnable alicalc;
	JTabbedPane masterPane;
	JTabbedPane tabPane;
	JProgressBar progress;


	private DBSearchGUI dbsearch;


	public static void main(String[] args){
		
		AlignmentGui.getInstance();

	}

	static final ResourceManager resourceManager = ResourceManager.getResourceManager("ce");

	private static final String MAIN_TITLE = "Pairwise Structure Alignment - Main - V." + resourceManager.getString("ce.version");;

	private static final AlignmentGui me = new AlignmentGui();

	public static AlignmentGui getInstance(){
				
		AbstractUserArgumentProcessor.printAboutMe();
		
		AligUIManager.setLookAndFeel();

		if (!  me.isVisible())
			me.setVisible(true);

		if ( ! me.isActive())
			me.requestFocus();


		return me;
	}

	public static AlignmentGui getInstanceNoVisibilityChange(){

		return me;
	}


	private AlignmentGui() {
		super();

		thread = null;

		JMenuBar menu = MenuCreator.initAlignmentGUIMenu(this);

		this.setJMenuBar(menu);

		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		this.setTitle(MAIN_TITLE);

		tab1 = new SelectPDBPanel();
		tab2 = new PDBUploadPanel();
		tab3 = new ScopSelectPanel();

		// setup tabPane
		tabPane = new JTabbedPane();

		tabPane.addTab("Select PDB ID", null, tab1, "Select PDB ID to align");

		tabPane.addTab("Domains",null, tab3,"Select domains to align.");
		
		tabPane.addTab("Custom files",null, tab2,"Align your own files.");

		

		Box hBoxAlgo = setupAlgorithm();

		Box vBox = Box.createVerticalBox();


		//vBox.add(hBoxAlgo);

		vBox.add(tabPane);
		vBox.add(Box.createGlue());

		//vBox.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));

		masterPane = new JTabbedPane();

		masterPane.addTab("Pairwise Comparison", vBox);

		dbsearch = new DBSearchGUI();

		masterPane.addTab("Database Search",dbsearch);

		//JPanel dir = tab1.getPDBDirPanel(pdbDir);

		Box vBoxMain = Box.createVerticalBox();
		vBoxMain.add(hBoxAlgo);

		// pairwise or db search

		vBoxMain.add(masterPane);

		// algorithm selection

		// PDB install config
		//vBoxMain.add(dir);
		// buttons
		vBoxMain.add(initButtons());

		this.getContentPane().add(vBoxMain);

		//SwingUtilities.updateComponentTreeUI( me);

		this.pack();
		this.setVisible(true);


	}

	private Box setupAlgorithm()
	{
		String[] algorithms = StructureAlignmentFactory.getAllAlgorithmNames();
		try {
			algorithm = StructureAlignmentFactory.getAlgorithm(algorithms[0]);
		} catch (StructureException e){
			e.printStackTrace();
		}

		JLabel algoLabel = new JLabel("Select alignment algorithm: ");

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


		algorithmList.addActionListener(actionAlgorithm);


		Action paramAction = new AbstractAction("Parameters") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				//System.out.println("calc structure alignment");
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

		//        Box hBox42 = Box.createHorizontalBox();
		progress =new JProgressBar();      
		progress.setIndeterminate(false);
		progress.setMaximumSize(new Dimension(10,100));
		progress.setVisible(false);

		//        hBox42.add(Box.createGlue());
		//        hBox42.add(progress);       
		//        hBox42.add(Box.createGlue());
		//        vBox.add(hBox42);
		Action action1 = new AbstractAction("Align") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				//System.out.println("calc structure alignment");
				int selectedIndex = masterPane.getSelectedIndex();
				if (selectedIndex == 0)
					calcAlignment();
				else if ( selectedIndex == 1)
					calcDBSearch();
				else {
					System.err.println("Unknown TAB: " + selectedIndex);
				}

			}
		};

		JButton submitB = new JButton(action1);

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
		StructureAlignment algorithm = getStructureAlignment();
		System.out.println("configure parameters for " + algorithm.getAlgorithmName());

		// show a new config GUI
		new ParameterGUI(algorithm);

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
			Structure s1 = tab.getStructure1();
			Structure s2 = tab.getStructure2();

			if ( s1 == null) {
				System.err.println("please select structure 1");
				return ;
			}

			if ( s2 == null) {
				System.err.println("please select structure 2");
				return;
			}

			String name1 = "custom1"; 
			String name2 = "custom2";

			if  ( pos == 0){
				name1 = tab1.getName1();
				name2 = tab1.getName2();
			} else {
				name1 = s1.getName();
				name2 = s2.getName();
			}
			
			System.out.println("aligning: " + name1 + " " + name2);


			alicalc = new AlignmentCalc(this,s1,s2, name1, name2);


			thread = new Thread(alicalc);
			thread.start();
			abortB.setEnabled(true);
			progress.setIndeterminate(true);
			ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
			drawer.start();
		} catch (StructureException e){
			JOptionPane.showMessageDialog(null,"Could not align structures. Exception: " + e.getMessage());
		}

	}


	private void calcDBSearch() {

		JTabbedPane tabPane = dbsearch.getTabPane();
		System.out.println("run DB search " + tabPane.getSelectedIndex());

		Structure s = null;
		boolean domainSplit = dbsearch.isDomainSplit();
		
		StructurePairSelector tab = null;
		int pos = tabPane.getSelectedIndex();

		if (pos == 0 ){

			tab = dbsearch.getSelectPDBPanel();

		}  else if (pos == 1){

			tab = dbsearch.getScopSelectPanel();


		} else if (pos == 2){

			tab = dbsearch.getPDBUploadPanel();

		}

		try {

			s = tab.getStructure1();

			if ( s == null) {
				System.err.println("please select structure 1");
				return ;
			}

		} catch (Exception e){
			e.printStackTrace();
		}

		String name1 = s.getName();
		if ( name1 == null || name1.equals(""))
			name1 = s.getPDBCode();

		
		
		System.out.println("name1 in alig gui:" + name1);
		String file = dbsearch.getOutFileLocation();
		if ( file == null || file.equals("")){
			JOptionPane.showMessageDialog(null,"Please select a directory to contain the DB search results.");
			return;
		}

		UserConfiguration config = WebStartMain.getWebStartConfig();

		int totalNrCPUs = Runtime.getRuntime().availableProcessors();

		int useNrCPUs = 1;
		if ( totalNrCPUs > 1){
			Object[] options = new Integer[totalNrCPUs];
			int posX = 0;
			for ( int i = totalNrCPUs; i> 0 ; i--){
				options[posX] = i;
				posX++;
			}
			int n = JOptionPane.showOptionDialog(null,
					"How many would you like to use for the calculations?",
					"We detected " + totalNrCPUs + " processors on your system.",    			  
					JOptionPane.OK_CANCEL_OPTION,
					JOptionPane.QUESTION_MESSAGE,
					null,
					options,
					options[0]);

			if ( n < 0)
				return;
			useNrCPUs = (Integer) options[n];
			System.out.println("will use " + useNrCPUs + " CPUs." );
		}
		System.out.println("using domainSplit data");
		alicalc = new AlignmentCalcDB(this, s,  name1,config,file, domainSplit);
		alicalc.setNrCPUs(useNrCPUs);
		abortB.setEnabled(true);
		progress.setIndeterminate(true);
		ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
		drawer.start();

		Thread t = new Thread(alicalc);
		t.start();
	}


	public DBSearchGUI getDBSearch(){
		return dbsearch;
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


	public StructureAlignment getStructureAlignment() {

		return algorithm;
	}

	private void updateAlgorithm(String algorithmName) {

		//String algorithmName = (String)algorithmList.getSelectedItem();
		try {
			algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
		} catch (StructureException ex){
			ex.printStackTrace();
		}

	}

}

class ProgressThreadDrawer extends Thread {

	JProgressBar progress;
	static int interval = 300;

	public ProgressThreadDrawer(JProgressBar progress) {
		this.progress = progress;
	}


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
