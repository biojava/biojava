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
 * created at May 24, 2008
 */
package org.biojava.nbio.structure.align.gui;


import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;


/**
 * Create the menus for structure alignment GUI windows (JFrames).
 * <p>
 * Examples: Text Frames, Alignment Panels, Jmol Panels.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 * @author Spencer Bliven
 * @since 1.7
 *
 */
public class MenuCreator {

	public static final String PRINT     = "Print";
	public static final String ALIGNMENT_PANEL = "Alignment Panel";
	public static final String TEXT_ONLY = "View Text Only";
	public static final String PAIRS_ONLY = "View Aligned Pairs";
	public static final String SELECT_EQR = "Select Equivalent Positions";
	public static final String SIMILARITY_COLOR = "Color By Similarity";
	public static final String EQR_COLOR = "Color By EQR";
	public static final String FATCAT_BLOCK = "Color By Alignment Block";
	public static final String LOAD_DB_RESULTS = "Load DB search results";
	public static final String SAVE_ALIGNMENT_XML = "Save Alignment XML";
	public static final String LOAD_ALIGNMENT_XML = "Load Alignment XML";
	public static final String FATCAT_TEXT = "View as FATCAT result";
	public static final String FASTA_FORMAT = "View FASTA Alignment";
	public static final String DIST_MATRICES = "Show Distance Matrices";
	public static final String DOT_PLOT = "Show Dot Plot";
	public static final String PAIRWISE_ALIGN = "New Pairwise Alignment";
	public static final String MULTIPLE_ALIGN = "New Multiple Alignment";
	public static final String PHYLOGENETIC_TREE = "Phylogenetic Tree";
	protected static final int keyMask =
			Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

	/**
	 * Provide a JMenuBar that can be added to a JFrame containing
	 * a JmolPanel. The alignment has to be either an AFPChain or a
	 * MultipleAlignment: set the other parameter to null.<p>
	 * Menus included:
	 * <ul><li>File: open, save, export, import, exit.
	 * <li>Align: new pairwise alignment, new multiple alignment.
	 * <li>View: aligment panel, aligned pairs, text format,
	 * FatCat format, distance matrices, dot plot.
	 * <li>Help
	 * </ul>
	 *
	 * @return a JMenuBar
	 */
	public static JMenuBar initJmolMenu(JFrame frame,
			AbstractAlignmentJmol parent, AFPChain afpChain,
			MultipleAlignment msa) {

		JMenuBar menu = new JMenuBar();

/// FILE MENU
		JMenu file= new JMenu("File");
		file.setMnemonic(KeyEvent.VK_F);
		file.getAccessibleContext().setAccessibleDescription("File Menu");
		//Load
		if (parent != null){
			JMenuItem loadF = getLoadMenuItem();
			loadF.addActionListener(new MyAlignmentLoadListener());
			file.add(loadF);
		}
		//Save
		JMenuItem saveF = getSaveAlignmentMenuItem(afpChain, msa);
		file.add(saveF);
		//Open PDB
		JMenuItem openPDB = getShowPDBMenuItem();
		file.add(openPDB);
		//Open Import
		JMenuItem openI = getOpenPDBMenuItem();
		file.add(openI);
		//Export
		if (parent != null){
			JMenuItem exportI =  getExportPDBMenuItem(parent);
			file.add(exportI);
		}
		//Open DBI
		JMenuItem openDBI = getDBResultMenuItem();
		file.add(openDBI);
		file.addSeparator();
		//Print
		if ( parent != null){
			JMenuItem print = getPrintMenuItem();
			print.addActionListener(parent.getJmolPanel());
			file.add(print);
		}
		file.addSeparator();
		//Close Frame
		JMenuItem closeI = getCloseMenuItem(frame);
		file.add(closeI);
		//Exit
		JMenuItem exitI = getExitMenuItem();
		file.add(exitI);

		menu.add(file);

/// ALIGN MENU
		JMenu align = new JMenu("Align");
		align.setMnemonic(KeyEvent.VK_A);
		//new Pairwise alignment
		JMenuItem pairI = getPairwiseAlignmentMenuItem();
		align.add(pairI);
		//new Multiple alignment
		JMenuItem multI = getMultipleAlignmentMenuItem();
		align.add(multI);

		menu.add(align);

/// VIEW MENU
		JMenu view = new JMenu("View");
		view.getAccessibleContext().setAccessibleDescription("View Menu");
		view.setMnemonic(KeyEvent.VK_V);

		if ( parent != null){
			//Alignment Panel
			JMenuItem apI = MenuCreator.getIcon(parent,ALIGNMENT_PANEL);
			apI.setMnemonic(KeyEvent.VK_M);
			apI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, keyMask));
			view.add(apI);
			//Text Format
			JMenuItem textI = MenuCreator.getIcon(parent,TEXT_ONLY);
			textI.setMnemonic(KeyEvent.VK_T);
			view.add(textI);
			//Alignment Pairs
			JMenuItem pairsI = MenuCreator.getIcon(parent,PAIRS_ONLY);
			pairsI.setMnemonic(KeyEvent.VK_P);
			view.add(pairsI);
			//FatCat Format
			JMenuItem textF = MenuCreator.getIcon(parent,FATCAT_TEXT);
			textF.setMnemonic(KeyEvent.VK_F);
			view.add(textF);
			//Distance Matrices
			JMenuItem distMax = new  JMenuItem(DIST_MATRICES);
			distMax.setMnemonic(KeyEvent.VK_D);
			distMax.addActionListener(new MyDistMaxListener(parent));
			view.add(distMax);
			//Dot Plot - only if the alignment was an afpChain
			if (afpChain != null){
				JMenuItem dotplot = new JMenuItem(DOT_PLOT);
				dotplot.setMnemonic(KeyEvent.VK_O);
				dotplot.addActionListener(new DotPlotListener(afpChain));
				view.add(dotplot);
			}
			//Phylogenetics - only if it is a MultipleAlignment
			if (afpChain == null){
				JMenuItem tree = getIcon(parent, PHYLOGENETIC_TREE);
				tree.setMnemonic(KeyEvent.VK_T);
				view.add(tree);
			}
		}
		menu.add(view);

/// HELP MENU
		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_H);

		JMenuItem helpM = getHelpMenuItem();
		about.add(helpM);

		JMenuItem aboutM = getAboutMenuItem();
		about.add(aboutM);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;
	}


	public static JMenuItem getDBResultMenuItem() {

		ImageIcon saveicon = createImageIcon("/icons/kpdf.png");
		JMenuItem saveI = null;

		if ( saveicon == null)
			saveI = new JMenuItem(LOAD_DB_RESULTS);
		else
			saveI = new JMenuItem(LOAD_DB_RESULTS, saveicon);

		saveI.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				final JFileChooser fc = new JFileChooser();

				//					In response to a button click:
				int returnVal = fc.showOpenDialog(null);
				if ( returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();

					UserConfiguration config = WebStartMain.getWebStartConfig();
					DBResultTable table = new DBResultTable();
					table.show(file,config);
				}

			}
		} );

		return saveI;
	}

	public static JMenuItem getShowPDBMenuItem() {
		ImageIcon loadI = createImageIcon("/icons/background.png");
		JMenuItem openI = null;

		if ( loadI == null)
			openI =new JMenuItem("Show By ID");
		else
			openI =new JMenuItem("Show By ID", loadI);
		openI.setMnemonic(KeyEvent.VK_O);
		openI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, keyMask));
		openI.addActionListener(new ShowPDBIDListener());
		return openI;
	}


	public static JMenuItem getOpenPDBMenuItem() {
		ImageIcon loadI = createImageIcon("/icons/background.png");
		JMenuItem openI = null;

		if ( loadI == null)
			openI =new JMenuItem("Open PDB file");
		else
			openI =new JMenuItem("Open PDB file", loadI);
		openI.setMnemonic(KeyEvent.VK_F);
		openI.addActionListener(new MyOpenPdbFileListener());
		return openI;
	}


	public static JMenuItem getLoadMenuItem() {

		JMenuItem loadF = null;
		ImageIcon loadI = createImageIcon("/icons/revert.png");
		if ( loadI == null)
			loadF = new JMenuItem(LOAD_ALIGNMENT_XML);
		else
			loadF = new JMenuItem(LOAD_ALIGNMENT_XML, loadI);
		loadF.setMnemonic(KeyEvent.VK_L);
		loadF.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_L, keyMask));

		return loadF;
	}


	/**
	 * Create the menu for the Alignment Panel representation of
	 * Structural Alignments. The alignment can be in AFPChain format
	 * or in the MultipleAlignment format.
	 *
	 * @param frame
	 * @param actionListener
	 * @param afpChain
	 * @param MultipleAlignment
	 * @return a JMenuBar
	 */
	public static JMenuBar getAlignmentPanelMenu(JFrame frame,
			ActionListener actionListener, AFPChain afpChain,
			MultipleAlignment msa){

		JMenuBar menu = new JMenuBar();

		JMenu file= new JMenu("File");
		file.getAccessibleContext().setAccessibleDescription("File Menu");
		menu.add(file);

		ImageIcon saveicon = createImageIcon("/icons/filesave.png");

		JMenuItem saveF = null;

		if (saveicon != null)
			saveF = new JMenuItem("Save text display", saveicon);
		else
			saveF = new JMenuItem("Save text display");

		saveF.setMnemonic(KeyEvent.VK_S);
		MySaveFileListener listener = new MySaveFileListener(afpChain, msa);
		listener.setTextOutput(true);
		saveF.addActionListener(listener);
		file.add(saveF);

		file.addSeparator();

		JMenuItem print = getPrintMenuItem();
		print.addActionListener(actionListener);
		file.add(print);

		file.addSeparator();

		JMenuItem closeI = MenuCreator.getCloseMenuItem(frame);
		file.add(closeI);
		JMenuItem exitI = MenuCreator.getExitMenuItem();
		file.add(exitI);

		JMenu edit = new JMenu("Edit");
		edit.setMnemonic(KeyEvent.VK_E);
		menu.add(edit);

		JMenuItem eqrI = MenuCreator.getIcon(actionListener,SELECT_EQR);
		edit.add(eqrI);

		JMenuItem eqrcI = MenuCreator.getIcon(actionListener,EQR_COLOR);
		edit.add(eqrcI);

		JMenuItem simI = MenuCreator.getIcon(actionListener, SIMILARITY_COLOR);
		edit.add(simI);

		JMenuItem fatcatI = MenuCreator.getIcon(actionListener, FATCAT_BLOCK );
		edit.add(fatcatI);

		JMenu view= new JMenu("View");
		view.getAccessibleContext().setAccessibleDescription("View Menu");
		view.setMnemonic(KeyEvent.VK_V);
		menu.add(view);

		JMenuItem textI = MenuCreator.getIcon(actionListener,TEXT_ONLY);
		view.add(textI);

		JMenuItem fastaI = MenuCreator.getIcon(actionListener,FASTA_FORMAT);
		view.add(fastaI);

		JMenuItem pairsI = MenuCreator.getIcon(actionListener,PAIRS_ONLY);
		view.add(pairsI);

		JMenuItem textF = MenuCreator.getIcon(actionListener,FATCAT_TEXT);
		view.add(textF);

		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_A);

		JMenuItem helpM = MenuCreator.getHelpMenuItem();
		about.add(helpM);

		JMenuItem aboutM = MenuCreator.getAboutMenuItem();
		about.add(aboutM);


		menu.add(Box.createGlue());
		menu.add(about);

		return menu;
	}

	/**
	 * Create the menu for the Text representations of Structural Alignments.
	 * @param frame
	 * @param actionListener
	 * @param afpChain
	 * @param msa
	 * @return a JMenuBar
	 */
	public static JMenuBar getAlignmentTextMenu(JFrame frame,
			ActionListener actionListener, AFPChain afpChain,
			MultipleAlignment msa){

		JMenuBar menu = new JMenuBar();

		JMenu file= new JMenu("File");
		file.getAccessibleContext().setAccessibleDescription("File Menu");
		menu.add(file);

		ImageIcon saveicon = createImageIcon("/icons/filesave.png");

		JMenuItem saveF = null;

		if (saveicon != null )
			saveF = new JMenuItem("Save text display", saveicon);
		else
			saveF = new JMenuItem("Save text display");

		saveF.setMnemonic(KeyEvent.VK_S);
		MySaveFileListener listener =
				new MySaveFileListener(afpChain, msa);
		listener.setTextOutput(true);
		saveF.addActionListener(listener);
		file.add(saveF);
		file.addSeparator();

		JMenuItem print = getPrintMenuItem();
		print.addActionListener(actionListener);
		file.add(print);

		file.addSeparator();

		JMenuItem closeI = MenuCreator.getCloseMenuItem(frame);
		file.add(closeI);
		JMenuItem exitI = MenuCreator.getExitMenuItem();
		file.add(exitI);

		JMenu view= new JMenu("View");
		view.getAccessibleContext().setAccessibleDescription("View Menu");
		view.setMnemonic(KeyEvent.VK_V);
		menu.add(view);

		JMenuItem textI = MenuCreator.getIcon(actionListener,TEXT_ONLY);
		view.add(textI);

		JMenuItem fastaI = MenuCreator.getIcon(actionListener,FASTA_FORMAT);
		view.add(fastaI);

		JMenuItem pairsI = MenuCreator.getIcon(actionListener,PAIRS_ONLY);
		view.add(pairsI);

		JMenuItem textF = MenuCreator.getIcon(actionListener,FATCAT_TEXT);
		view.add(textF);

		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_A);

		JMenuItem helpM = MenuCreator.getHelpMenuItem();
		about.add(helpM);

		JMenuItem aboutM = MenuCreator.getAboutMenuItem();
		about.add(aboutM);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;
	}

	protected static JMenuItem getIcon(ActionListener actionListener, String text) {
		JMenuItem to = new JMenuItem(text);

		to.setMnemonic(KeyEvent.VK_E);
		to.addActionListener(actionListener);

		return to;
	}

	public static JMenuItem getPrintMenuItem() {

		//ImageIcon printIcon = createImageIcon("/icons/print_printer.png");
		ImageIcon printIcon = createImageIcon("/icons/fileprint.png");
		JMenuItem print ;

		if ( printIcon == null)
			print = new JMenuItem(PRINT);
		else
			print= new JMenuItem(PRINT,printIcon);


		print.setMnemonic(KeyEvent.VK_P);
		print.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_P, keyMask));

		return print;
	}

	public static JMenuItem getExportPDBMenuItem(AbstractAlignmentJmol parent) {
		ImageIcon saveicon = createImageIcon("/icons/compfile.png");
		JMenuItem exportI = null;

		if ( saveicon == null)
			exportI = new JMenuItem("Export PDB file");
		else
			exportI = new JMenuItem("Export PDB file", saveicon);

		exportI.setMnemonic(KeyEvent.VK_E);
		exportI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, keyMask));
		exportI.addActionListener(new MyExportListener(parent));
		return exportI;
	}

	public static JMenuItem getSaveAlignmentMenuItem(AFPChain afpChain,
			MultipleAlignment msa){

		ImageIcon saveicon = createImageIcon("/icons/filesave.png");
		JMenuItem saveF = null;

		if (saveicon == null)
			saveF = new JMenuItem(SAVE_ALIGNMENT_XML);
		else
			saveF = new JMenuItem(SAVE_ALIGNMENT_XML, saveicon);

		saveF.setMnemonic(KeyEvent.VK_S);
		saveF.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, keyMask));
		saveF.addActionListener(new MySaveFileListener(afpChain, msa));

		return saveF;
	}

	public static JMenuItem getAboutMenuItem() {

		ImageIcon helpIcon = createImageIcon("/icons/help.png");

		JMenuItem aboutM = null;

		if ( helpIcon == null)
			aboutM = new  JMenuItem("About this Software");
		else
			aboutM = new  JMenuItem("About this Software", helpIcon);

		aboutM.setMnemonic(KeyEvent.VK_A);
		aboutM.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				MenuCreator.showAboutDialog();

			}
		});
		return aboutM;
	}


	private static JMenuItem getSystemInfoItem()
	{

		ImageIcon helpIcon = createImageIcon("/icons/help.png");

		JMenuItem aboutM = null;

		if ( helpIcon == null)
			aboutM = new  JMenuItem("System Info");
		else
			aboutM = new  JMenuItem("System Info", helpIcon);

		aboutM.setMnemonic(KeyEvent.VK_S);
		aboutM.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				MenuCreator.showSystemInfo();

			}
		});
		return aboutM;
	}

	public static JMenuItem getExitMenuItem(){

		ImageIcon exitIcon = createImageIcon("/icons/exit.png");


		JMenuItem exitI;

		if ( exitIcon == null)
			exitI = new JMenuItem("Quit");
		else
			exitI = new JMenuItem("Quit",exitIcon);
		exitI.setMnemonic(KeyEvent.VK_Q);
		exitI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, keyMask));
		exitI.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("Quit")){
					System.exit(0);
				}
			}
		});
		return exitI;
	}





	public static JMenuItem getHelpMenuItem(){
		ImageIcon helpIcon = createImageIcon("/icons/help.png");
		JMenuItem helpM ;

		if ( helpIcon == null )
			helpM = new JMenuItem("Help");
		else
			helpM = new JMenuItem("Help", helpIcon);

		helpM.setMnemonic(KeyEvent.VK_H);
		helpM.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				HelpDialog d = new HelpDialog();
				d.showDialog();

			}
		});

		return helpM;
	}

	public static JMenuItem getCloseMenuItem(JFrame frame){

		ImageIcon closeIcon = createImageIcon("/icons/editdelete.png");

		JMenuItem closeI ;

		if ( closeIcon == null )
			closeI = new JMenuItem("Close Frame");
		else
			closeI = new JMenuItem("Close Frame", closeIcon);

		closeI.setMnemonic(KeyEvent.VK_C);
		closeI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_W, keyMask));

		class MyCloseListener implements ActionListener{
			JFrame f;
			MyCloseListener(JFrame frame){
				f=frame;
			}
			@Override
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("Close Frame")){
					f.dispose();
				}
			}

		}

		closeI.addActionListener(new MyCloseListener(frame));
		return closeI;
	}


	/**
	 * Provide a display for the pairwise structure alignment.
	 */
	private static void showPairDialog(){
		AlignmentGui gui =  AlignmentGui.getInstance();
		gui.setVisible(true);
	}

	/**
	 * Provide a display for the multiple structure alignment.
	 */
	private static void showMultipleDialog(){
		MultipleAlignmentGUI gui =  MultipleAlignmentGUI.getInstance();
		gui.setVisible(true);
	}

	/**
	 * Show some info about this GUI
	 *
	 */
	public static void showAboutDialog(){

		AboutDialog dialog = new AboutDialog();
		dialog.showDialog();
	}

	public static void showSystemInfo(){
		SystemInfo dialog = new SystemInfo();
		dialog.showDialog();
	}

	/**
	 * Returns an ImageIcon, or null if the path was invalid.
	 * @param path the path to the icon
	 * @return ImageIcon object
	 */
	public static ImageIcon createImageIcon(String path) {

		java.net.URL imgURL = MenuCreator.class.getResource(path);

		if (imgURL != null) {
			return new ImageIcon(imgURL);
		} else {
			System.err.println("Couldn't find file: " + path);
			return null;
		}
	}

	protected static JMenuItem getPairwiseAlignmentMenuItem() {
		ImageIcon alignIcon = createImageIcon("/icons/window_new.png");

		JMenuItem pairI ;
		if ( alignIcon == null)
			pairI = new JMenuItem(PAIRWISE_ALIGN);
		else
			pairI = new JMenuItem(PAIRWISE_ALIGN, alignIcon);
		pairI.setMnemonic(KeyEvent.VK_N);
		pairI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, keyMask));
		pairI.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals(PAIRWISE_ALIGN)){
					MenuCreator.showPairDialog();
				}
			}
		});
		return pairI;
	}

	protected static JMenuItem getMultipleAlignmentMenuItem() {
		ImageIcon alignIcon = createImageIcon("/icons/window_new.png");

		JMenuItem multipleI ;
		if ( alignIcon == null)
			multipleI = new JMenuItem(MULTIPLE_ALIGN);
		else
			multipleI = new JMenuItem(MULTIPLE_ALIGN, alignIcon);
		multipleI.setMnemonic(KeyEvent.VK_N);
		multipleI.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, keyMask));
		multipleI.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals(MULTIPLE_ALIGN)){
					MenuCreator.showMultipleDialog();
				}
			}
		});
		return multipleI;
	}

	/**
	 * Creates a frame to display a DotPlotPanel.
	 * Used by the 'View>Show Dot Plot' menu item
	 *
	 */
	public static class DotPlotListener implements ActionListener {
		private final AFPChain afpChain;
		public DotPlotListener(AFPChain afpChain) {
			this.afpChain = afpChain;
		}
		@Override
		public void actionPerformed(ActionEvent e) {
			String title = String.format("%s vs. %s",
					afpChain.getName1(),afpChain.getName2());

			// Create window
			JFrame frame = new JFrame(title);
			frame.addWindowListener(new WindowAdapter(){
				@Override
				public void windowClosing(WindowEvent e){
					JFrame f = (JFrame) e.getSource();
					f.setVisible(false);
					f.dispose();
				}
			});

			DotPlotPanel dotplot = new DotPlotPanel(afpChain);

			frame.getContentPane().add(dotplot);
			frame.pack();
			frame.setVisible(true);
		}
	}

	public static JMenuBar initAlignmentGUIMenu(JFrame frame) {

		JMenu file= new JMenu("File");
		file.getAccessibleContext().setAccessibleDescription("File Menu");

		JMenuItem loadF = MenuCreator.getLoadMenuItem();
		loadF.addActionListener(new MyAlignmentLoadListener());
		file.add(loadF);

		JMenuItem openPDB = MenuCreator.getShowPDBMenuItem();
		file.add(openPDB);

		JMenuItem openI = MenuCreator.getOpenPDBMenuItem();
		file.add(openI);

		JMenuItem dbI = MenuCreator.getDBResultMenuItem();
		file.add(dbI);
		file.addSeparator();

		JMenuItem configI = MenuCreator.getConfigMenuItem();
		file.add(configI);
		file.addSeparator();

		JMenuItem closeI = MenuCreator.getCloseMenuItem(frame);
		file.add(closeI);
		JMenuItem exitI = MenuCreator.getExitMenuItem();
		file.add(exitI);

		JMenuBar menu = new JMenuBar();
		menu.add(file);

		//JMenu alig = new JMenu("Align");
		//menu.add(alig);

		//JMenuItem pw = MenuCreator.getPairwiseAlignmentMenuItem();
		//alig.add(pw);


		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_A);

		JMenuItem aboutM = MenuCreator.getAboutMenuItem();
		about.add(aboutM);


		JMenuItem techM = MenuCreator.getSystemInfoItem();
		about.add(techM);

		JMenuItem memM = MenuCreator.getMemoryMonitorItem();
		about.add(memM);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;

	}

	private static JMenuItem getMemoryMonitorItem() {
		ImageIcon helpIcon = createImageIcon("/icons/help.png");

		JMenuItem aboutM = null;

		if ( helpIcon == null)
			aboutM = new  JMenuItem("Memory Monitor");
		else
			aboutM = new  JMenuItem("Memory Monitor", helpIcon);

		aboutM.setMnemonic(KeyEvent.VK_M);
		aboutM.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				MenuCreator.showMemoryMonitor();

			}
		});
		return aboutM;
	}

	protected static void showMemoryMonitor() {
		final MemoryMonitor demo = new MemoryMonitor();

		final JFrame f = new JFrame("MemoryMonitor");
		f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		WindowListener l = new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {}
			@Override
			public void windowDeiconified(WindowEvent e) { demo.surf.start(); }
			@Override
			public void windowIconified(WindowEvent e) { demo.surf.stop(); }
		};

		f.addWindowListener(l);


		Box vBox = Box.createVerticalBox();

		vBox.add(demo);

		Box b = Box.createHorizontalBox();
		JButton b1 = new JButton("Run Garbage Collector");

		b1.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				System.gc();

			}
		});

		b.add(b1);
		b.add(Box.createGlue());

		vBox.add(b);
		f.getContentPane().add("Center",vBox);
		f.pack();
		f.setSize(new Dimension(200,200));
		f.setVisible(true);
		demo.surf.start();

	}

	private static JMenuItem getConfigMenuItem() {

		ImageIcon configIcon = createImageIcon("/icons/configure.png");

		JMenuItem configI;

		if ( configIcon == null)
			configI = new JMenuItem("Settings");
		else
			configI = new JMenuItem("Settings",configIcon);
		configI.setMnemonic(KeyEvent.VK_S);
		configI.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("Settings")){
					ConfigPDBInstallPanel.showDialog();
				}
			}
		});
		return configI;
	}
}
