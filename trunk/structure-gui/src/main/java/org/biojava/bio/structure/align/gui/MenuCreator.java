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
package org.biojava.bio.structure.align.gui;


import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;

import javax.swing.Box;
import javax.swing.ImageIcon;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.jama.Matrix;



/**
 *  Create the menu for Fatcat/CE structure alignment stuff
 * @author Andreas Prlic
 * @since 1.7
 */
public class MenuCreator {

	public static final String PRINT     = "Print";
	public static final String ALIGNMENT_PANEL = "Alignment Panel";
	public static final String TEXT_ONLY = "View Text Only";
	public static final String SELECT_EQR = "Select Equivalent Positions";
	public static final String SIMILARITY_COLOR = "Color By Similarity";
	public static final String EQR_COLOR = "Color By EQR";
	public static final String FATCAT_BLOCK = "Color By Alignment Block";
	public static final String LOAD_DB_RESULTS = "Load DB search results";
	public static final String SAVE_ALIGNMENT_XML = "Save Alignment XML";
	public static final String LOAD_ALIGNMENT_XML = "Load Alignment XML";
	/** provide a JMenuBar that can be added to a JFrame
	 * 
	 * @return a JMenuBar
	 */
	public static JMenuBar initMenu(JFrame frame, StructureAlignmentJmol parent, AFPChain afpChain){

		// show a menu

		JMenuBar menu = new JMenuBar();

		JMenu file= new JMenu("File");
		file.getAccessibleContext().setAccessibleDescription("File Menu");

		if ( parent != null){
			JMenuItem loadF = getLoadMenuItem();
			loadF.addActionListener(new MyAlignmentLoadListener(parent));
			file.add(loadF);
		}

		JMenuItem saveF = getSaveAlignmentMenuItem(afpChain);
		file.add(saveF);


		JMenuItem openI = getOpenPDBMenuItem();

		file.add(openI);

		if ( parent != null){
			JMenuItem exportI =  getExportPDBMenuItem(parent);

			file.add(exportI);
		}

		JMenuItem openDBI = getDBResultMenuItem();
		file.add(openDBI);
		file.addSeparator();

		if ( parent != null){
			JMenuItem print = getPrintMenuItem();
			print.addActionListener(parent.getJmolPanel());

			file.add(print);
		}
		file.addSeparator();

		JMenuItem closeI = getCloseMenuItem(frame);

		file.add(closeI);

		JMenuItem exitI = getExitMenuItem();		
		file.add(exitI);
		menu.add(file);


		JMenu align = new JMenu("Align");
		JMenuItem pairI = getPairwiseAlignmentMenuItem();
		align.add(pairI);
		JMenuItem dbI = getDBSearchMenuItem();
		align.add(dbI);
		menu.add(align);

		JMenu view = new JMenu("View");
		view.getAccessibleContext().setAccessibleDescription("View Menu");
		view.setMnemonic(KeyEvent.VK_V);

		if ( parent != null){
			JMenuItem aligpI = MenuCreator.getIcon(parent,ALIGNMENT_PANEL);
			view.add(aligpI);

			JMenuItem textI = MenuCreator.getIcon(parent,TEXT_ONLY);
			view.add(textI);
		}
		
		if ( afpChain != null){
			JMenuItem distMax = new  JMenuItem("Show Distance Matrices");
			distMax.addActionListener(new MyDistMaxListener(afpChain));
			view.add(distMax);
			
			JMenuItem dotplot = new JMenuItem("Show Dot Plot");

			dotplot.addActionListener(new DotPlotListener(afpChain,null));
			view.add(dotplot);
		}

		menu.add(view);


		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_A);

		JMenuItem helpM = getHelpMenuItem();
		about.add(helpM);

		JMenuItem aboutM = getAboutMenuItem();
		about.add(aboutM);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;

	}


	public static JMenuItem getDBResultMenuItem() {
		JMenuItem item = getIcon(new ActionListener() {

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
		}, LOAD_DB_RESULTS );

		return item;
	}


	public static JMenuItem getOpenPDBMenuItem() {
		ImageIcon loadI = createImageIcon("/icons/background.png");
		JMenuItem openI = null;

		if ( loadI == null)
			openI =new JMenuItem("Open PDB file");
		else 
			openI =new JMenuItem("Open PDB file", loadI);
		openI.setMnemonic(KeyEvent.VK_O);
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

		return loadF;
	}


	public static JMenuBar getAlignmentTextMenu(JFrame frame, ActionListener actionListener,AFPChain afpChain){


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
		MySaveFileListener listener = new MySaveFileListener(afpChain);
		listener.setFatCatOutput(true);
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




	private static JMenuItem getIcon(ActionListener actionListener, String text) {
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

		return print;
	}

	public static JMenuItem getExportPDBMenuItem(StructureAlignmentJmol parent) {
		ImageIcon saveicon = createImageIcon("/icons/compfile.png");
		JMenuItem exportI = null;

		if ( saveicon == null)
			exportI = new JMenuItem("Export PDB file");
		else 
			exportI = new JMenuItem("Export PDB file", saveicon);

		exportI.setMnemonic(KeyEvent.VK_E);

		exportI.addActionListener(new MyExportListener(parent));
		return exportI;
	}

	public static JMenuItem getSaveAlignmentMenuItem(AFPChain afpChain) {
		ImageIcon saveicon = createImageIcon("/icons/filesave.png");
		JMenuItem saveF = null;

		if (saveicon == null) 
			saveF = new JMenuItem(SAVE_ALIGNMENT_XML);
		else 
			saveF = new JMenuItem(SAVE_ALIGNMENT_XML, saveicon);

		saveF.setMnemonic(KeyEvent.VK_S);
		saveF.addActionListener(new MySaveFileListener(afpChain));

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

			public void actionPerformed(ActionEvent e) {
				MenuCreator.showAboutDialog();

			}			
		});
		return aboutM;
	}

	public static JMenuItem getExitMenuItem(){

		ImageIcon exitIcon = createImageIcon("/icons/exit.png");


		JMenuItem exitI;

		if ( exitIcon == null)
			exitI = new JMenuItem("Exit");
		else 
			exitI = new JMenuItem("Exit",exitIcon);
		exitI.setMnemonic(KeyEvent.VK_X);
		exitI.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("Exit")){
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

		class MyCloseListener implements ActionListener{
			JFrame f;
			MyCloseListener(JFrame frame){
				f=frame;
			}
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


	/** provide a display for the pairwise structure alignment
	 * 
	 */
	private static void showPairDialog(){
		AlignmentGui gui =  AlignmentGui.getInstance();
		gui.setVisible(true);
	}

	/** show some info about this gui
	 * 
	 */
	public static void showAboutDialog(){

		AboutDialog dialog = new AboutDialog();
		dialog.showDialog();


	}

	/** Returns an ImageIcon, or null if the path was invalid. 
	 * @param path the path to the icon
	 * @return ImageIcon object*/
	public static ImageIcon createImageIcon(String path) {


		java.net.URL imgURL = MenuCreator.class.getResource(path);

		if (imgURL != null) {
			return new ImageIcon(imgURL);
		} else {
			System.err.println("Couldn't find file: " + path);
			return null;
		}

	}


	private static JMenuItem getPairwiseAlignmentMenuItem() {
		ImageIcon alignIcon = createImageIcon("/icons/window_new.png");

		JMenuItem pairI ;
		if ( alignIcon == null) 
			pairI = new JMenuItem("Pairwise Alignment");
		else 
			pairI = new JMenuItem("Pairwise Alignment", alignIcon);
		pairI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("Pairwise Alignment")){
					MenuCreator.showPairDialog();
				}				
			}
		});
		return pairI;
	}
	
	/**
	 * Creates a frame to display a DotPlotPanel.
	 * 
	 * Used by the 'View>Show Dot Plot' menu item
	 * @author Spencer Bliven
	 *
	 */
	private static class DotPlotListener implements ActionListener {
		private final AFPChain afpChain;
		private final Matrix background;
		public DotPlotListener(AFPChain afpChain, Matrix background) {
			this.afpChain = afpChain;
			this.background = background;
		}
		public void actionPerformed(ActionEvent e) {
			String title = String.format("%s vs. %s", afpChain.getName1(),afpChain.getName2());

			// Create window
			JFrame frame = new JFrame(title);
			frame.addWindowListener(new WindowAdapter(){
				public void windowClosing(WindowEvent e){
					JFrame f = (JFrame) e.getSource();
					f.setVisible(false);
					f.dispose();
				}
			});

			DotPlotPanel dotplot = new DotPlotPanel(afpChain, background);			

			frame.getContentPane().add(dotplot);

			frame.pack();
			frame.setVisible(true);
		}
	}


	public static JMenuItem getDBSearchMenuItem() {
		JMenuItem dbI = null;

		ImageIcon dbSearchIcon = createImageIcon("/icons/kpdf.png");

		if ( dbSearchIcon == null )
			dbI = new JMenuItem("DB search");
		else {
			dbI = new JMenuItem("DB search", dbSearchIcon);
		}

		dbI.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				DBSearchGUI.getInstance();

			}			
		});


		return dbI;

	}


	public static JMenuBar initAlignmentGUIMenu(JFrame frame) {


		JMenu file= new JMenu("File");
		file.getAccessibleContext().setAccessibleDescription("File Menu");

		JMenuItem loadF = MenuCreator.getLoadMenuItem();
		loadF.addActionListener(new MyAlignmentLoadListener(null));
		file.add(loadF);

		JMenuItem openI = MenuCreator.getOpenPDBMenuItem();		
		file.add(openI);

		JMenuItem dbI = MenuCreator.getDBResultMenuItem();
		file.add(dbI);
		file.addSeparator();

		JMenuItem closeI = MenuCreator.getCloseMenuItem(frame);
		file.add(closeI);
		JMenuItem exitI = MenuCreator.getExitMenuItem();		
		file.add(exitI);

		JMenuBar menu = new JMenuBar();
		menu.add(file);

		JMenu alig = new JMenu("Align");
		menu.add(alig);

		JMenuItem pw = MenuCreator.getPairwiseAlignmentMenuItem();
		alig.add(pw);

		JMenuItem dbF = MenuCreator.getDBSearchMenuItem();
		alig.add(dbF);

		JMenu about = new JMenu("Help");
		about.setMnemonic(KeyEvent.VK_A);

		JMenuItem aboutM = MenuCreator.getAboutMenuItem();
		about.add(aboutM);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;

	}








}
