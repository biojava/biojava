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
package org.biojava.bio.structure.gui.util;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;

import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.gui.BiojavaJmol;


/**
 *  Create the menu for BiojavaJmol
 * @author Andreas Prlic
 * @since 1.7
 */
public class MenuCreator {

	/** provide a JMenuBar that can be added to a JFrame
	 * 
	 * @return a JMenuBar
	 */
	public static JMenuBar initMenu(){

		// show a menu

		JMenuBar menu = new JMenuBar();

		JMenu file= new JMenu("File");
		file.getAccessibleContext().setAccessibleDescription("File Menu");

		JMenuItem openI = new JMenuItem("Open");
		openI.setMnemonic(KeyEvent.VK_O);
		openI.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();			        
				if ( cmd.equals("Open")){
					final JFileChooser fc = new JFileChooser();

//					In response to a button click:
					int returnVal = fc.showOpenDialog(null);
					if ( returnVal == JFileChooser.APPROVE_OPTION) {
						File file = fc.getSelectedFile();
						
						PDBFileReader reader = new PDBFileReader();
						try {
							Structure s = reader.getStructure(file);
							BiojavaJmol jmol = new BiojavaJmol();
							jmol.setStructure(s);
							jmol.evalString("select * ; color chain;");
							jmol.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");

						} catch (Exception ex){
							ex.printStackTrace();
						}
						

					}
				}				
			}		
		});
		file.add(openI);

		JMenuItem exitI = new JMenuItem("Exit");
		exitI.setMnemonic(KeyEvent.VK_X);
		exitI.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("Exit")){
					System.exit(0);
				}				
			}			
		});

		file.add(exitI);
		menu.add(file);

		JMenu align = new JMenu("Align");
		JMenuItem pairI = new JMenuItem("2 protein structures");
		pairI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("2 protein structures")){
					MenuCreator.showPairDialog();
				}				
			}
		});

		align.add(pairI);

		menu.add(align);

		JMenu about = new JMenu("About");
		JMenuItem aboutI = new JMenuItem("PDBview");
		aboutI.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				String cmd = e.getActionCommand();

				if ( cmd.equals("PDBview")){
					MenuCreator.showAboutDialog();
				}				
			}			
		});

		about.add(aboutI);

		menu.add(Box.createGlue());
		menu.add(about);

		return menu;

	}


	/** provide a display for the pairwise structure alignment
	 * 
	 */
	private static void showPairDialog(){
		AlignmentGui gui = AlignmentGui.getInstance();
		gui.setVisible(true);
	}

	/** show some info about this gui
	 * 
	 */
	private static void showAboutDialog(){

		JDialog dialog = new JDialog();

		dialog.setSize(new Dimension(300,300));

		String msg = "This viewer is based on <b>BioJava</b> and <b>Jmol</>. <br>Author: Andreas Prlic <br> ";
		msg += "Structure Alignment algorithm based on a variation of the PSC++ algorithm by Peter Lackner.";


		JEditorPane txt = new JEditorPane("text/html", msg);
		txt.setEditable(false);


		JScrollPane scroll = new JScrollPane(txt);

		Box vBox = Box.createVerticalBox();
		vBox.add(scroll);

		JButton close = new JButton("Close");

		close.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent event) {
				Object source = event.getSource();

				JButton but = (JButton)source;
				Container parent = but.getParent().getParent().getParent().getParent().getParent().getParent() ;

				JDialog dia = (JDialog) parent;
				dia.dispose();
			}
		});

		Box hBoxb = Box.createHorizontalBox();
		hBoxb.add(Box.createGlue());
		hBoxb.add(close,BorderLayout.EAST);

		vBox.add(hBoxb);

		dialog.getContentPane().add(vBox);
		dialog.setVisible(true);


	}


}
