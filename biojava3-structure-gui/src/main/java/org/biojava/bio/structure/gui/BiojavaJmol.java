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
 * Created on 24.05.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.bio.structure.gui;


import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JTextField;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.align.gui.jmol.JmolPanel;
import org.biojava.bio.structure.gui.BiojavaJmol;

import org.biojava.bio.structure.gui.util.MenuCreator;


/** A class that provides a simple GUI for Jmol
 * 
 * @author Andreas Prlic
 * @since 1.6
 *
 *
 *
 */
public class BiojavaJmol  {

	public static final String viewer       = "org.jmol.api.JmolSimpleViewer";
	public static final String adapter      = "org.jmol.api.JmolAdapter";
	public static final String smartAdapter = "org.jmol.adapter.smarter.SmarterJmolAdapter";

	Structure structure; 

	JmolPanel jmolPanel;
	JFrame frame ;


	public static void main(String[] args){
		try {

			PDBFileReader pdbr = new PDBFileReader();

			pdbr.setAutoFetch(true);
			pdbr.setPath("/tmp/");

			String pdbCode = "5pti";

			Structure struc = pdbr.getStructureById(pdbCode);

			BiojavaJmol jmolPanel = new BiojavaJmol();

			jmolPanel.setStructure(struc);

			// send some RASMOL style commands to Jmol
			jmolPanel.evalString("select * ; color chain;");
			jmolPanel.evalString("select *; spacefill off; wireframe off; backbone 0.4;  ");
			jmolPanel.evalString("save STATE state_1");
		} catch (Exception e){
			e.printStackTrace();
		}
	}




	public BiojavaJmol() {		

		frame = new JFrame();

		JMenuBar menu = MenuCreator.initMenu();

		frame.setJMenuBar(menu);

		frame.addWindowListener( new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				frame.dispose();
				//System.exit(0);
			}
		});

		Container contentPane = frame.getContentPane();

		Box vBox = Box.createVerticalBox();

		jmolPanel = new JmolPanel();
	
		jmolPanel.setPreferredSize(new Dimension(500,500));
		vBox.add(jmolPanel);


		JTextField field = new JTextField();

		field.setMaximumSize(new Dimension(Short.MAX_VALUE,30));   
		field.setText("enter RASMOL like command...");
		org.biojava.bio.structure.align.gui.jmol.RasmolCommandListener listener = new org.biojava.bio.structure.align.gui.jmol.RasmolCommandListener(jmolPanel,field) ;

		field.addActionListener(listener);
		field.addMouseListener(listener);
		field.addKeyListener(listener);
		vBox.add(field);


		/// COMBO BOXES 
		Box hBox1 = Box.createHorizontalBox();


		String[] styles = new String[] { "Cartoon", "Backbone", "CPK", "Ball and Stick", "Ligands","Ligands and Pocket"};
		JComboBox style = new JComboBox(styles);

		hBox1.add(new JLabel("Style"));
		hBox1.add(style);
		vBox.add(hBox1);
		

		style.addActionListener(jmolPanel);

		String[] colorModes = new String[] { "Secondary Structure", "By Chain", "Rainbow", "By Element", "By Amino Acid", "Hydrophobicity" };
		JComboBox colors = new JComboBox(colorModes);
		colors.addActionListener(jmolPanel);
		hBox1.add(Box.createGlue());
		hBox1.add(new JLabel("Color"));
		hBox1.add(colors);
		
		// Check boxes
		Box hBox2 = Box.createHorizontalBox();
		
		
		JButton resetDisplay = new JButton("Reset Display");
		
		resetDisplay.addActionListener(new ActionListener() {
			
			public void actionPerformed(ActionEvent e) {
				System.out.println("reset!!");
				jmolPanel.executeCmd("restore STATE state_1");
				
			}
		});
		
		hBox2.add(resetDisplay); hBox2.add(Box.createGlue());
		
		JCheckBox toggleSelection = new JCheckBox("Show Selection");
		toggleSelection.addItemListener(
			    new ItemListener() {
					
					public void itemStateChanged(ItemEvent e) {
						  boolean showSelection = (e.getStateChange() == ItemEvent.SELECTED);
						  
						  if (showSelection){
							  jmolPanel.executeCmd("set display selected");
						  } else {
							  jmolPanel.executeCmd("set display off");
						  }
						
					}
				}
			);
		
		
		
		hBox2.add(toggleSelection);
		
		hBox2.add(Box.createGlue());
		vBox.add(hBox2);
		
	
		// finish up
		contentPane.add(vBox);
		frame.pack();
		frame.setVisible(true); 
	
	
	}



	/** returns true if Jmol can be found in the classpath, otherwise false.
	 * 
	 * @return true/false depending if Jmol can be found
	 */
	public static boolean jmolInClassPath(){
		try {
			Class.forName(viewer);		
		} catch (ClassNotFoundException e){
			e.printStackTrace();			
			return false;
		}
		return true;
	}

	public void evalString(String rasmolScript){
		if ( jmolPanel == null ){
			System.err.println("please install Jmol first");
			return;
		}
		jmolPanel.evalString(rasmolScript);
	}

	public void setStructure(Structure s) {

		if ( jmolPanel == null ){
			System.err.println("please install Jmol first");
			return;
		}

		setTitle(s.getPDBCode());

		// actually this is very simple
		// just convert the structure to a PDB file

		String pdb = s.toPDB();	
		//System.out.println(s.isNmr());

		//System.out.println(pdb);
		// Jmol could also read the file directly from your file system
		//viewer.openFile("/Path/To/PDB/1tim.pdb");

		//System.out.println(pdb);
		jmolPanel.openStringInline(pdb);

		// send the PDB file to Jmol.
		// there are also other ways to interact with Jmol, e.g make it directly
		// access the biojava structure object, but they require more
		// code. See the SPICE code repository for how to do this.

		


	}

	public void setTitle(String label){
		frame.setTitle(label);
		frame.repaint();
	}

	public JFrame getFrame(){
		return frame;
	}


}