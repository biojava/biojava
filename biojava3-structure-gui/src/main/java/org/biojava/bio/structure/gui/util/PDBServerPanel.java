
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
 * Created on Feb 9, 2007
 *
 */
package org.biojava.bio.structure.gui.util;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.io.IOException;
import java.util.logging.Logger;
import javax.swing.BorderFactory;
import javax.swing.Box;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;


import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.io.PDBFileReader;

/** A class to define where a structure for the alignment is coming from
 *
 * @author Andreas Prlic
 * @since 1.7
 * @version %I% %G%
 */
public class PDBServerPanel
extends JPanel
implements StructurePairSelector{

	/**
	 *
	 */
	private static final long serialVersionUID = -5682120627824627408L;

	boolean debug = true;
	JTextField pdbDir;
	JTextField f1;
	JTextField f2;
	JTextField c1;
	JTextField c2;



	public static Logger logger =  Logger.getLogger("org.biojava");

	/** load the PDB files from a local directory
	 *
	 */
	public PDBServerPanel() {

		Box vBox = Box.createVerticalBox();

		int pdbfSize = 4;

		f1 = new JTextField(pdbfSize);
		c1 = new JTextField(1);
		JPanel p1 = getPDBFilePanel(1,f1,c1);
		vBox.add(p1);

		f2 = new JTextField(pdbfSize);
		c2 = new JTextField(1);
		JPanel p2 = getPDBFilePanel(2, f2,c2);
		vBox.add(p2);


		this.add(vBox);

	}


	private Structure fromPDB(JTextField f, JTextField c) throws StructureException{
		String pdb = f.getText();


		if ( pdb.length() < 4) {
			f.setText("!!!");
			return null;
		}

		String chain = c.getText();
		if ( debug )
			System.out.println("file :" + pdb + " " +  chain);
		/// prepare structures

		// load them from the file system

		PDBFileReader reader = new PDBFileReader();

		reader.setPath(".");
		reader.setAutoFetch(true);

		Structure tmp1 = new StructureImpl();

		try {
			Structure structure1 = reader.getStructureById(pdb);

			// no chain has been specified
			// return whole structure
			if (( chain == null) || (chain.length()==0)){
				return structure1;
			}
			if ( debug)
				System.out.println("using chain " + chain +  " for structure " + structure1.getPDBCode());
			Chain c1 = structure1.findChain(chain);
			tmp1.setPDBCode(structure1.getPDBCode());
			tmp1.setPDBHeader(structure1.getPDBHeader());
			tmp1.setPDBCode(structure1.getPDBCode());
			tmp1.addChain(c1);
			System.out.println("ok");

		} catch (IOException e){
			logger.warning(e.getMessage());
			throw new StructureException(e);
		}
		return tmp1;
	}



	public Structure getStructure1() throws StructureException{
		return fromPDB(f1,c1);
	}

	public Structure getStructure2() throws StructureException{
		return fromPDB(f2,c2);
	}


//	private JPanel getPDBDirPanel(JTextField f){
//		JPanel panel = new JPanel();
//		panel.setBorder(BorderFactory.createLineBorder(Color.black));
//
//
//		JLabel l01 = new JLabel("Select PDB directory");
//		panel.add(l01);
//
//		panel.add(f);
//
//		Action action = new ChooseDirAction(pdbDir);
//
//		JButton chooser = new JButton(action);
//		panel.add(chooser);
//		return panel;
//	}

	private JPanel getPDBFilePanel(int pos ,JTextField f, JTextField c){

		JPanel panel = new JPanel();
		panel.setBorder(BorderFactory.createLineBorder(Color.black));



		JLabel l01 = new JLabel("PDB code ");

		panel.add(l01);
		Box hBox11 = Box.createHorizontalBox();

		JLabel l11 = new JLabel(pos + ":");



		f.setMaximumSize(new Dimension(Short.MAX_VALUE,30));


		hBox11.add(l11);
		hBox11.add(Box.createVerticalGlue());
		hBox11.add(f, BorderLayout.CENTER);
		hBox11.add(Box.createVerticalGlue());

		panel.add(hBox11);

		Box hBox21 = Box.createHorizontalBox();
		JLabel l21 = new JLabel("Chain" + pos + ":");

		c.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		hBox21.add(l21);
		hBox21.add(Box.createGlue());
		hBox21.add(c, BorderLayout.CENTER);
		hBox21.add(Box.createGlue());

		panel.add(hBox21);

		return panel;
	}

}


