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
package org.biojava.nbio.structure.align.gui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.IOException;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;
import org.biojava.nbio.structure.gui.util.StructurePairSelector;


/** A Panel that allows user to specify PDB & chain ID, as well as sub-ranges
 *
 * @author Andreas
 *
 */
public class SelectPDBPanel
extends JPanel
implements StructurePairSelector{

	boolean debug = true;

	JTextField f1;
	JTextField f2;
	JTextField c1;
	JTextField c2;
	JTextField r1;
	JTextField r2;

	UserConfiguration config;
	JTabbedPane configPane;

	/**
	 *
	 */
	private static final long serialVersionUID = 4002475313717172193L;



	public SelectPDBPanel(){
		this(true);
	}
	public SelectPDBPanel(boolean show2PDBs) {

		Box vBox = Box.createVerticalBox();

		JLabel help = new JLabel("Optional: specify chain ID or range.");
		Box hBox1 = Box.createHorizontalBox();
		hBox1.add(Box.createGlue());
		hBox1.add(help);
		vBox.add(hBox1);


		//pdbDir = new JTextField(20);

		int pdbfSize = 4;

		f1 = new JTextField(pdbfSize);
		c1 = new JTextField(1);
		r1 = new JTextField(5);
		Box p1 = getPDBFilePanel(1,f1,c1,r1);
		vBox.add(p1);

		f2 = new JTextField(pdbfSize);
		c2 = new JTextField(1);
		r2 = new JTextField(5);
		Box p2 = getPDBFilePanel(2, f2,c2,r2);

		if ( show2PDBs)
			vBox.add(p2);

		//vBox.setBorder(BorderFactory.createLineBorder(Color.black));
		this.add(vBox);
	}

	public StructureIdentifier getName1() {
		String pdbId = f1.getText().trim();
		String chainId = c1.getText().trim();
		String range = r1.getText().trim();

		// Prefer range over chain
		if( ! range.isEmpty() ) {
			return new SubstructureIdentifier(pdbId, ResidueRange.parseMultiple(range));
		} else if ( ! chainId.isEmpty() ){
			return new SubstructureIdentifier(pdbId, ResidueRange.parseMultiple(chainId));
		}
		return new SubstructureIdentifier(pdbId);
	}
	public StructureIdentifier getName2() {
		String pdbId = f2.getText().trim();
		String chainId = c2.getText().trim();
		String range = r2.getText().trim();

		// Prefer range over chain
		if( ! range.isEmpty() ) {
			return new SubstructureIdentifier(pdbId, ResidueRange.parseMultiple(range));
		} else if ( ! chainId.isEmpty() ){
			return new SubstructureIdentifier(pdbId, ResidueRange.parseMultiple(chainId));
		}
		return new SubstructureIdentifier(pdbId);
	}
	@Override
	public Structure getStructure1() throws StructureException, IOException{
		return getStructure(getName1());
	}

	@Override
	public Structure getStructure2() throws StructureException, IOException{
		return getStructure(getName2());
	}

	private Structure getStructure(StructureIdentifier name) throws IOException, StructureException {
		UserConfiguration config = WebStartMain.getWebStartConfig();
		AtomCache cache = new AtomCache(config);
		return cache.getStructure(name);
	}

	private Box getPDBFilePanel(int pos ,JTextField f, JTextField c, JTextField r){

		//JPanel panel = new JPanel();
		//panel.setBorder(BorderFactory.createLineBorder(Color.black));

		JLabel l01 = new JLabel("PDB code ");

		//panel.add(l01);
		Box hBox = Box.createHorizontalBox();
		hBox.add(Box.createGlue());
		hBox.add(l01);

		JLabel l11 = new JLabel(pos + ":");
		f.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		f.setToolTipText("Provide 4-character PDB code here. Example: 4hhb");
		hBox.add(l11);
		hBox.add(Box.createVerticalGlue());
		hBox.add(f, BorderLayout.CENTER);
		hBox.add(Box.createGlue());

		//panel.add(hBox11);

		//Box hBox21 = Box.createHorizontalBox();
		JLabel l21 = new JLabel("Chain" + pos + ":");
		hBox.add(l21);

		c.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		//hBox.add(Box.createGlue());
		hBox.add(c, BorderLayout.CENTER);

		String msg1 = "Both chainID and range specification are optional. If both are provided, range has preference.";
		l21.setToolTipText(msg1);
		c.setToolTipText(msg1);

		JLabel rangeL = new JLabel(" Range " + pos + ":");
		hBox.add(Box.createGlue());
		hBox.add(rangeL);
		r.setMaximumSize(new Dimension(Short.MAX_VALUE,30));

		// set help text:
		String msg ="Syntax example: A:407-495,A:582-686";
		rangeL.setToolTipText(msg);
		r.setToolTipText(msg);

		//hBox.add(Box.createGlue());
		hBox.add(r,BorderLayout.CENTER);

		//hBox21.add(Box.createGlue());

		//panel.add(hBox21);



		return hBox;
	}





}
