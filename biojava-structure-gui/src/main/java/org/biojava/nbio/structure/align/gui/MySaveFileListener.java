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

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;

import javax.swing.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

/**
 * Save an alignment to a specified File by the user. The alignment to be
 * saved depends on the constructor used to instantiate this class, so that
 * AFPChains and MultipleAlignments can be saved.
 * <p>
 * The format to save the alignment depends on the Frame that generated the
 * Action: from a sequence alignment a FatCat or FASTA formats are saved,
 * and from a Jmol view an XML format is saved.
 *
 * @author Aleix Lafita
 * @version 2.0 - adapted for MultipleAligments
 *
 */
public class MySaveFileListener implements ActionListener {

	private AFPChain afpChain;
	private MultipleAlignment msa;
	private boolean printText;

	public MySaveFileListener (AFPChain afpChain){
		this.afpChain = afpChain;
		this.msa = null;
		printText = false;
	}

	public MySaveFileListener (MultipleAlignment msa){
		this.afpChain = null;
		this.msa = msa;
		printText = false;
	}

	/**
	 * Constructor to avoid checking which of the two is null before
	 * instantiating this class. One of the two, or both, have to be null.
	 * If both are different than null the MultipleAlignment will be saved
	 * only.
	 *
	 * @param afpChain
	 * @param msa
	 */
	public MySaveFileListener (AFPChain afpChain, MultipleAlignment msa){
		this.afpChain = afpChain;
		this.msa = msa;
		printText = false;
	}

	@Override
	public void actionPerformed(ActionEvent evt) {

		//Return if nothing to save
		if (afpChain == null && msa == null) {
			JOptionPane.showMessageDialog(null,
					"Could not save alignment, no alignment being displayed.");
			return;
		}

		//Choose the file path from the user input
		JFileChooser fc = new JFileChooser();
		int returnVal = fc.showSaveDialog(null);

		if (returnVal != JFileChooser.APPROVE_OPTION) {
			System.err.println("User canceled file save.");
			return;
		}
		File selFile = fc.getSelectedFile();
		if (selFile == null) return;

		System.out.println("Saving alignment to file: " + selFile.getName());

		//XML serialization of the alignment
		try {

			String output = "";

			if (msa!=null){
				if (!printText) {
					output = MultipleAlignmentWriter.toXML(msa.getEnsemble());
				} else {
					output = MultipleAlignmentWriter.toFASTA(msa);
				}
			}
			else if (afpChain!=null){
				UserConfiguration config = WebStartMain.getWebStartConfig();
				AtomCache cache = new AtomCache(config);

				//TODO use the right ca atoms, fails for a custom file!
				//This is a bad solution, solved for MultipleAlignments
				Atom[] ca1 =cache.getAtoms(afpChain.getName1());
				Atom[] ca2 =cache.getAtoms(afpChain.getName2());

				if (!printText) {
					output = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
				} else {
					output = afpChain.toFatcat(ca1,ca2);
				}
			}

			//Write to the file
			BufferedWriter out = new BufferedWriter(new FileWriter(selFile));
			out.write(output);
			out.close();

		} catch (Exception e){
			e.printStackTrace();
			JOptionPane.showMessageDialog(null,
					"Could not save file. Exception: " + e.getMessage());
		}
	}

	/**
	 * If true, the alignment format saved will be a text output (FASTA for
	 * MultipleAlignments and FatCat for pairwise alignments)
	 *
	 * @param text if true the output will be text format
	 */
	public void setTextOutput(boolean text) {
		printText = text;
	}

}
