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
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.biojava.nbio.structure.align.xml.MultipleAlignmentXMLParser;
import org.biojava.nbio.core.util.InputStreamProvider;

import javax.swing.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

/**
 * Loads an alignment in an XML format and displays its content in a
 * new Jmol panel. Can handle both alignment formats: AFPChain and
 * MultipleAlignment.
 * <p>
 * All the alignments stored in the File are displayed, not only the first
 * one. However, all the alignments in the same file have to be in the same
 * format (either AFPChain or MultipleAlignment).
 * Multiple Jmol panels can be created for that purpose.
 *
 * @author Aleix Lafita
 * @version 2.0 - adapted for MultipleAlignments
 *
 */
public class MyAlignmentLoadListener implements ActionListener {

	@Override
	public void actionPerformed(ActionEvent evt) {

		final JFileChooser fc = new JFileChooser();

		//in response to a button click
		int returnVal = fc.showOpenDialog(null);

		if (returnVal == JFileChooser.APPROVE_OPTION) {

			File file = fc.getSelectedFile();
			try {

				InputStreamProvider ip = new InputStreamProvider();
				InputStream stream = ip.getInputStream(file);
				BufferedReader in = new BufferedReader(
						new InputStreamReader(stream));

				StringBuffer input = new StringBuffer();
				String str;
				while ((str = in.readLine()) != null) {
					input.append(str);
				}
				in.close();

				String xml = input.toString();

				//Determine the format of the file
				if (xml.contains("MultipleAlignmentEnsemble")){

					List<MultipleAlignmentEnsemble> ensembles =
							MultipleAlignmentXMLParser.parseXMLfile(xml);

					//Display all ensembles, and all its alignments
					for (MultipleAlignmentEnsemble e:ensembles){
						for (MultipleAlignment msa:e.getMultipleAlignments()){
							MultipleAlignmentJmolDisplay.display(msa);
						}
					}

				}
				else {

					AFPChain[] afps = AFPChainXMLParser.parseMultiXML(xml);

					UserConfiguration conf = WebStartMain.getWebStartConfig();
					AtomCache cache = new AtomCache(
							conf.getPdbFilePath(),conf.getCacheFilePath());

					for (AFPChain afpChain:afps){
						Atom[] ca1 = cache.getAtoms(afpChain.getName1());
						Atom[] ca2 = cache.getAtoms(afpChain.getName2());

						AFPChainXMLParser.rebuildAFPChain(afpChain, ca1, ca2);
						StructureAlignmentJmol jmol =
								StructureAlignmentDisplay.display(
										afpChain, ca1, ca2);

						DisplayAFP.showAlignmentPanel(afpChain, ca1,ca2,jmol);
					}
				}
			} catch (Exception e){
				e.printStackTrace();
				JOptionPane.showMessageDialog(null,"Could not load alignment "
						+ "file. Exception: " + e.getMessage());
			}
		}
	}

}
