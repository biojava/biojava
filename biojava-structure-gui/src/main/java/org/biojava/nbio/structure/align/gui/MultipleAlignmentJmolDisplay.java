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

import java.awt.Dimension;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JScrollPane;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.aligpanel.MultipleAligPanel;
import org.biojava.nbio.structure.align.gui.aligpanel.MultipleStatusDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentDisplay;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

/**
 * Utility Class that provides helper methods for the visualization of
 * {@link MultipleAlignment}s.
 * <p>
 * Currently supported: Alignment Panel Display, select aligned
 * residues in Jmol by their PDB name, show a text Frame for any sequence
 * alignment format, basic Jmol display from a MultipleAlignment, generate
 * an artificial PDB structure with a new model for every aligned structure.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class MultipleAlignmentJmolDisplay  {

	//private static final Logger logger =
	//		LoggerFactory.getLogger(MultipleAlignmentJmolDisplay.class);

	/**
	 * Utility method used in the {@link MultipleAlignmentJmol} Frame,
	 * when the aligned residues of a structure in the alignment have
	 * to be selected for formatting them (coloring and style).
	 *
	 * @param structNum the structure index (row) of the alignment
	 * @param multAln the MultipleAlignment that contains the equivalent
	 * 			positions
	 * @param ca the atom array of the structure specified
	 * 			(corresponding to the structure index)
	 * @return List of pdb Strings corresponding to the aligned positions
	 * 			of the structure.
	 */
	public static List<String> getPDBresnum(int structNum,
			MultipleAlignment multAln, Atom[] ca){

		List<String> lst = new ArrayList<String>();

		for(Block block : multAln.getBlocks() ) {

			for (int i=0; i<block.length(); i++){
				Integer pos = block.getAlignRes().get(structNum).get(i);
				if (pos==null) continue; //gap
				else if (pos < ca.length) {
					String pdbInfo = JmolTools.getPdbInfo(ca[pos]);
					lst.add(pdbInfo);
				}
			}
		}
		return lst;
	}

	/**
	 * Creates a new Frame with the MultipleAlignment Sequence Panel.
	 * The panel can communicate with the Jmol 3D visualization by
	 * selecting the aligned residues of every structure.
	 *
	 * @param multAln
	 * @param jmol

	 * @throws StructureException
	 */
	public static void showMultipleAligmentPanel(MultipleAlignment multAln,
			AbstractAlignmentJmol jmol) throws StructureException {

		MultipleAligPanel me = new MultipleAligPanel(multAln, jmol);
		JFrame frame = new JFrame();

		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setTitle(jmol.getTitle());
		me.setPreferredSize(new Dimension(
				me.getCoordManager().getPreferredWidth() ,
				me.getCoordManager().getPreferredHeight()));

		JMenuBar menu = MenuCreator.getAlignmentPanelMenu(
				frame,me,null, multAln);
		frame.setJMenuBar(menu);

		JScrollPane scroll = new JScrollPane(me);
		scroll.setAutoscrolls(true);

		MultipleStatusDisplay status = new MultipleStatusDisplay(me);
		me.addAlignmentPositionListener(status);

		Box vBox = Box.createVerticalBox();
		vBox.add(scroll);
		vBox.add(status);
		frame.getContentPane().add(vBox);

		frame.pack();
		frame.setVisible(true);

		frame.addWindowListener(me);
		frame.addWindowListener(status);
	}

	/**
	 * Creates a new Frame with the String output representation of the
	 * {@link MultipleAlignment}.
	 *
	 * @param multAln
	 * @param result String output
	 */
	public static void showAlignmentImage(MultipleAlignment multAln,
			String result) {

		JFrame frame = new JFrame();

		String title = multAln.getEnsemble().getAlgorithmName() +
				" V."+multAln.getEnsemble().getVersion();
		frame.setTitle(title);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		AlignmentTextPanel txtPanel = new AlignmentTextPanel();
		txtPanel.setText(result);

		JMenuBar menu = MenuCreator.getAlignmentTextMenu(
				frame,txtPanel,null,multAln);

		frame.setJMenuBar(menu);
		JScrollPane js = new JScrollPane();
		js.getViewport().add(txtPanel);
		js.getViewport().setBorder(null);

		frame.getContentPane().add(js);
		frame.pack();
		frame.setVisible(true);
	}

	/**
	 * Display a MultipleAlignment with a JmolPanel.
	 * New structures are downloaded if they were
	 * not cached in the alignment and they are entirely
	 * transformed here with the superposition information
	 * in the Multiple Alignment.
	 *
	 * @param multAln
	 * @return MultipleAlignmentJmol instance
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol display(MultipleAlignment multAln)
			throws StructureException {

		List<Atom[]> rotatedAtoms = MultipleAlignmentDisplay.getRotatedAtoms(multAln);

		MultipleAlignmentJmol jmol =
				new MultipleAlignmentJmol(multAln, rotatedAtoms);

		jmol.setTitle(jmol.getStructure().getPDBHeader().getTitle());
		return jmol;
	}

}
