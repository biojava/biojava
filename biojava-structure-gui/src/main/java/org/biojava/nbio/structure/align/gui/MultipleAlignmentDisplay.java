package org.biojava.nbio.structure.align.gui;

import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JScrollPane;
import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.aligpanel.MultipleAligPanel;
import org.biojava.nbio.structure.align.gui.aligpanel.MultipleStatusDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.ReferenceSuperimposer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
public class MultipleAlignmentDisplay {

	private static final Logger logger = 
			LoggerFactory.getLogger(MultipleAlignmentDisplay.class);

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
	public static final List<String> getPDBresnum(int structNum, 
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
	 * @param colors
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

		int size = multAln.size();

		List<Atom[]> atomArrays = multAln.getAtomArrays();
		for (int i=0; i<size; i++){
			if (atomArrays.get(i).length < 1) 
				throw new StructureException(
						"Length of atoms arrays is too short! Size: " 
								+ atomArrays.get(i).length);
		}

		List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();

		//TODO implement independent BlockSet superposition of the structure
		List<Matrix4d> transf = multAln.getBlockSet(0).getTransformations();

		if(transf == null) {

			logger.error("Alignment Transformations are not calculated. "
					+ "Superimposing to first structure as reference.");

			multAln = multAln.clone();
			MultipleSuperimposer imposer = new ReferenceSuperimposer();
			imposer.superimpose(multAln);
			transf = multAln.getBlockSet(0).getTransformations();
			assert(transf != null);
		}

		//Rotate the atom coordinates of all the structures
		for (int i=0; i<size; i++){
			//TODO handle BlockSet-level transformations
			//make sure this method has the same behavior as the other display.
			//-SB 2015-06

			//Assume all atoms are from the same structure
			Structure displayS = atomArrays.get(i)[0].getGroup().
					getChain().getParent().clone();
			//Get all the atoms and include ligands and hetatoms
			Atom[] rotCA = StructureTools.getRepresentativeAtomArray(displayS);
			List<Group> hetatms = StructureTools.getUnalignedGroups(rotCA);
			for (Group g:hetatms){
				rotCA = Arrays.copyOf(rotCA, rotCA.length + 1);
				rotCA[rotCA.length - 1] = g.getAtom(0);
			}

			//Transform the structure to ensure a full rotation in the display
			Calc.transform(displayS, transf.get(i));
			rotatedAtoms.add(rotCA);
		}

		MultipleAlignmentJmol jmol = 
				new MultipleAlignmentJmol(multAln, rotatedAtoms);
		jmol.setTitle(jmol.getStructure().getPDBHeader().getTitle());
		return jmol;
	}

	/** 
	 * Get an artifical Structure containing a different model for every
	 * input structure, so that the alignment result can be viewed in Jmol.
	 * The Atoms have to be rotated beforehand.
	 * 
	 * @param atomArrays an array of Atoms for every aligned structure
	 * @return a structure object containing a set of models, 
	 * 			one for each input array of Atoms.
	 * @throws StructureException
	 */
	public static final Structure getAlignedStructure(List<Atom[]> atomArrays) 
			throws StructureException {

		Structure s = new StructureImpl();
		for (int i=0; i<atomArrays.size(); i++){
			List<Chain> model = DisplayAFP.getAlignedModel(atomArrays.get(i));
			s.addModel(model);
		}
		return s;
	}

}
