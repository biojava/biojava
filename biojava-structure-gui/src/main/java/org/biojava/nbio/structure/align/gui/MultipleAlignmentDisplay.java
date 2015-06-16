package org.biojava.nbio.structure.align.gui;

import java.awt.Dimension;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JScrollPane;
import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.aligpanel.MultipleAligPanel;
import org.biojava.nbio.structure.align.gui.aligpanel.MultipleStatusDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleSuperimposer;
import org.biojava.nbio.structure.align.multiple.ReferenceSuperimposer;
import org.jcolorbrewer.ColorBrewer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Utility Class that provides methods for the visualization of {@link MultipleAlignment}s.
 * <p>
 * Currently supported: Alignment Panel Display, select aligned residues by PDB code, 
 * JmolPanel display.
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentDisplay {
	
	private static final Logger logger = LoggerFactory.getLogger(MultipleAlignmentDisplay.class);

	/**
	 * Utility method used in the {@link MultipleAlignmentJmol} panel, when the aligned residues of a 
	 * structure in the alignment have to be selected for formatting (coloring and style).
	 * 
	 * @param structNum the structure index (row) of the alignment
	 * @param multAln the MultipleAlignment that contains the equivalent positions
	 * @param ca the atom array of the structure specified (corresponding to the structure index)
	 * @return List of pdb Strings corresponding to the aligned positions of the molecule.
	 */
	public static final List<String> getPDBresnum(int structNum, MultipleAlignment multAln, Atom[] ca){
		
		List<String> lst = new ArrayList<String>();

		//Loop through all the Blocks in the alignment
		for(Block block : multAln.getBlocks() ) {
			//Loop though all the residues in the Block
			for (int i=0; i<block.length(); i++){
				Integer pos = block.getAlignRes().get(structNum).get(i);
				if (pos==null) continue; //It means a GAP
				else if (pos < ca.length) {
					String pdbInfo = JmolTools.getPdbInfo(ca[pos]);
					lst.add(pdbInfo);
				}
			}
		}
		return lst;
	}
	
	/**
	 * Creates a new Frame with the MultipleAlignment Sequence Panel. The panel can communicate with
	 * the Jmol 3D visualization by selecting the aligned residues of all structures.
	 * 
	 * @param multAln
	 * @param jmol
	 * @param colors
	 * @throws StructureException
	 */
	public static void showMultipleAligmentPanel(MultipleAlignment multAln, AbstractAlignmentJmol jmol, ColorBrewer colorPattelete) throws StructureException {
		
		MultipleAligPanel me = new MultipleAligPanel(multAln, jmol);
		JFrame frame = new JFrame();

		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);		
		frame.setTitle(jmol.getTitle());
		me.setPreferredSize(new Dimension(me.getCoordManager().getPreferredWidth() , me.getCoordManager().getPreferredHeight()));

		JMenuBar menu = MenuCreator.getAlignmentPanelMenu(frame,me,null);
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
		//make sure they get cleaned up correctly:
		frame.addWindowListener(me);
		frame.addWindowListener(status);
	}
	
	/**
	 * Creates a new Frame with the String output representation of the {@link MultipleAlignment}.
	 * 
	 * @param multAln
	 * @param result String output
	 */
	public static void showAlignmentImage(MultipleAlignment multAln, String result) {

		JFrame frame = new JFrame();

		String title = multAln.getEnsemble().getAlgorithmName() + " V."+multAln.getEnsemble().getVersion();
		frame.setTitle(title);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		AlignmentTextPanel txtPanel = new AlignmentTextPanel();
		txtPanel.setText(result);

		JMenuBar menu = MenuCreator.getAlignmentTextMenu(frame,txtPanel,null);

		frame.setJMenuBar(menu);
		JScrollPane js = new JScrollPane();
		js.getViewport().add(txtPanel);
		js.getViewport().setBorder(null);

		frame.getContentPane().add(js);
		frame.pack();      
		frame.setVisible(true);
	}
	
   /** 
    * Display a MultipleAlignment with a JmolPanel. New structures are downloaded if they were 
    * not cached in the alignment and they are entirely transformed here with the cached
    * superposition information.
    * 
    * @param multAln
    * @return MultipleAlignmentJmol instance
    * @throws StructureException
    */
   public static MultipleAlignmentJmol display(MultipleAlignment multAln) throws StructureException {

	   int size = multAln.size();

	   List<Atom[]> atomArrays = multAln.getEnsemble().getAtomArrays();
	   for (int i=0; i<size; i++){
		   if (atomArrays.get(i).length < 1) 
			   throw new StructureException("Length of atoms arrays is too short! " + atomArrays.get(i).length);
	   }

	   List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();

	   List<Matrix4d> transformations = multAln.getTransformations();
	   if( transformations == null ) {
		   //TODO temporary hack for missing transformations
		   logger.error("BlockSet transformations are unimplemented. Superimposing to first structure.");
		   // clone input, since we're about to re-superimpose it
		   multAln = multAln.clone();
		   MultipleSuperimposer imposer = new ReferenceSuperimposer();
		   imposer.superimpose(multAln);
		   transformations = multAln.getTransformations();
		   assert(transformations != null);
	   }

	   //Rotate the atom coordinates of all the structures
	   for (int i=0; i<size; i++){
		   //TODO handle BlockSet-level transformations for flexible alignments.
		   // In general, make sure this method has the same behavior as the other display. -SB 2015-06

		   // Assume all atoms are from the same structure
		   Structure displayS = atomArrays.get(i)[0].getGroup().getChain().getParent().clone();
		   Atom[] rotCA = StructureTools.getRepresentativeAtomArray(displayS);
		   //Rotate the structure to ensure a full rotation in the display
		   Calc.transform(rotCA[0].getGroup().getChain().getParent(), multAln.getTransformations().get(i));
		   rotatedAtoms.add(rotCA);
	   }


	   MultipleAlignmentJmol jmol = new MultipleAlignmentJmol(multAln, rotatedAtoms);
	   jmol.setTitle(jmol.getStructure().getPDBHeader().getTitle());
	   return jmol;
   }
}
