package org.biojava.nbio.structure.align.gui;

import java.awt.Color;
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
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.StructureAlignmentException;

/**
 * Utility Class that provides methods for the visualization of {@link MultipleAlignment}s.
 * <p>
 * Currently supported: Alignment Panel Display, select aligned residues by PDB code
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentDisplay {

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
	 * @throws StructureAlignmentException
	 * @throws StructureException
	 */
	public static void showMultipleAligmentPanel(MultipleAlignment multAln, AbstractAlignmentJmol jmol, Color[] colors) throws StructureAlignmentException, StructureException {
		
		MultipleAligPanel me = new MultipleAligPanel(multAln, colors, jmol);
		JFrame frame = new JFrame();

		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);		
		frame.setTitle("Alignment Panel for Multiple Structure Alignments - alpha Version");
		me.setPreferredSize(new Dimension(me.getCoordManager().getPreferredWidth() , me.getCoordManager().getPreferredHeight()));

		JMenuBar menu = MenuCreator.getAlignmentTextMenu(frame,me,null);
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
}
