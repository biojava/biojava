package org.biojava.nbio.structure.symmetry.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;

/**
 * Action Listener for the symmetry menu.
 * Trigger various symmetry analysis.
 * 
 * @author Aleix Lafita
 *
 */
public class SymmetryListener implements ActionListener{

	private MultipleAlignmentJmol jmol;
	private MultipleAlignment msa;
	private SymmetryAxes axes;

	public SymmetryListener(MultipleAlignmentJmol jmol, SymmetryAxes axes) {
		this.jmol = jmol;
		this.msa = jmol.getMultipleAlignment();
		this.axes = axes;
	}

	@Override
	public void actionPerformed(ActionEvent ae) {

		String cmd = ae.getActionCommand();
		if (cmd.equals("Subunit Superposition")){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				MultipleAlignmentJmol j = SymmetryDisplay.displaySubunits(msa);
				String s = SymmetryDisplay.printSymmetryAxes(msa, axes, true);
				j.evalString(s);
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Multiple Structure Alignment")){
			if (msa == null) {
				System.err.println("Currently not displaying a symmetry!");
				return;
			}
			try {
				MultipleAlignmentJmol j = SymmetryDisplay.displayFull(msa);
				String s = SymmetryDisplay.printSymmetryAxes(msa, axes, false);
				j.evalString(s);
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Point Group Symmetry")){
			if (msa != null) {
				String script = SymmetryDisplay.printPointGroupAxes(msa);
				jmol.evalString(script);
				return;
			}
			
		} else if (cmd.equals("Show Symmetry Axes")){
			if (axes != null) {
				String s = SymmetryDisplay.printSymmetryAxes(msa, axes, false);
				jmol.evalString(s);
				return;
			} else System.err.println("No axes for this symmetry");
			
		} else if (cmd.equals("New Symmetry Analysis")){
			SymmetryGui.getInstance();
		}
	}
	
}