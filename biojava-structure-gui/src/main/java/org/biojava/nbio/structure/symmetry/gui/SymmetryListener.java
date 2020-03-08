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
package org.biojava.nbio.structure.symmetry.gui;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Action Listener for the symmetry menu. Trigger various internal symmetry
 * specific analysis.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryListener implements ActionListener {

	private final MultipleAlignmentJmol jmol;
	private final CeSymmResult symm;

	private static final Logger logger = LoggerFactory
			.getLogger(SymmetryListener.class);

	public SymmetryListener(MultipleAlignmentJmol jmol, CeSymmResult symm) {
		this.jmol = jmol;
		this.symm = symm;
	}

	@Override
	public void actionPerformed(ActionEvent ae) {
		String cmd = ae.getActionCommand();
		if (cmd.equals("New Symmetry Analysis"))
			SymmetryGui.getInstance();

		if (symm == null)
			logger.error("Currently not displaying a symmetry!");

		try {
			switch (cmd) {
				case "Repeats Superposition": {
					MultipleAlignmentJmol j = SymmetryDisplay.displayRepeats(symm);
					String s = SymmetryDisplay.printSymmetryAxes(symm, false);
					j.evalString(s);
					j.evalString("save STATE state_1");


					break;
				}
				case "Multiple Structure Alignment": {
					MultipleAlignmentJmol j = SymmetryDisplay.displayFull(symm);
					String s = SymmetryDisplay.printSymmetryAxes(symm);
					j.evalString(s);
					j.evalString("save STATE state_1");

					break;
				}
				case "Optimal Self Alignment": {
					Atom[] cloned = StructureTools.cloneAtomArray(symm.getAtoms());
					AbstractAlignmentJmol jmol = StructureAlignmentDisplay.display(
							symm.getSelfAlignment(), symm.getAtoms(), cloned);
					RotationAxis axis = new RotationAxis(symm.getSelfAlignment());
					jmol.evalString(axis.getJmolScript(symm.getAtoms()));
					jmol.setTitle(SymmetryDisplay.getSymmTitle(symm));

					break;
				}
				case "Show Symmetry Group":
					jmol.evalString(SymmetryDisplay.printSymmetryGroup(symm));
					break;
				case "Show Symmetry Axes": {
					jmol.evalString(SymmetryDisplay.printSymmetryAxes(symm));
					break;
				}
			}

		} catch (Exception e) {
			logger.error("Could not complete display option", e);
		}
	}
}
