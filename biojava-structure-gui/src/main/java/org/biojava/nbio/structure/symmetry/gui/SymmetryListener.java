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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Action Listener for the symmetry menu. Trigger various symmetry analysis.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryListener implements ActionListener {

	private MultipleAlignmentJmol jmol;
	private CeSymmResult symm;

	private static final Logger logger = LoggerFactory
			.getLogger(SymmetryListener.class);

	public SymmetryListener(MultipleAlignmentJmol jmol, CeSymmResult symm) {
		this.jmol = jmol;
		this.symm = symm;
	}

	@Override
	public void actionPerformed(ActionEvent ae) {

		String cmd = ae.getActionCommand();
		if (cmd.equals("Subunit Superposition")) {
			if (symm == null) {
				logger.error("Currently not displaying a symmetry!");
				return;
			}
			try {
				MultipleAlignmentJmol j = SymmetryDisplay.displaySubunits(symm
						.getMultipleAlignment());
				String s = SymmetryDisplay.printSymmetryAxes(symm, true);
				j.evalString(s);
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Multiple Structure Alignment")) {
			if (symm == null) {
				logger.error("Currently not displaying a symmetry!");
				return;
			}
			try {
				MultipleAlignmentJmol j = SymmetryDisplay.displayFull(symm
						.getMultipleAlignment());
				String s = SymmetryDisplay.printSymmetryAxes(symm, false);
				j.evalString(s);
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Optimal Self-Alignment")) {
			if (symm == null) {
				logger.error("Currently not displaying a symmetry!");
				return;
			}
			try {
				StructureAlignmentDisplay.display(symm.getSelfAlignment(),
						symm.getAtoms(), symm.getAtoms());
			} catch (StructureException e) {
				e.printStackTrace();
			}

		} else if (cmd.equals("Point Group Symmetry")) {
			if (symm == null) {
				logger.error("Currently not displaying a symmetry!");
				return;
			}
			String script = SymmetryDisplay.printPointGroupAxes(symm);
			jmol.evalString(script);
			return;

		} else if (cmd.equals("Show Symmetry Axes")) {
			if (symm != null) {
				String s = SymmetryDisplay.printSymmetryAxes(symm, false);
				jmol.evalString(s);
				return;
			} else
				logger.error("Currently not displaying a symmetry!");

		} else if (cmd.equals("New Symmetry Analysis")) {
			SymmetryGui.getInstance();
		}
	}
}