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
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Action Listener for the symmetry menu.
 * Trigger various symmetry analysis.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryListener implements ActionListener{

	private MultipleAlignmentJmol jmol;
	private MultipleAlignment msa;
	private SymmetryAxes axes;
	
	private static final Logger logger = 
			LoggerFactory.getLogger(SymmetryListener.class);

	public SymmetryListener(MultipleAlignmentJmol jmol, SymmetryAxes axes) {
		this.jmol = jmol;
		if (jmol != null) this.msa = jmol.getMultipleAlignment();
		this.axes = axes;
	}

	@Override
	public void actionPerformed(ActionEvent ae) {

		String cmd = ae.getActionCommand();
		if (cmd.equals("Subunit Superposition")){
			if (msa == null) {
				logger.error("Currently not displaying a symmetry!");
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
				logger.error("Currently not displaying a symmetry!");
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
			if (msa == null) {
				logger.error("Currently not displaying a symmetry!");
				return;
			}
			String script = SymmetryDisplay.printPointGroupAxes(msa);
			jmol.evalString(script);
			return;
			
		} else if (cmd.equals("Show Symmetry Axes")){
			if (axes != null) {
				String s = SymmetryDisplay.printSymmetryAxes(msa, axes, false);
				jmol.evalString(s);
				return;
			} else logger.error("No axes found for this symmetry result");
			
		} else if (cmd.equals("New Symmetry Analysis")){
			SymmetryGui.getInstance();
		}
	}
	
}