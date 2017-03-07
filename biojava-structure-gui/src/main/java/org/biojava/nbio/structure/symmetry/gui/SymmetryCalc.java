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
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.AlignmentCalculationRunnable;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Calculates a symmetry analysis and displays the results. Linked to the
 * SymmetryGUI. Does not generalize, uses CeSymm class directly to allow for the
 * symmetry axis recovery.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryCalc implements AlignmentCalculationRunnable {

	private static final Logger logger = LoggerFactory
			.getLogger(SymmetryCalc.class);

	boolean interrupted = false;

	private Structure structure;
	private SymmetryGui parent;

	/**
	 * Requests for a structure to analyze.
	 */
	public SymmetryCalc(SymmetryGui p, Structure s) {
		parent = p;
		structure = s;
	}

	@Override
	public void run() {

		CESymmParameters params = parent.getParameters();

		try {

			Atom[] atoms = SymmetryTools.getRepresentativeAtoms(structure);

			CeSymmResult result = CeSymm.analyze(atoms, params);
			SymmetryDisplay.display(result);

		} catch (StructureException e) {
			logger.warn(e.getMessage());
		}
		parent.notifyCalcFinished();
	}

	@Override
	public void interrupt() {
		interrupted = true;
	}

	@Override
	public void cleanup() {

		parent.notifyCalcFinished();
		parent = null;
		structure = null;
	}

	@Override
	public void setNrCPUs(int useNrCPUs) {
	}
}
