/*
 *                  BioJava development code
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
 * Created on Jul 16, 2006
 *
 */
package org.biojava.nbio.structure.align.gui;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.MultipleStructureAligner;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *  A class that obtains structures via DAS and aligns them.
 *  This is done in a separate thread.
 *  It is possible to register Event listeners to get notification of when
 *  the download has finished.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class MultipleAlignmentCalc implements AlignmentCalculationRunnable {

	private static final Logger logger =
			LoggerFactory.getLogger(MultipleAlignmentCalc.class);

	private List<StructureIdentifier> names;
	private List<Structure> structures;

	private MultipleAlignmentGUI parent;

	/**
	 * Requests an alignment of the pdbs.
	 * If they are empty strings, they are ignored.
	 *
	 * @param parent the gui frame that interacts with this class
	 * @param structures
	 * @param names
	 */
	public MultipleAlignmentCalc(MultipleAlignmentGUI parent,
			List<Structure> structures, List<StructureIdentifier> names) {

		this.parent= parent;
		this.structures = structures;
		this.names = names;
	}

	@Override
	public void run() {

		MultipleStructureAligner algorithm =
				parent.getMultipleStructureAligner();
		try {

			List<Atom[]> atomArrays = new ArrayList<Atom[]>();
			for (Structure s:structures){
				Atom[] ca = StructureTools.getRepresentativeAtomArray(s);
				atomArrays.add(ca);
			}

			MultipleAlignment msa = algorithm.align(atomArrays);
			msa.getEnsemble().setStructureIdentifiers(names);

			MultipleAlignmentJmolDisplay.display(msa);

		} catch (StructureException e) {
			e.printStackTrace();
			logger.warn(e.getMessage());
		}

		parent.notifyCalcFinished();
	}

	@Override
	public void interrupt() {}

	@Override
	public void cleanup() {

		parent.notifyCalcFinished();
		parent=null;
		structures = null;
		names = null;
	}

	@Override
	public void setNrCPUs(int useNrCPUs) {
		// TODO Auto-generated method stub
	}
}
