package org.biojava.nbio.structure.symmetry.gui;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.AlignmentCalculationRunnable;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** 
 * Calculates a symmetry analysis and displays the results.
 * Linked to the SymmetryGUI.
 * Does not generalize, uses CeSymm class directly to allow
 * for the symmetry axis recovery.
 *  
 * @author Aleix Lafita
 * @since 4.2.0
 * 
 */
public class SymmetryCalc implements AlignmentCalculationRunnable {
	
	private static final Logger logger = 
			LoggerFactory.getLogger(SymmetryCalc.class);
	
	boolean interrupted = false;
	
	private String name;
	private Structure structure;
	private SymmetryGui parent;

	/** Requests for a structure to analyze.
	 */
	public SymmetryCalc(SymmetryGui p, Structure s, String n) {
		parent = p;
		structure = s;
		name = n;
	}

	@Override
	public void run() {

		//The structure has been downloaded, now calculate the alignment ...
		CeSymm ceSymm = parent.getSymmetryAlgorithm();
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		
		try {

			Atom[] atoms = StructureTools.getRepresentativeAtomArray(structure);
			
			MultipleAlignment msa = ceSymm.analyze(atoms);

			List<String> names = new ArrayList<String>();
			for (int su=0; su<msa.size(); su++){
				names.add(name);
			}
			msa.getEnsemble().setStructureNames(names);

			MultipleAlignmentJmol jmol = 
					SymmetryDisplay.display(msa, ceSymm.getSymmetryAxes());
			String title = jmol.getTitle();
			
			if (params != null) 
				title += " | OrderDetector=" + params.getOrderDetectorMethod()+
				" Refiner: "+params.getRefineMethod();
			jmol.setTitle(title);

		} catch (StructureException e){
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
	public void setNrCPUs(int useNrCPUs) {}
}
