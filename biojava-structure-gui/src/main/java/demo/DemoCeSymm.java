package demo;

import java.io.IOException;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;

/**
 * Quick demo of how to call CE-Symm programmatically.
 * Some examples of different symmetry types are proposed.
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class DemoCeSymm {

	public static void main(String[] args) 
			throws IOException, StructureException {

		/* 
		 * Some examples:
		 * 
		 * CLOSED
		 * 2-fold: 1hiv.A, 
		 * 3-fold: 4i4q, 4dou
		 * 5-fold: 2jaj.A
		 * 6-fold: 1u6d
		 * 7-fold: 1jof.A
		 * 8-fold: 1vzw, d1i4na_
		 * 
		 * OPEN
		 * ankyrin: 1n0r.A, 3ehq.A
		 * leucine repeats: 2bnh.A, 3o6n
		 * helical: 1d0b.A
		 * 
		 * MULTIPLE AXES
		 * dihedral: 4hhb, 1vym, 1mmi, 1hiv
		 * hierarchical: 4gcr, 1ppr.O, 1hiv
		 * monoclonal Ab: 4NZU
		 * 
		 * - For more examples see the symmetry benchmark
		 */

		//Set the name of the protein structure to analyze
		String name = "1vym";

		//Download the atoms
		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure(name);
		Atom[] atoms = StructureTools.getRepresentativeAtomArray(s);

		CeSymm ceSymm = new CeSymm();

		//Choose some parameters
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setSymmetryType(SymmetryType.AUTO);
		params.setOptimization(true);
		params.setMultipleAxes(true);

		//Run the alignment
		MultipleAlignment symmetry = ceSymm.analyze(atoms, params);

		//Display the results in jmol
		SymmetryDisplay.display(symmetry, ceSymm.getSymmetryAxes());
	}
	
}
