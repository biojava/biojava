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
package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;

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
		 * dihedral: 4hhb, 1vym
		 * hierarchical: 4gcr, 1ppr.O, 1hiv
		 * monoclonal Ab: 4NZU
		 * 
		 * - For more examples see the symmetry benchmark
		 */

		//Set the name of the protein structure to analyze
		String name = "1u6d";

		//Download the atoms
		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure(name);
		Atom[] atoms = StructureTools.getRepresentativeAtomArray(s);

		CeSymm ceSymm = new CeSymm();

		//Choose some parameters
		CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
		params.setRefineMethod(RefineMethod.SINGLE);
		params.setSymmType(SymmetryType.AUTO);
		params.setOptimization(true);
		params.setSymmLevels(0);
		params.setSSEThreshold(2);

		//Run the alignment
		MultipleAlignment symmetry = ceSymm.analyze(atoms, params);
		
		//Display the results in FatCat format
		System.out.println(MultipleAlignmentWriter.toFatCat(symmetry));
		
		//Obtain the point group symmetry
		QuatSymmetryResults pg = SymmetryTools.getQuaternarySymmetry(symmetry);
		System.out.println("Point group internal symmetry: "+pg.getSymmetry());
	}
	
}
