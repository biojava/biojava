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

import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;

public class DemoBioAssemblies {

	public static void main(String[] args) throws Exception {
		
		// 1st method: get 1 bioassembly at a time, parses the file each time
		System.out.println("Getting one bioassembly at a time");
		Structure asymUnit = StructureIO.getStructure("2trx");
		System.out.println("Number of bioassemblies: "+asymUnit.getPDBHeader().getNrBioAssemblies());

		for (int id = 1; id<=asymUnit.getPDBHeader().getNrBioAssemblies(); id++) {
			Structure bioAssembly = StructureIO.getBiologicalAssembly("2trx", id);
			findQuatSym(bioAssembly);
		}

		
		// 2nd method: get all bioassemblies at once, parses the file only once
		System.out.println("Getting all bioassemblies");
		List<Structure> bioAssemblies = StructureIO.getBiologicalAssemblies("2trx");
		
		for (Structure bioAssembly : bioAssemblies) {			
			findQuatSym(bioAssembly);			
		}
		

	}

	private static void findQuatSym(Structure bioAssembly) {
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				bioAssembly, symmParams, clusterParams);

		// C2 symmetry non pseudosymmetric
		System.out.println(symmetry.getSymmetry() +" "+ symmetry.getStoichiometry());

	}
}
