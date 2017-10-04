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
