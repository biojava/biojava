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
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;

/**
 * A demo on how to use the quaternary symmetry detection algorithms.
 * 
 * @author Jose Duarte
 *
 */
public class DemoSymmetry {

	public static void main(String[] args) throws Exception {
		
		System.out.println("Getting all bioassemblies");
		List<Structure> bioAssemblies = StructureIO.getBiologicalAssemblies("4hhb");
		
		for (Structure bioAssembly : bioAssemblies) {			
			findQuatSym(bioAssembly);			
		}
		

	}

	private static void findQuatSym(Structure bioAssembly) throws StructureException {

		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();		
		
		System.out.println("GLOBAL SYMMETRY, NO CLUSTERING");
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		clusterParams.setSequenceIdentityThreshold(0.95);
		clusterParams.setRMSDThreshold(0.0);
		clusterParams.setClustererMethod(SubunitClustererMethod.SEQUENCE);

		QuatSymmetryResults globalResults = QuatSymmetryDetector.calcGlobalSymmetry(bioAssembly, symmParams, clusterParams);	
		
		
		
		System.out.println(globalResults.getSymmetry() + (globalResults.isPseudoStoichiometric()?"(pseudo)":""));
		
		System.out.println("There are "+globalResults.getSubunitClusters().size()+" subunit clusters");
		int i = 1;
		for (SubunitCluster suc : globalResults.getSubunitClusters()) {
			//System.out.println(suc.getClustererMethod());
			MultipleAlignment ma = suc.getMultipleAlignment();
						
			System.out.printf("Cluster %d (clustered by %s), RMSD = %4.2f\n", i, suc.getClustererMethod(), ma.getScore("RMSD"));
			
			i++;
		}

		System.out.println("\nGLOBAL SYMMETRY, WITH CLUSTERING (PSEUDO-SYMMETRY)");
		clusterParams = new SubunitClustererParameters();
		// we can either set sequence identity to 40% or rmsd to 2, both would have the same effect of clustering the alpha/beta subunits together
		clusterParams.setSequenceIdentityThreshold(0.4);
		clusterParams.setClustererMethod(SubunitClustererMethod.STRUCTURE);
		clusterParams.setRMSDThreshold(3.0);

		globalResults = QuatSymmetryDetector.calcGlobalSymmetry(bioAssembly, symmParams, clusterParams);



		System.out.println(globalResults.getSymmetry() + (globalResults.isPseudoStoichiometric()?"(pseudo)":""));

		System.out.println("There are "+globalResults.getSubunitClusters().size()+" subunit clusters");
		i = 1;
		for (SubunitCluster suc : globalResults.getSubunitClusters()) {
			//System.out.println(suc.getClustererMethod());
			MultipleAlignment ma = suc.getMultipleAlignment();

			System.out.printf("Cluster %d (clustered by %s), RMSD = %4.2f\n", i, suc.getClustererMethod(), ma.getScore("RMSD"));

			i++;
		}


		System.out.println("\n\nLOCAL SYMMETRIES");
		List<QuatSymmetryResults> localResults = QuatSymmetryDetector.calcLocalSymmetries(bioAssembly, symmParams, clusterParams);

		System.out.println("Number of local symmetries: "+localResults.size());

		for (QuatSymmetryResults results : localResults) {
			System.out.println(results.getSymmetry()+" "+results.getStoichiometry());
		}


	}
}
