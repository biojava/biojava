package org.biojava.bio.structure.quaternary.io;

import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;

/** Provides access to the data that is needed in order to rebuild the correct biological assembly of a protein.
 * 
 * This is probably the simpler approach of accessing the necessary information. There is a second access layer, which is 
 * closer to the way the PDB is representing the files.
 * 
 * @author Andreas Prlic
 *
 */
public interface BioUnitDataProvider {
	
	/** get the data for a particular assembly, counting starts at 1...
	 * 
	 * @param pdbId the PDB ID. E.g. 1STP
	 * @param biolAssemblyNr the number of the assembly, the first one is nr 1. 0 refers to the asym unit
	 * @return list of transformations.
	 */
	public List<ModelTransformationMatrix>  getBioUnitTransformationList(String pdbId, int biolAssemblyNr);
	
	/** Returns the number of available biological assemblies.
	 *  @param pdbId the PDB ID. E.g. 1STP
	 * @return nr of available assemblies.
	 */
	public int getNrBiolAssemblies(String pdbId);
	
	
	/** Does the PDB ID have biological assembly information?
	 * 
	 * @param pdbId the PDB ID. E.g. 1STP
	 * @return boolean flag
	 */
	public boolean hasBiolAssembly(String pdbId);
}
