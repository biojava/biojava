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
package org.biojava.nbio.structure.quaternary.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;

import java.util.List;

/** Provides access to the data that is needed in order to rebuild the correct biological assembly of a protein.
 * 
 * This is probably the simpler approach of accessing the necessary information. There is a second access layer, which is 
 * closer to the way the PDB is representing the files, it is defined by the interface RawBioUnitDataProvider.
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
	public List<BiologicalAssemblyTransformation>  getBioUnitTransformationList(String pdbId, int biolAssemblyNr);
	
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
	
	
	/** load the asym unit, but set the info how to re-create the bio unit in the PdbHeader object
	 * 
	 * @param pdbId
	 * @return
	 */
	public Structure getAsymUnit(String pdbId);
	
	public void setAsymUnit(Structure asymUnit);
	
	public void setAtomCache(AtomCache cache);
	
	public AtomCache getAtomCache();
	
}
