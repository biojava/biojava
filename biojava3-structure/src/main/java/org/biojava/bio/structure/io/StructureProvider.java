package org.biojava.bio.structure.io;

import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

/** A class that can provide a protein structure object from somewhere.
 * 
 * @author Andreas Prlic
 * @since 3.2
 */
public interface StructureProvider {

	/** get the structure for a PDB ID
	 * 
	 * @param pdbId
	 * @return
	 */
	public Structure getStructureById(String pdbId) throws StructureException,IOException;
	
	/** get the biological unit for a file
	 * 
	 * @param pdbId
	 * @return
	 */
	public Structure getBiologicalUnit(String pdbId) throws StructureException, IOException;

	
	/** Set the parameters that should be used for file parsing
	 * 
	 * @param params FileParsingParameters
	 */
	public void setFileParsingParameters(FileParsingParameters params);
	   
	
	/** Get the parameters that should be used for file parsing
     * 
     * @return the FileParsingParameters that are configuring the behavior of the parser
     */
    public FileParsingParameters getFileParsingParameters();
    
    
}
