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
package org.biojava.bio.structure.io;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import java.io.IOException;

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
	 * @deprecated Better to use {@link StructureIO#getBiologicalAssembly(String)}
	 * or a {@link BioUnitDataProvider}
	 */
	//public Structure getBiologicalUnit(String pdbId) throws StructureException, IOException;

	
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
