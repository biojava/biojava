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

import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList;

import java.util.List;


/** Defines the methods that have to be implemented by a class that provides the data that is necessary to recreate the correct biological assembly of a protein.
 * This is very close to the way PDB is representing biological assemblies. For the outside it is probably easier to use the other way of accessing the data as defined in 
 * BioUnitDataProvider
 * 
 * @author Andreas Prlic
 * @since 3.0.5
 */
public interface RawBioUnitDataProvider {
	
	/** Tell the provider for which PDB ID the quaternary structure should be returned.
	 * 
	 * @param pdbId
	 */
	public void setPdbId(String pdbId);
	
	/** Data access method for list describing all assemblies
	 * 
	 * @return
	 */
	public List<PdbxStructAssembly> getPdbxStructAssemblies();
	
	/** Data access method for list describing all assemblies
	 * 
	 * @return
	 */
	public List<PdbxStructAssemblyGen> getPdbxStructAssemblyGens();
	
	/** Get all the possible operators
	 * 
	 * @return
	 */
	public List<PdbxStructOperList> getPdbxStructOperList();
	
	
	/** Returns the number of available biological assemblies.
	 * 
	 * @return
	 */
	public int getNrBiolAssemblies();
	
	
	/** Does the PDB ID have biological assembly information?
	 * 
	 * @return boolean flag
	 */
	public boolean hasBiolAssembly();
	
	/** get the data for a particular pdbxStructAssembly. We start counting at 0.
	 * 
	 * @param biolAssemblyNr
	 * @return
	 */
	public PdbxStructAssembly getPdbxStructAssembly(int biolAssemblyNr);
	
	
	/** get the data for a particular pdbxStructAssemblyGen. We start counting at 0.
	 * 
	 * @param biolAssemblyNr
	 * @return
	 */
	public List<PdbxStructAssemblyGen> getPdbxStructAssemblyGen(int biolAssemblyNr);
	
	
}
