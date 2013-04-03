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
 * created at Sep 7, 2007
 */
package org.biojava.bio.structure.server;

import java.util.List;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;

/** an interface that defines different access methods to PDB installations.
 * Installations can be Flat file based ones, database ones, Installations that download all files from the net, etc.
 * 
 * 
 * @author Andreas Prlic
 * @deprecated
 */
public interface PDBInstallation{
	
	/** get all PDBHeaders that pass the added Filters, if no filters have been added returns all available PDBs
	 * 
	 * @return a list of PDBHeader objects
	 */
	public List<PDBHeader> getAll();
    
    /** get the PDB header for a single protein structure
     * 
     * @param pdbId
     * @return the PDB header object
     */
    public PDBHeader getPDBHeader(String pdbId);
	
	/** add a filter for PDB files. THis can be used to request, e.g. all X-ray structures, or all structures with
	 * a given resolution, all proteins with a certain function, etc.
	 * @see #getAll()
	 * @param filter the filter to apply when getAll is being called.
	 */
	public void addPDBFilter(PDBFilter filter);
	
	/** remove all filters, next time getAll is called, it will return all available PDBs
	 * 
	 */
	public void clearFilters();
	
	/** request a structure by its PDB identifier
	 * 
	 * @param pdbId
	 * @return the structure for the pdbId
	 */
	public Structure getStructure(String pdbId);

	/** iterate over all structures in this Installation that pass the provided filters and 
	 * return the next one in the list.
	 * @return the next structure 
	 */
	public Structure next();
	
	/** return if the iteration over all structures will return another structure
	 * 
	 * @return true if there is another structure that has not been iterated over yet
	 */
	public boolean hasNext();
	
	
}

