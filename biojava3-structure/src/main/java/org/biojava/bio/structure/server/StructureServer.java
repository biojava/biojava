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


/** a server that pre-loads PDB files in memory, for quicker processing by StructureListeners
 * 
 * @author Andreas Prlic
 * @deprecated
 */
public interface StructureServer {

	public void addStructureListener(StructureListener listener);
	public void clearStructureListeners();
	public void setPDBInstallation(PDBInstallation installation);
	public PDBInstallation getPDBInstallation();
	
	/** set how many structures should be kept in memory - for quicker access
	 * 
	 * @param nr
	 */
	public void setCacheSize(int  nr);
	public int getCacheSize();
	
	public void requestNextStructure(StructureListener listener); 
	public boolean hasNextStructure();
}
