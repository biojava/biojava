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
package org.biojava.nbio.structure;

import java.io.Serializable;

/** 
 * An interface implemented by all classes that represent PDB records.
 *
 * @author Andreas Prlic
 * @since 1.6
 */
public interface PDBRecord extends Serializable {

	/** Returns a PDB file like representation of this record.
	 *
	 * @return a String providing a PDB file like representation of the record.
	 */
	public String toPDB();


	/** Appends a PDB file like representation of this record to the provided StringBuffer.
	 *
	 */
	public void toPDB(StringBuffer buf);


}
