package org.biojava.bio.structure;

/** An interface implemented by all classes that represent PDB records
 *
 * @author Andreas Prlic
 * @since 1.6
 */
public interface PDBRecord {

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
