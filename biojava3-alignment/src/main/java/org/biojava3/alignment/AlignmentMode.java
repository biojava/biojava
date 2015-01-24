package org.biojava3.alignment;

/**
 * Local, global, or semiglobal.
 * @author dmyersturnbull
 */
public enum AlignmentMode {

	GLOBAL,
	LOCAL,
	/**
	 * Global for the target sequence but local for the query.
	 * Great for database searches of a whole sequence, where the database is the query and the sequence being searched for is the target.
	 */
	SEMIGLOBAL

}
