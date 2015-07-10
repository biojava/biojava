package org.biojava.nbio.ws.uniprot;

/**
 * Just a factory wrapper to prevent users from directly accessing the 
 * exact implementations object.
 * @author pbansal
 */
public class UniprotUtilities {
	public static BatchRetrieve getBatchUtility() {
		return new BatchRetrieveUniprotEntries(); 
	}
} // UniprotUtilities
