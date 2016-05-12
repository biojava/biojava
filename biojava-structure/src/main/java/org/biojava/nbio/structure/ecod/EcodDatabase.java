/*
 * BioJava development code
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
 */

package org.biojava.nbio.structure.ecod;

import java.io.IOException;
import java.util.List;

/** General API for interacting with ECOD.
 *
 * @author Spencer Bliven
 */
public interface EcodDatabase {

	/** Return the release version.
	 *
	 * @return version
	 * @throws IOException
	 */
	public String getVersion() throws IOException;

	/**
	 * Get a particular ECOD domain by the domain ID (e.g. "e4hhbA1")
	 * @param ecodId
	 * @return
	 * @throws IOException
	 */
	public EcodDomain getDomainsById(String ecodId) throws IOException;

	/**
	 * Get a list of all ECOD domains for a particular PDB ID
	 * @param pdbId
	 * @return the list of domains, or null if no matching domains were found
	 * @throws IOException
	 */
	public List<EcodDomain> getDomainsForPdb(String pdbId) throws IOException;

	/**
	 * Get a list of domains within a particular level of the hierarchy
	 * @param hierarchy A dot-separated list giving the X-group, H-group, and/or
	 *  T-group (e.g. "1.1" for all members of the RIFT-related H-group)
	 * @return
	 * @throws IOException
	 */
	public List<EcodDomain> filterByHierarchy(String hierarchy) throws IOException;

	/**
	 * Get all ECOD domains
	 * @return
	 * @throws IOException
	 */
	public List<EcodDomain> getAllDomains() throws IOException;
}
