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
package org.biojava.nbio.structure.domain;

import java.io.IOException;
import java.util.SortedSet;

import org.biojava.nbio.structure.StructureException;

/**
 * Decomposes a structure from the PDB into representative domains
 */
public interface DomainProvider {

	/**
	 * Get a list of constituent domain identifiers
	 * @param name a structure identifier
	 * @return A list of domain names
	 * @throws IOException For IO errors getting the domains
	 * @throws StructureException For errors converting name to a valid identifier
	 */
	public SortedSet<String> getDomainNames(String name) throws IOException, StructureException;
	/**
	 * Get the full list of representative domains for the PDB.
	 *
	 * The exact definition representatives is implementation-specific, but
	 * should cover as many structures as possible.
	 * @return A full list of all representative domains recognized by this provider
	 * @throws IOException For IO errors getting the representatives
	 */
	public SortedSet<String> getRepresentativeDomains() throws IOException;
}
