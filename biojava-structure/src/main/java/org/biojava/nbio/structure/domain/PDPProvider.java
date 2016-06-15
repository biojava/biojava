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
 * Created on May 1, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.domain;

import java.io.IOException;
import java.util.SortedSet;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Decomposes a structure into representative PDP domains.
 *
 * Implementations will probably want to also implement {@link DomainProvider},
 * which provides a very similar set of methods for general structure domain
 * decomposition.
 * @author Andreas Prlic
 * @since 3.0.2
 */
public interface PDPProvider {

	/**
	 * Get a list of all PDP domains for a given PDB entry
	 * @param pdbId PDB ID
	 * @return Set of domain names, e.g. "PDP:4HHBAa"
	 * @throws IOException
	 */
	public SortedSet<String> getPDPDomainNamesForPDB(String pdbId) throws IOException;
	/**
	 * Get the structure for a particular PDP domain
	 * @param pdpDomainName PDP identifier, e.g. "PDP:4HHBAa"
	 * @param cache AtomCache, responsible for fetching and storing the coordinates
	 * @return Structure representing the PDP domain
	 * @throws IOException For IO errors, e.g. when parsing PDP information
	 * @throws StructureException For errors creating the structure
	 */
	public Structure getDomain(String pdpDomainName, AtomCache cache) throws IOException, StructureException;
	/**
	 * Get a StructureIdentifier representing the specified PDP domain.
	 *
	 * @param pdpDomainName PDP domain name
	 * @return a PDPDomain representing this domain name
	 * @throws IOException
	 */
	public PDPDomain getPDPDomain(String pdpDomainName) throws IOException;
}
