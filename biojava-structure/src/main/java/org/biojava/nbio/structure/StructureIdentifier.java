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
 * Created on December 19, 2013
 * Author: Douglas Myers-Turnbull
 */

package org.biojava.nbio.structure;

import java.io.IOException;
import java.io.Serializable;

import org.biojava.nbio.structure.align.util.AtomCache;


/**
 * An identifier that <em>uniquely</em> identifies a whole {@link Structure} or
 * arbitrary substructure. Common examples would be reducing a structure to a
 * single chain, domain, or residue range.
 *
 * StructureIdentifiers are represented by unique strings. The getId() and fromId()
 * methods convert to and from the string representation.
 *
 * Implementations should provide a constructor which takes a String. A static
 * <tt>fromId(String)</tt> method is also recommended.
 *
 * @author dmyersturnbull
 * @author Spencer Bliven
 */
public interface StructureIdentifier extends Serializable {

	/**
	 * Get the String form of this identifier.
	 *
	 * It is recommended that the {@link #toString()} method also return the
	 * identifier, for consistency during serialization.
	 * @return The String form of this identifier
	 */
	String getIdentifier();


	/**
	 * Loads a structure encompassing the structure identified.
	 * The Structure returned should be suitable for passing as
	 * the input to {@link #reduce(Structure)}.
	 *
	 * It is recommended that the most complete structure available be returned
	 * (e.g. the full PDB) to allow processing of unselected portions where
	 * appropriate.
	 * @param AtomCache A potential sources of structures
	 * @return A Structure containing at least the atoms identified by this,
	 *  or null if Structures are not applicable.
	 * @throws StructureException For errors loading and parsing the structure
	 * @throws IOException Errors reading the structure from disk
	 */
	Structure loadStructure(AtomCache cache) throws StructureException, IOException;

	/**
	 * Convert to a canonical SubstructureIdentifier.
	 *
	 * <p>This allows all domains to be converted to a standard format String.
	 * @return A SubstructureIdentifier equivalent to this
	 * @throws StructureException Wraps exceptions that may be thrown by individual
	 *  implementations. For example, a SCOP identifier may require that the
	 *  domain definitions be available for download.
	 */
	SubstructureIdentifier toCanonical() throws StructureException;

	/**
	 * Takes a complete structure as input and reduces it to the substructure
	 * represented by this StructureIdentifier.
	 *
	 * <p>The returned structure may be a shallow copy of the input, with shared
	 * Chains, Residues, etc.
	 * @param input A full structure, e.g. as loaded from the PDB. The structure
	 * ID should match that returned by getPdbId(), if applicable.
	 * @return
	 * @throws StructureException
	 * @see StructureTools#getReducedStructure(Structure, String)
	 */
	Structure reduce(Structure input) throws StructureException;

}
