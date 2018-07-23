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
package org.biojava.nbio.structure.align.quaternary;

/**
 * The Quaternary Structure Relation describes the pairwise relation between two
 * quaternary structures.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public enum QsRelation {

	/**
	 * All the Subunits of one Structure have an equivalent in the other
	 * Structure.
	 * <p>
	 * An example of this relation is comparing an A2B2 complex with C2 symmetry
	 * against a similar A2B2 complex with the same C2 symmetry axis (same
	 * topology).
	 **/
	EQUIVALENT,

	/**
	 * All the Subunits of one Structure have an equivalent in the other
	 * Structure, but the other Structure contains additional non-matched
	 * Subunits.
	 * <p>
	 * An example of this relation is comparing an A2B2 complex with C2 symmetry
	 * against a similar A2 complex with the same C2 symmetry axis. Another
	 * example is comparing an A4 complex with D2 symmetry against an A2 complex
	 * with the same C2 symmetry axis as one of the D2 axes.
	 */
	PARTIAL_COMPLETE,

	/**
	 * Only a subset of Subunits of one Structure have an equivalent in the
	 * other Structure, and the other Structure also contains additional
	 * non-matched Subunits.
	 * <p>
	 * An example of this relation is comparing an A2B2 complex with C2 symmetry
	 * against an A2C2 complex, where As are similar but B and C differ, with
	 * the same C2 symmetry axis. Only the A2 part of the structure matches.
	 */
	PARTIAL_INCOMPLETE,

	/**
	 * None of the Subunits of one Structure have an equivalent in the other
	 * Structure.
	 * <p>
	 * Two quaternary structures are different if they are not
	 * {@link #EQUIVALENT}, {@link #PARTIAL_COMPLETE} or
	 * {@link #PARTIAL_INCOMPLETE}.
	 */
	DIFFERENT;
}
