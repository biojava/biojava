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
package org.biojava.nbio.structure.cluster;

/**
 * The SubunitClustererMethod ennummerates all methods that can be used to
 * cluster {@link Subunit} in the {@link SubunitCluster}.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public enum SubunitClustererMethod {

	/**
	 * The SEQUENCE clustering method uses the residue sequence of the
	 * {@link Subunit} to calculate sequence alignments.
	 * <p>
	 * Two {@link Subunit} with sufficient sequence identity and coverage are
	 * clustered together.
	 */
	SEQUENCE,

	/**
	 * The STRUCTURE clustering method uses the residue sequence and the
	 * coordinates of its Atom representatives of the {@link Subunit} to
	 * calculate sequence and structure alignments.
	 * <p>
	 * Two {@link Subunit} with sufficient structural similarity and coverage
	 * are clustered together.
	 */
	STRUCTURE,
	/**
	 * The SEQUENCE_STRUCTURE clustering method uses the residue sequence and the
	 * coordinates of its Atom representatives of the {@link Subunit} to
	 * calculate sequence and structure alignments.
	 * <p>
	 * Two {@link Subunit} with sufficient sequence identity and coverage are
	 * clustered together. Additionally, two {@link Subunit} with sufficient
	 * structural similarity and coverage are clustered together. If the
	 * sequence and structure clustering differ, the cluster contains
	 * pseudosymmetry (by definition).
	 */
	SEQUENCE_STRUCTURE
}

