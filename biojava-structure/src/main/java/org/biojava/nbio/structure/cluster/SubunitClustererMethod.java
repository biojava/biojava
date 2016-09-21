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
	 * The IDENTITY clustering method uses the residue sequence of the
	 * {@link Subunit}. Two {@link Subunit} with exactly the same sequence will
	 * be clustered together.
	 */
	IDENTITY,

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
	 * Two {@link Subunit} with sufficient sequence identity and coverage are
	 * clustered together. Additionally, two {@link Subunit} with sufficient
	 * structural similarity and coverage are clustered together. If the
	 * sequence and structure clustering differ, the cluster contains
	 * pseudosymmetry (by definition).
	 */
	STRUCTURE;

}
