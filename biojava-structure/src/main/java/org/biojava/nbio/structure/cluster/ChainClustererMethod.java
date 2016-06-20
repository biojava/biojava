package org.biojava.nbio.structure.cluster;

public enum ChainClustererMethod {

	/**
	 * The SEQUENCE clustering method uses the SEQRES of the Chains in the
	 * Structure to calculate a sequence alignment.
	 * <p>
	 * Chains with sufficient sequence identity and coverage are clustered
	 * together.
	 */
	SEQUENCE,

	/**
	 * The STRUCTURE clustering method uses the SEQRES and the representative
	 * Atoms of the residues of the Chains in the Structure to calculate a
	 * sequence and a structure alignment.
	 * <p>
	 * Chains with sufficient sequence identity and coverage are clustered
	 * together. Additionally, chains with sufficient structural similarity and
	 * coverage are clustered together. If the sequence and structure clustering
	 * differ, the structure contains pseudosymmetry (by definition).
	 */
	STRUCTURE,

	/**
	 * The INTERNAL SYMMETRY clustering method uses the SEQRES and the
	 * representative Atoms of the residues of the Chains in the Structure to
	 * calculate a sequence and a structure alignment.
	 * <p>
	 * Chains with sufficient sequence identity and coverage are clustered
	 * together. Additionally, chains with sufficient structural similarity and
	 * coverage are clustered together. If the sequence and structure clustering
	 * differ, the structure contains pseudosymmetry (by definition).
	 * <p>
	 * In a final step, internal symmetry of each chain cluster is analyzed. If
	 * the chain cluster is internally symmetric, chains are divided into its
	 * internally symmetric domains.
	 */
	INTERNAL_SYMMETRY;
}
