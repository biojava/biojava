package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A utility object that describes Stoichiometry (composition of a protein assembly),
 * determined via clustering procedure {@link org.biojava.nbio.structure.cluster.SubunitClusterer},
 * and implements human-readable representation using various strategies.
 *
 * @author Dmytro Guzenko
 * @since 5.0.0
 */

public class Stoichiometry {

	/**
	 * What to do when the number of {@link SubunitCluster} exceeds the length of the alphabet.
	 */
	public enum StringOverflowStrategy {
		/**
		 * Put '?' symbol for every (alphabet.length+i)-th cluster
		 */
		QUESTIONMARK,
		/**
		 * Cycle through the alphabet (e.g., ...xyzABC...)
		 */
		CYCLE,
		/**
		 * Represent every cluster with two symbols from the alphabet,
		 * this forces us to specify number of subunits for every subunit (e.g., AA1AB1AC1...).
		 * This strategy will not work correctly if there are more than alphabet.length^2 subunit clusters.
		 */
		DOUBLE
	}

	/**
	 * Alphabet (a sequence of characters) used in this stoichiometry to construct human-readable representation.
	 */
	private String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

	/**
	 * Strategy determines how this stoichiometry will construct human-readable representation in case number
	 * of clusters exceeds number of letters in the alphabet.
	 */
	private StringOverflowStrategy strategy = StringOverflowStrategy.CYCLE;

	/**
	 * Subunit clusters that define this stoichiometry.
	 */
	private List<SubunitCluster> orderedClusters = new ArrayList<>();

	/** Prevent instantiation **/
	private Stoichiometry() {
	}

	/**
	 * Constructor for Stoichiometry. The default strategy is CYCLE,
	 * the letters assigned for each cluster will be reset.
	 *
	 * @param clusters
	 *            List of {@link SubunitCluster} that defines assembly composition.
	 */
	public Stoichiometry(List<SubunitCluster> clusters) {
		this(clusters,StringOverflowStrategy.CYCLE,true);
	}

	/**
	 * Constructor for Stoichiometry. The default strategy is CYCLE.
	 *
	 * @param clusters
	 *            List of {@link SubunitCluster} that defines assembly composition.
	 * @param resetAlphas
	 *            Whether to keep alphas assigned to {@link SubunitCluster} object (useful for local symmetry detection)
	 *            or to generate them anew.
	 */
	public Stoichiometry(List<SubunitCluster> clusters, boolean resetAlphas) {
		this(clusters,StringOverflowStrategy.CYCLE,resetAlphas);
	}

	/**
	 * Constructor for Stoichiometry. The alphas assigned to {@link SubunitCluster} objects will be reset.
	 *
	 * @param clusters
	 *            List of {@link SubunitCluster} that defines assembly composition.
	 * @param strategy
	 *            What to do if number of {@link SubunitCluster} exceeds the alphabet length.
	 */
	public Stoichiometry(List<SubunitCluster> clusters, StringOverflowStrategy strategy) {
		this(clusters,strategy,true);
	}

	/**
	 * Constructor for Stoichiometry.
	 *
	 * @param clusters
	 *            List of {@link SubunitCluster} that defines assembly composition.
	 * @param strategy
	 *            What to do if number of {@link SubunitCluster} exceeds the alphabet length.
	 * @param resetAlphas
	 *            Whether to keep alphas assigned to {@link SubunitCluster} object (useful for local symmetry detection)
	 *            or to generate them anew.
	 */
	public Stoichiometry(List<SubunitCluster> clusters, StringOverflowStrategy strategy, boolean resetAlphas) {
		this.strategy = strategy;
		this.orderedClusters =
				clusters.stream().
					sorted(Comparator.
						comparing(SubunitCluster::size).
						reversed()).
					collect(Collectors.toList());
		if (resetAlphas) {
			doResetAlphas();
		}
	}

	private void doResetAlphas() {
		for (int i = 0; i < this.orderedClusters.size(); i++) {
			this.orderedClusters.get(i).setAlpha(generateAlpha(i));
		}
	}

	/**
	 * Produce a string ("alpha") that describes each component depending on the current strategy.
	 * @param clusterInd
	 *          component index
	 * @return alphanumeric string.
	 */
	private String generateAlpha(int clusterInd) {
		String key;
		int alphabetInd;
		switch (strategy) {
			case CYCLE:
				alphabetInd = clusterInd % alphabet.length();
				key = alphabet.substring(alphabetInd, alphabetInd + 1);
				break;

			case DOUBLE:
				if (orderedClusters.size()>alphabet.length()) {
					int alphabetInd1 = clusterInd / alphabet.length();
					int alphabetInd2 = clusterInd % alphabet.length();
					key = alphabet.substring(alphabetInd1, alphabetInd1 + 1);
					key+=alphabet.substring(alphabetInd2, alphabetInd2 + 1);
				} else {
					key = alphabet.substring(clusterInd, clusterInd + 1);
				}
				break;

			case QUESTIONMARK:
				key = "?";
				if(clusterInd<alphabet.length()) {
					key = alphabet.substring(clusterInd, clusterInd + 1);
				}
				break;

			default:
				key = "?";
				if(clusterInd<alphabet.length()) {
					key = alphabet.substring(clusterInd, clusterInd + 1);
				}
				break;
		}
		return key;
	}
	/**
	 * @return list of {@link SubunitCluster}, ordered by the number of subunits (decreasing).
	 */
	public List<SubunitCluster> getClusters() {
		return orderedClusters;
	}

	/**
	 * @return Number of distinct components in this stoichiometry.
	 */
	public int numberOfComponents() {
		return orderedClusters.size();
	}

	/**
	 * Make a combined Stoichiometry object of <i>this</> and the <i>other</>.
	 * All clusters are assumed to be distinct, so the alphas will be reset in case of repeats.
	 * The combined list of clusters will be ordered by the number of subunits.
	 * @return new {@link Stoichiometry} object.
	 */
	public Stoichiometry combineWith(Stoichiometry other) {
		List<SubunitCluster> combinedClusters = new ArrayList<>();
		combinedClusters.addAll(this.orderedClusters);
		combinedClusters.addAll(other.orderedClusters);

		// check that there is no alpha duplication
		Set<String> thisAlpha = this.orderedClusters.stream().map(SubunitCluster::getAlpha).collect(Collectors.toSet());
		Set<String> otherAlpha = other.orderedClusters.stream().map(SubunitCluster::getAlpha).collect(Collectors.toSet());
		boolean resetAlphas = true;
		if(Collections.disjoint(thisAlpha,otherAlpha)) {
			resetAlphas = false;
		}
		Stoichiometry combinedStoichiometry = new Stoichiometry(combinedClusters,this.strategy,resetAlphas);

		return combinedStoichiometry;
	}

	/**
	 * Make a Stoichiometry object that corresponds to a single component.
	 * @param i component index
	 * @return new {@link Stoichiometry} object.
	 */
	public Stoichiometry getComponent(int i) {
		return new Stoichiometry(Collections.singletonList(orderedClusters.get(i)),this.strategy,false);
	}

	/**
	 * @return {@link StringOverflowStrategy} used in this stoichiometry
	 *          to construct human-readable representation in case number
	 *          of clusters exceeds number of letters in the alphabet.
	 */
	public StringOverflowStrategy getStrategy() {
		return strategy;
	}

	/**
	 * Change string representation of a stoichiometry in case number of clusters exceeds number of letters in the alphabet.
	 * This action may invalidate alphas already assigned to the clusters.
	 * @param strategy
	 *          {@link StringOverflowStrategy} used in this stoichiometry
	 *          to construct human-readable representation in case number
	 *          of clusters exceeds number of letters in the alphabet.
	 */
	public void setStrategy(StringOverflowStrategy strategy) {
		if(this.strategy != strategy) {
			this.strategy = strategy;
			if(orderedClusters.size()>alphabet.length())
				doResetAlphas();
		}
	}

	/**
	 * @return Alphabet (a sequence of characters) used in this stoichiometry to construct human-readable representation.
	 */
	public String getAlphabet() {
		return alphabet;
	}

	/**
	 * Change alphabet used for string representation of a stoichiometry.
	 * This action invalidates alphas already assigned to the clusters.
	 * @param alphabet
	 *          a sequence of characters used in this stoichiometry to construct human-readable representation.
	 */
	public void setAlphabet(String alphabet) {
		this.alphabet = alphabet;
		doResetAlphas();
	}

	/**
	 * @return Human-readable representation of this stoichiometry.
	 */
	@Override
	public String toString() {
		StringBuilder formula = new StringBuilder();

		orderedClusters.forEach((SubunitCluster r) -> {
			formula.append(r.getAlpha());
			if(r.getAlpha().length()>1 || r.size()>1)
				formula.append(r.size());
		});

		return formula.toString();
	}

	/**
	 * A pseudostoichiometric {@link SubunitCluster} was obtained using the
	 * {@link SubunitClustererMethod#STRUCTURE} similarity,
	 * or {@link SubunitClustererMethod#SEQUENCE} similarity with low scores.
	 *
	 * @return true if any of the clusters is pseudostoichiometric, false
	 *         otherwise
	 */
	public boolean isPseudoStoichiometric() {
		for (SubunitCluster c : orderedClusters) {
			if(c.isPseudoStoichiometric())
				return true;
		}
		return false;
	}

}
