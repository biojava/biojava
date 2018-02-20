package org.biojava.nbio.structure.cluster;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Some utility methods commonly used for Subunit clusters. They take as input a
 * List of {@link SubunitCluster}, the output of the {@link SubunitClusterer}.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 * @deprecated Use {@link org.biojava.nbio.structure.symmetry.core.Stoichiometry} instead.
 */
@Deprecated
public class SubunitClusterUtils {

	/** Prevent instantiation **/
	private SubunitClusterUtils() {
	}

	/**
	 * Return a canonical String representation of the stoichiometry of a group
	 * of Subunits.
	 * <P>
	 * This method only uses alphabetic charaters for Subunit ids.
	 * For stoichiometries with a large number of entities,
	 * when all alphabetic characters are used, they will restart from the beginning.
	 *
	 * @param clusters
	 *            List of Subunit clusters
	 * @return String representation of the stoichiometry
	 */
	public static String getStoichiometryString(List<SubunitCluster> clusters) {

		// List number of members in each cluster
		List<Integer> stoichiometries =
			clusters.stream().
				map(SubunitCluster::size).
				sorted().
				collect(Collectors.toList());

		Collections.reverse(stoichiometries);

		// build formula string
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		StringBuilder formula = new StringBuilder();
		for (int i = 0; i < stoichiometries.size(); i++) {
			String key = clusters.get(i).getAlpha();
			if (key.isEmpty()) {
				int ind = i % alpha.length();
				key = alpha.substring(ind, ind + 1);
			}

			formula.append(key);
			if (stoichiometries.get(i) > 1)
				formula.append(stoichiometries.get(i));
		}

		return formula.toString();
	}


	/**
	 * Order clusters by the number of subunits (decreasing) and assign alphabetic characters in order.
	 * <P>
	 * This method makes an alpha-identifier a property of a cluster, which is useful for local symmetries,
	 * whose stoichiometries will be consistent with respect to the global clustering.
	 *
	 * @param clusters
	 *            List of Subunit clusters
	 * @return Ordered list of Subunit clusters with alpha-identifiers.
	 */
	public static List<SubunitCluster> orderByStoichiometry(List<SubunitCluster> clusters) {
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

		List<SubunitCluster> clustersSorted =
			clusters.stream().
				sorted(Comparator.
					comparing(SubunitCluster::size).
					reversed()).
				collect(Collectors.toList());

		for (int i = 0; i < clustersSorted.size(); i++) {
			int ind = i % alpha.length();
			String key = alpha.substring(ind, ind + 1);
			clustersSorted.get(i).setAlpha(key);
		}
		return clustersSorted;
	}

	/**
	 * A pseudostoichiometric {@link SubunitCluster} was obtained using the
	 * {@link SubunitClustererMethod#STRUCTURE} similarity,
	 * or {@link SubunitClustererMethod#SEQUENCE} similarity with low scores.
	 *
	 * @param clusters
	 * @return true if any of the clusters is pseudostoichiometric, false
	 *         otherwise
	 */
	public static boolean isPseudoStoichiometric(List<SubunitCluster> clusters) {

		for (SubunitCluster c : clusters) {
			if(c.isPseudoStoichiometric())
				return true;
		}
		return false;
	}

}
