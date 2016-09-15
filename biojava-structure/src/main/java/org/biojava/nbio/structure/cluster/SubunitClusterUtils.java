package org.biojava.nbio.structure.cluster;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Some utility methods commonly used for Subunit clusters. They take as input a
 * List of {@link SubunitCluster}, the output of the {@link SubunitClusterer}.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class SubunitClusterUtils {

	/** Prevent instantiation **/
	private SubunitClusterUtils() {
	}

	/**
	 * Return a canonical String representation of the stoichiometry of a group
	 * of Subunits.
	 * <P>
	 * This method only uses alphabetic charaters for Subunit ids. For
	 * stoichiometries with a large number of entities, ? will be used when all
	 * alphabetic characters are used.
	 * 
	 * @param clusters
	 *            List of Subunit clusters
	 * @return String representation of the stoichiometry
	 */
	public static String getStoichiometryString(List<SubunitCluster> clusters) {

		// List number of members in each cluster
		List<Integer> stoichiometries = clusters.stream().map(c -> c.size())
				.collect(Collectors.toList());
		Collections.sort(stoichiometries);
		Collections.reverse(stoichiometries);

		// build formula string
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		StringBuilder formula = new StringBuilder();
		for (int i = 0; i < stoichiometries.size(); i++) {
			String key = "?";
			if (i < alpha.length())
				key = alpha.substring(i, i + 1);

			formula.append(key);
			if (stoichiometries.get(i) > 1)
				formula.append(stoichiometries.get(i));
		}

		return formula.toString();
	}

	/**
	 * A pseudostoichiometric {@link SubunitCluster} was obtained using the
	 * {@link SubunitClustererMethod#STRUCTURE} similarity.
	 * 
	 * @param clusters
	 * @return true if any of the clusters is pseudostoichiometric, false
	 *         otherwise
	 */
	public static boolean isPseudoStoichiometric(List<SubunitCluster> clusters) {

		for (SubunitCluster c : clusters) {
			if (c.getClustererMethod() == SubunitClustererMethod.STRUCTURE)
				return true;
		}
		return false;
	}

}
