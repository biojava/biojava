package org.biojava.nbio.structure.align.multiple;

import java.util.Set;

/**
 * Interface for classes which implement a temporary cache for various numeric scores,
 * e.g. RMSD, TM-Score, etc.
 * 
 * Implementations can treat this as a cache, with as much or as little persistence
 * as required. If getScore() returns null, callers should recalculate the score
 * and then store it for future calls. Implementations can therefor reset the
 * cache if things change which might impact the scores.
 * 
 * A list of common scores are included to coordinate property names, but
 * additional scores can be stored as well.
 * 
 * @author Spencer Bliven
 *
 */
public interface ScoresCache {

	/**
	 * Add a score to the list of scores.
	 * @param property A string identifying the score and suitable for printing in headers.
	 *  E.g. {@link #SCORE_TMSCORE}
	 * @param score Value of the score
	 */
	public void putScore(String property, Double score);
	/**
	 * Get the value for a particular score. Scores which return null
	 * should be recalculated and then stored using {@link #putScore(String, Double)};
	 * @param property Name of the score to fetch
	 * @return Value of the score, or null if it is not set.
	 */
	public Double getScore(String property);
	
	/**
	 * Get a collection of all scores that have been set
	 * @return Set of all score names
	 */
	public Set<String> getScores();
}
