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
package org.biojava.nbio.structure.align.multiple;

import java.util.Set;

import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;

/**
 * Interface for classes which implement a temporary cache for various numeric 
 * scores, e.g. RMSD, TM-Score, etc.
 * <p> 
 * Implementations can treat this as a cache, with as much or as little 
 * persistence as required. If getScore() returns null, callers should 
 * recalculate the score and then store it for future calls. Implementations
 * can therefore reset the cache if things which might impact the scores 
 * change.
 * <p>
 * A list of common scores are included to coordinate property names, but
 * additional scores can be stored as well. As a result, this is a flexible 
 * cache that can store any score as a key-value object.
 * 
 * @author Spencer Bliven
 * @since 4.1.0
 *
 */
public interface ScoresCache {

	/**
	 * Add a score to the list of scores.
	 * 
	 * @param property A string identifying the score and suitable for printing
	 * in headers. Example names found in: {@link MultipleAlignmentScorer}.
	 * @param score Value of the score
	 */
	public void putScore(String property, Double score);

	/**
	 * Get the value for a particular score. Scores which return null
	 * should be recalculated and then stored using 
	 * {@link #putScore(String, Double)}.
	 * 
	 * @param property Name of the score to fetch
	 * @return Value of the score, or null if it is not set.
	 */
	public Double getScore(String property);

	/**
	 * Get a collection of all scores that have been set.
	 * 
	 * @return Set of all score names
	 */
	public Set<String> getScores();
}
