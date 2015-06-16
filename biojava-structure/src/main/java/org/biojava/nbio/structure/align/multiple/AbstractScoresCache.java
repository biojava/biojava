package org.biojava.nbio.structure.align.multiple;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class AbstractScoresCache implements ScoresCache {
	private Map<String,Double> scores = null;
	
	protected AbstractScoresCache() {
		scores = null;
	}
	protected AbstractScoresCache(AbstractScoresCache cache) {
		this.scores = cache.scores;
	}
	/**
	 * Add a score to the list of scores.
	 * @param property A string identifying the score and suitable for printing in headers.
	 *  E.g. {@link #SCORE_TMSCORE}
	 * @param score Value of the score
	 * @param persistent Boolean indicating whether this score should be persisted
	 *  across serialization, cloning, and other changes.
	 */
	@Override
	public void putScore(String property, Double score) {
			if(scores == null) {
				scores = new HashMap<String, Double>();
			}
			scores.put(property, score);
	}
	/**
	 * Get the value for a particular score
	 * @param property Name of the score to fetch
	 * @return Value of the score, or null if it is not set.
	 */
	@Override
	public Double getScore(String property) {
		if(scores != null && scores.containsKey(property)) {
			return scores.get(property);
		}
		return null;
	}
	
	/**
	 * Get a collection of all scores which have been set
	 * @return Set of all score names
	 */
	@Override
	public Set<String> getScores() {
		return Collections.unmodifiableSet(scores.keySet());
	}
	
	/**
	 * Subclasses should override clone and use the copy constructor.
	 * @param e
	 * @return
	 * @throws CloneNotSupportedException
	 */
	protected Object clone(Object e) throws CloneNotSupportedException {
		throw new CloneNotSupportedException("Subclasses must override clone");
	}
	
	/**
	 * Resets all scores
	 */
	public void clear() {
		scores = null;
	}
}
