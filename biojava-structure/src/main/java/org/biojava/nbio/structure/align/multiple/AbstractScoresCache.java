package org.biojava.nbio.structure.align.multiple;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Abstact implementation of the {@link ScoresCache} with the shared code used
 * in all objects with a variables cache.
 * 
 * @author Spencer Bliven
 * @since 4.1.0
 *
 */
public abstract class AbstractScoresCache implements ScoresCache {
	
	private Map<String,Double> scores = null;

	protected AbstractScoresCache() {
		scores = null;
	}
	
	protected AbstractScoresCache(AbstractScoresCache cache) {
		this.scores = cache.scores;
	}

	@Override
	public void putScore(String property, Double score) {
		if(scores == null) {
			scores = new TreeMap<String, Double>();
		}
		scores.put(property, score);
	}

	@Override
	public Double getScore(String property) {
		if(scores != null && scores.containsKey(property)) {
			return scores.get(property);
		}
		return null;
	}

	@Override
	public Set<String> getScores() {
		if(scores == null) return Collections.emptySet();
		return Collections.unmodifiableSet(scores.keySet());
	}

	/**
	 * Subclasses should override clone and use the copy constructor.
	 * 
	 * @param e
	 * @return
	 * @throws CloneNotSupportedException
	 */
	protected Object clone(Object e) throws CloneNotSupportedException {
		throw new CloneNotSupportedException("Subclasses must override clone");
	}

	/**
	 * Clear the cached scores. This frees memory after the alignment changed.
	 */
	public void clear() {
		scores = null;
	}
}
