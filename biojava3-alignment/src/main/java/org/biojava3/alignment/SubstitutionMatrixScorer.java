package org.biojava3.alignment;

import org.biojava3.alignment.template.AbstractScorer;
import org.biojava3.alignment.template.PairwiseSequenceScorer;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Scores using a substitution matrix. Specifically, the score is the sum of the substitution matrix entries
 * corresponding to the alignment. Gaps are scored according to the substitution matrix, just as matches and mismatches.
 * @author dmyersturnbull
 *
 * @param <S>
 * @param <C>
 */
public class SubstitutionMatrixScorer<S extends Sequence<C>, C extends Compound> extends AbstractScorer
implements PairwiseSequenceScorer<S, C> {
	
	private final SubstitutionMatrix<C> matrix;
	
	private S query;
	private S target;
	private double score;
	
	public SubstitutionMatrixScorer(SequencePair<S, C> pair, SubstitutionMatrix<C> matrix) {
		super();
		this.query = pair.getQuery().getOriginalSequence();
		this.target = pair.getTarget().getOriginalSequence();
		this.matrix = matrix;
		for (int i = 1; i <= pair.getLength(); i++) {
			C query = pair.getCompoundAt(1, i);
			C target = pair.getCompoundAt(2, i);
			score += matrix.getValue(query, target);
		}
	}

	/**
	 * @return The maximum score the query could be assigned when aligned against any target sequence.
	 */
	@Override
	public double getMaxScore() {
		// assume nothing about the matrix
		double score = 0;
		for (C queryC : query.getAsList()) {
			short max = Short.MIN_VALUE;
			for (Short value : matrix.getRow(queryC).values()) {
				if (value > max) max = value;
			}
			score += max;
		}
		return score;
	}

	/**
	 * @return The minimum score the query could be assigned when aligned against any target sequence.
	 */
	@Override
	public double getMinScore() {
		// assume nothing about the matrix
		double score = 0;
		for (C queryC : query.getAsList()) {
			short min = Short.MAX_VALUE;
			for (Short value : matrix.getRow(queryC).values()) {
				if (value < min) min = value;
			}
			score += min;
		}
		return score;
	}

	@Override
	public double getScore() {
		return score;
	}

	@Override
	public S getQuery() {
		return query;
	}

	@Override
	public S getTarget() {
		return target;
	}
	
}
