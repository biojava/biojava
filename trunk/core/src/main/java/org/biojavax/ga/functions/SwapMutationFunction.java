package org.biojavax.ga.functions;

import java.util.Random;

import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.functions.AbstractMutationFunction;

/**
 * This class does a sort of mutation by exchanging two positions on the
 * chromosome. Thus it can be used for implementations where a change of the
 * amount of one symbol is undesired, e.g. some TSP implementations
 *
 * @author Susanne Merz
 */
public class SwapMutationFunction extends AbstractMutationFunction {

	/**
	 * Sets the mutation probabilities to the designated values.
	 *
	 * @param probabilities
	 *          An array, which contains the mutation probabilities.
	 */
	public SwapMutationFunction(double[] probabilities) {
		this.setMutationProbs(probabilities);
		// as we don't do a random mutation, but rather a swapping we don't need a
		// distribution
		this.setMutationSpectrum(null);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see org.biojavax.ga.functions.MutationFunction#mutate(org.biojava.bio.symbol.SymbolList)
	 */
	public SymbolList mutate(SymbolList seq) throws IllegalAlphabetException,
	    ChangeVetoException, IllegalSymbolException {
		int[] pos = this.generateMutationPositions(seq.length());
		if (pos[0] == pos[1] || pos[0] == 0.0) {
			return seq;
		}
		Edit ed1 = new Edit(pos[0], seq.getAlphabet(), seq.symbolAt(pos[1]));
		Edit ed2 = new Edit(pos[1], seq.getAlphabet(), seq.symbolAt(pos[0]));
		seq.edit(ed1);
		seq.edit(ed2);
		return seq;
	}

	/**
	 * generate an array of two positions on the chromosome that get the highest
	 * values in a random choosing relative to their mutationProbability
	 *
	 * @param seqLength
	 *          the legth of the sequence used
	 * @return an int[] of two positions on the sequence used.
	 */
	private int[] generateMutationPositions(int seqLength) {
		Random n = new Random();
		int[] positions = new int[2];
		// step 1: create an array containing a value for every position
		// on the chromosome. The higher the number the more likely it
		// will be choosen as position for a mutation.
		double[] tempprobs = new double[seqLength];
		int maxIndex = getMutationProbs().length - 1;

		for (int i = 1; i <= seqLength; i++) {
			int index = Math.min(i - 1, maxIndex);
			double mutProb = getMutationProbs()[index] - n.nextDouble();
			if (mutProb >= 0) {
				tempprobs[i - 1] = mutProb;
			} else {
				tempprobs[i - 1] = 0;
			}
		}

		// step 2 find the two highest numbers and return their positions
		double highest = 0;
		double second = 0;
		for (int j = 0; j < tempprobs.length; j++) {
			double current = tempprobs[j];
			if (current > second) {
				// allways put highest number in pos[1], second in pos[0]
				if (current > highest) {
					positions[0] = positions[1];
					second = highest;
					positions[1] = j + 1;
					highest = current;
				} else {
					second = current;
					positions[0] = j + 1;
				}

			}
		}
		// switch places so the lower number is in pos[0]
		if (positions[0] > positions[1]) {
			int temp = positions[0];
			positions[0] = positions[1];
			positions[1] = temp;
		}

		return positions;

	}

}
