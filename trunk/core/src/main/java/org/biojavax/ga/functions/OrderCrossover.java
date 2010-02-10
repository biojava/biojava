/**
 *
 */
package org.biojavax.ga.functions;

import java.util.Iterator;
import java.util.Random;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.functions.AbstractCrossOverFunction;
import org.biojavax.ga.functions.GACrossResult;
import org.biojavax.ga.functions.SimpleGACrossResult;

/**
 * This does a 2-point-crossover on two chromosomes keeping the Symbols in each
 * chromosome constant. The method is commonly named OX - operator
 *
 * @author Susanne Merz
 */
public class OrderCrossover extends AbstractCrossOverFunction {

	/**
	 * Sets the maximal number of crossover points to two and the crossover
	 * probability to 0.5 and initializes this object.
	 */
	public OrderCrossover() {
		this.setMaxCrossOvers(2);
		double[] pro = {0.5};
		this.setCrossOverProbs(pro);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see org.biojavax.ga.functions.CrossOverFunction#performCrossOver(org.biojava.bio.symbol.SymbolList,
	 *      org.biojava.bio.symbol.SymbolList)
	 */
	public GACrossResult performCrossOver(SymbolList chromA, SymbolList chromB)
	    throws ChangeVetoException {
		int[] crosslocations = this.generateCrossoverLocations(chromA.length());
		// only do a chrossover if there are valid positions.
		if (crosslocations[1] != 0.0) {
			// if we found only one crossover position, we assume the
			// second position to be the end of teh chromosome
			if (crosslocations[0] == 0.0) {
				crosslocations[0] = crosslocations[1];
				crosslocations[1] = chromA.length();
			}
			Alphabet alphA = chromA.getAlphabet();
			Alphabet alphB = chromB.getAlphabet();
			SymbolList snipA = chromA.subList(crosslocations[0], crosslocations[1]);
			SymbolList snipB = chromB.subList(crosslocations[0], crosslocations[1]);
			SymbolList snipAwoB = snipA.subList(1, snipA.length());
			SymbolList snipBwoA = snipB.subList(1, snipB.length());
			// remove all elements of snipA from snipBwoA and vice versa
			for (int i = 1; i <= snipB.length(); i++) {
				int j = 1;
				Iterator it = snipAwoB.iterator();
				boolean notfound = true;
				while (it.hasNext() && notfound) {

					if (it.next().equals(snipB.symbolAt(i))) {
						try {
							snipAwoB.edit(new Edit(j, 1, SymbolList.EMPTY_LIST));
						} catch (Exception e) {
							e.printStackTrace();
						}
						break;
					}
					j++;
				}
				Iterator it2 = snipBwoA.iterator();
				j = 1;
				while (it2.hasNext()) {
					if (it2.next().equals(snipA.symbolAt(i))) {
						try {
							snipBwoA.edit(new Edit(j, 1, SymbolList.EMPTY_LIST));
							break;
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
					j++;
				}
			}
			// remove all elements of snipA and snipB from both chromsomes
			int sniplength = (crosslocations[1] - crosslocations[0]) + 1;
			try {
				// remove snipA from chromB and snipB from chromB
				chromA.edit(new Edit(crosslocations[0], sniplength,
				    SymbolList.EMPTY_LIST));
				chromB.edit(new Edit(crosslocations[0], sniplength,
				    SymbolList.EMPTY_LIST));
				// remove elements of snipA from chromB
				Iterator it1 = snipA.iterator();
				while (it1.hasNext()) {
					Object current = it1.next();
					Iterator it2 = chromB.iterator();
					Iterator it3 = chromA.iterator();
					int position = 1;
					while (it2.hasNext()) {
						if (it2.next().equals(current)) {
							chromB.edit(new Edit(position, 1, SymbolList.EMPTY_LIST));
							break;
						}
						position++;
					}
				}
				// remove elements of snipB from chromA
				it1 = snipB.iterator();
				while (it1.hasNext()) {
					Object current = it1.next();
					Iterator it2 = chromA.iterator();
					Iterator it3 = chromB.iterator();
					int position = 1;
					boolean notfound = true;
					while (it2.hasNext() && notfound) {
						if (it2.next().equals(current)) {
							chromA.edit(new Edit(position, 1, SymbolList.EMPTY_LIST));
							notfound = false;
							break;
						}
						position++;
					}
				}
			} catch (IllegalAlphabetException e) {
				System.out.println("Sorry, you used an illegal alphabet");
				e.printStackTrace();
			}
			// put parts together in the order snipXwoY-snipY-X
			try {
				chromA.edit(new Edit(1, 0, snipB));
				chromA.edit(new Edit(1, 0, snipAwoB));
				chromB.edit(new Edit(1, 0, snipA));
				chromB.edit(new Edit(1, 0, snipBwoA));
			} catch (IllegalAlphabetException e) {
				e.printStackTrace();
			}
			// if the chromosomes are not cycle-invariant we need to
			// adjust positions
			int newposition = 1;
			Iterator it = chromA.iterator();
			while (it.hasNext()) {
				if (it.next().equals(snipB.symbolAt(1))) {
					// newposition++;
					break;
				}
				newposition++;
			}
			int versatz = crosslocations[0] - newposition;

			if (versatz < 0) {
				versatz = Math.abs(versatz);
				// remove from front and add to end
				try {
					SymbolList temp = chromA.subList(1, versatz);
					chromA.edit(new Edit(1, versatz, SymbolList.EMPTY_LIST));
					chromA.edit(new Edit(chromA.length() + 1, 0, temp));
					temp = chromB.subList(1, versatz);
					chromB.edit(new Edit(1, versatz, SymbolList.EMPTY_LIST));
					chromB.edit(new Edit(chromB.length() + 1, 0, temp));

				} catch (IllegalAlphabetException e) {
					e.printStackTrace();
				}
			}

			else if (versatz > 0) {
				// remove from end and add to front
				try {
					versatz--;
					int start = chromA.length() - versatz;
					SymbolList temp = chromA.subList(start, chromA.length());
					chromA.edit(new Edit(start, chromA.length() - (start - 1),
					    SymbolList.EMPTY_LIST));
					chromA.edit(new Edit(1, 0, temp));
					temp = chromB.subList(start, chromB.length());
					chromB.edit(new Edit(start, chromB.length() - (start - 1),
					    SymbolList.EMPTY_LIST));
					chromB.edit(new Edit(1, 0, temp));

				} catch (IllegalAlphabetException e) {
					e.printStackTrace();
				}
			}

			PointLocation[] points = {new PointLocation(crosslocations[0]),
			    new PointLocation(crosslocations[1])};
			GACrossResult result = new SimpleGACrossResult(points, new SymbolList[] {
			    chromA, chromB});

			return result;
		} else return null;
	}

	/**
	 * generate a pair of locations where the crossover should be made do this in
	 * two steps for retracabilities sake. the location array will always list the
	 * earlier position first.
	 *
	 * @param chromlength
	 *          the lenghth of the chromosomes, needed to decide if the
	 *          probabilityarray needs to be extended.
	 * @return an int[] containing two positions on the chomosome
	 */
	private synchronized int[] generateCrossoverLocations(int chromlength) {
		Random n = new Random();
		int[] locations = new int[2];
		// step 1 create an array containing a value for every position
		// on the chromosome. The higher the number the more likely it
		// will be choosen as crossoverposition.
		double[] crosses = new double[chromlength];
		int maxIndex = getCrossOverProbs().length - 1;

		for (int i = 1; i <= chromlength; i++) {
			int index = Math.min(i - 1, maxIndex);
			double crossProb = getCrossOverProbs()[index] - n.nextDouble();
			if (crossProb >= 0) {
				crosses[i - 1] = crossProb;
			} else {
				crosses[i - 1] = 0;
			}
		}

		// step 2 find the two highest numbers and return their positions
		double highest = 0;
		double second = 0;
		for (int j = 0; j < crosses.length; j++) {
			double current = crosses[j];
			if (current > second) {
				// allways put highest number in pos[1], second in pos[0]
				if (current > highest) {
					locations[0] = locations[1];
					second = highest;
					locations[1] = j + 1;
					highest = current;
				} else {
					second = current;
					locations[0] = j + 1;
				}

			}
		}
		// switch places so the lower number is in pos[0]
		if (locations[0] > locations[1]) {
			int temp = locations[0];
			locations[0] = locations[1];
			locations[1] = temp;
		}

		return locations;
	}

}
