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
package org.biojava.nbio.survival.cox;

/**
 * Not used and probably should be deleted
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvivalInfoIndex implements Comparable<SurvivalInfoIndex> {

	private double time;
	private int event;
	private int index;
	private double[] data;

	/**
	 *
	 * @param t
	 * @param e
	 * @param i
	 */
	public SurvivalInfoIndex(double t, int e, int i) {
		time = t;
		event = e;
		index = i;
	}

	/**
	 *
	 * @param t
	 * @param e
	 * @param i
	 * @param d
	 */
	public SurvivalInfoIndex(double t, int e, int i, double[] d) {
		time = t;
		event = e;
		index = i;
		data = d;
	}

	/**
	 *
	 * @param t
	 * @param e
	 * @param i
	 * @param d
	 */
	public SurvivalInfoIndex(double t, int e, int i, double d) {
		time = t;
		event = e;
		index = i;
		data = new double[1];
		data[0] = d;
	}

	@Override
	public String toString() {
		return "t=" + time + " e=" + event + " o=" + index;
	}
	//    double CompNum4Sort(double[] a, double[] b) {
	//(time - time - (event -event) /1024)
	//    return (a[0] - b[0] - (a[1] - b[1]) / 1024);
	// }

	@Override
	public int compareTo(SurvivalInfoIndex o) {
	//    double compare = (this.time - o.time - (this.event - o.event) / 1024);
		if (time < o.time) {
			return -1;
		} else if (time > o.time) {
			return 1;
		} else {
			if (this.event == o.event) {
				return 0;
			} else if (event == 1) {
				return -1;
			} else {
				return 1;
			}
		}

	}
}
