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

package org.biojava.nbio.structure.symmetry.core;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Peter
 */
public class HelixLayers {
	private List<Helix> helices = new ArrayList<Helix>();
	private double symmetryDeviation = 0;

	public int size() {
		return helices.size();
	}

	public void addHelix(Helix helix) {
		helices.add(helix);
	}

	public Helix getHelix(int index) {
		return helices.get(index);
	}

	/*
	 * Returns Helix with lowest twist angle
	 */
	public Helix getByLowestAngle() {
		double angle = Double.MAX_VALUE;
		Helix lowest = null;
		for (Helix helix: helices) {
			if (helix.getAngle() < angle) {
				angle = helix.getAngle();
				lowest = helix;
			}
		}
		return lowest;
	}

	/*
	 * Returns Helix with largest number of intermolecular contacts
	 * between repeat units
	 */
	public Helix getByLargestContacts() {
		double contacts = 0;
		Helix largest = null;
		for (Helix helix: helices) {
			if (helix.getContacts() > contacts) {
				contacts = helix.getContacts();
				largest = helix;
			}
		}
		return largest;
	}

	/*
	 * Returns Helix that has the largest number of contacts, besides
	 * the Helix with the lowest twist angle
	 */
	public Helix getByLargestContactsNotLowestAngle() {
		double contacts = 0;
		Helix lowest = getByLowestAngle();
		// TODO why are there helices with almost identical helix parameters??
		double angle = lowest.getAngle() + 0.05;
		Helix largest = null;
		for (Helix helix: helices) {
			if (helix == lowest) {
				continue;
			}
			if (helix.getContacts() > contacts && helix.getAngle() > angle) {
				contacts = helix.getContacts();
				largest = helix;
			}
		}
		if (largest == null) {
			return lowest;
		}
		return largest;
	}

	/**
	 * Returns QuatSymmetryScores averaged over all rotations
	 * (except the first rotation, which is the unit operation E)
	 * @return mean scores average over rotations
	 */
	public QuatSymmetryScores getScores() {
		QuatSymmetryScores scores = new QuatSymmetryScores();

		double[] values = new double[helices.size()];

		// minRmsd
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getMinRmsd();
		}
		scores.setMinRmsd(minScores(values));

		// maxRmsd
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getMaxRmsd();
		}
		scores.setMaxRmsd(maxScores(values));

		// Rmsd
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getRmsd();
		}
		scores.setRmsd(averageScores(values));

		// minTm
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getMinTm();
		}
		scores.setMinTm(minScores(values));

		// maxTm
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getMaxTm();
		}
		scores.setMaxTm(maxScores(values));

		// Tm
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getTm();
		}
		scores.setTm(averageScores(values));

		// Rmsd subunit centers
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getRmsdCenters();
		}
		scores.setRmsdCenters(averageScores(values));

		// TmIntra
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getTmIntra();
		}
		scores.setTmIntra(averageScores(values));

		// RmsdIntra
		for (int i = 0; i < helices.size(); i++) {
			values[i] = helices.get(i).getScores().getRmsdIntra();
		}
		scores.setRmsdIntra(averageScores(values));

		// SymDeviation
		scores.setSymDeviation(symmetryDeviation);

		return scores;
	}

	/**
	 * @param symmetryDeviation the symmetryDeviation to set
	 */
	public void setSymmetryDeviation(double symmetryDeviation) {
		this.symmetryDeviation = symmetryDeviation;
	}

	private double averageScores(double[] scores) {
		double sum = 0;
		for (double s: scores) {
			sum += s;
		}
		return sum/scores.length;
	}

	private double minScores(double[] scores) {
		double score = Double.MAX_VALUE;
		for (double s: scores) {
			score = Math.min(score, s);
		}
		return score;
	}

	private double maxScores(double[] scores) {
		double score = Double.MIN_VALUE;
		for (double s: scores) {
			score = Math.max(score, s);
		}
		return score;
	}

	public void clear() {
		helices.clear();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Helices: ").append(size()).append("\n");
		for (Helix s: helices) {
			sb.append(s.toString()).append("\n");
		}
		return sb.toString();
	}



}
