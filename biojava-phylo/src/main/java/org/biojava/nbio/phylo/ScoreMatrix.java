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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.phylo;

/**
 * This source file is derived from jalview.schemes.ScoreMatrix to provide tree
 * construction in the biojava package to minimize the number of external
 * JalView classes required.
 *
 * Jalview - A Sequence Alignment Editor and Viewer (Version 2.4) Copyright (C)
 * 2008 AM Waterhouse, J Procter, G Barton, M Clamp, S Searle
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 */
public class ScoreMatrix {

	/** Reference to the substitution matrix. */
	private int[][] matrix;

	private ScoreMatrixType type;

	public ScoreMatrix(int[][] matrix, ScoreMatrixType type) {
		this.matrix = matrix;
		this.type = type;
	}

	public boolean isProtein() {
		return ScoreMatrixType.isProtein(type);
	}

	public boolean isNucleotide() {
		return ScoreMatrixType.isNucleotide(type);
	}

	public int[][] getMatrix() {
		return matrix;
	}

	/**
	 * Return the pairwise score between two aligned characters.
	 * 
	 * @param A1 character if sequence 1
	 * @param A2 character of sequence 2
	 * @return score for substituting first char in A1 with first char in A2
	 */
	public int getPairwiseScore(String A1, String A2) {
		return getPairwiseScore(A1.charAt(0), A2.charAt(0));
	}

	public int getPairwiseScore(char c, char d) {

		int a = (isNucleotide()) ? ResidueProperties.aaIndex[c]
				: ResidueProperties.nucleotideIndex[c];
		int b = (isProtein()) ? ResidueProperties.aaIndex[d]
				: ResidueProperties.nucleotideIndex[d];

		return matrix[a][b];
	}

}
