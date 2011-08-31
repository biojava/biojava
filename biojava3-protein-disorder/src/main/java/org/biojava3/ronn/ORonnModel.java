/* 
 * @(#)ORonnModel.java	1.0 June 2010
 * 
 * Copyright (c) 2010 Peter Troshin
 *  
 *        BioJava development code
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
package org.biojava3.ronn;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

import org.biojava3.ronn.ModelLoader.Model;
import org.biojava3.ronn.ModelLoader.Threshold;



/**
 * Fully re-factored version of RONN model. Based on the code in C version of
 * RONN.
 * 
 * @author Peter Troshin
 * @version 1.0
 * @since 3.0.2
 */
public final class ORonnModel {

    /**
     * Order probability, corresponds to disorder as 1-order
     */
    private final float disorder_weight;

    private final static int AA_ALPHABET = 19;
    private final static int maxR = 110;
    private final static float coef = 1.0f;
    /**
     * Holds encoded query sequence
     */
    private final short[] seqAA;
    /**
     * Holds query sequence
     */
    private final char[] query;

    private final Model model;

    /**
     * Disorder scores for all residues
     */
    private float[] scores = null;

    final float[] detect() {

	scores = new float[query.length];
	int sResidue;
	int dIndex;
	int r;
	float est, fOrder, pDisor, fDisor;
	final float[][] Z = new float[seqAA.length][ORonnModel.maxR];
	final int[] Q = new int[seqAA.length];
	final Threshold thold = new ModelLoader.Threshold(model.modelNum);

	/*
	 * 19 looks like a size of the sliding window. So for any sequences
	 * shorted than 19 AA the score will be NaN. Original RONN segfault in
	 * such condition
	 */
	for (sResidue = 0; sResidue <= query.length - ORonnModel.AA_ALPHABET; sResidue++) {
	    est = 0.0f;

	    for (dIndex = 0; dIndex < model.numOfDBAAseq; dIndex++) {
		final float[] rho = align(sResidue, dIndex);// search for the
		// maximum alignment between ith peptide from the
		// query and the dIndex-th database sequence
		est += model.W[dIndex] * Math.exp((rho[1] - rho[0]) / rho[0]);
	    }

	    fOrder = (float) (Math.exp(-0.5 * Math.pow(est - thold.mu0, 2.0)
		    / thold.sigma0) / (Math.sqrt(6.28) * thold.sigma0));

	    fDisor = (float) (Math.exp(-0.5 * Math.pow(est - thold.mu1, 2.0)
		    / thold.sigma1) / (Math.sqrt(6.28) * thold.sigma1));

	    pDisor = (float) (disorder_weight * fDisor / ((1.0 - disorder_weight)
		    * fOrder + disorder_weight * fDisor));
	    for (r = sResidue; r < sResidue + ORonnModel.AA_ALPHABET; r++) {
		Z[r][Q[r]] = pDisor;
		Q[r]++;
	    }
	}

	for (sResidue = 0; sResidue < query.length; sResidue++) {
	    est = 0.0f;
	    float[] zRow = Z[sResidue];
	    int numOfIterations = Q[sResidue];
	    for (r = 0; r < numOfIterations; r++) {
		est += zRow[r];
	    }
	    scores[sResidue] = est / numOfIterations;
	}
	return scores;
    }

    public void getScores(final File outfile) throws FileNotFoundException {
	final PrintWriter output = new PrintWriter(outfile);
	if (scores == null) {
	    synchronized (this) {
		if (scores == null) {
		    detect();
		}
	    }
	}
	for (int i = 0; i < scores.length; i++) {
	    output.printf("%c\t%f\n", query[i], scores[i]);
	}
	output.close();
    }

    // sResidue query sequence index and dIndex database sequence index
    private final float[] align(final int sResidue, final int dIndex) {
	int dResidue, r;
	float maxScore = -1000000;
	float rho1 = 0;
	int maxIdx = 0;
	float rho0 = 0;
	short[] dbAARow = model.dbAA[dIndex];
	int numOfIterations = model.Length[dIndex] - ORonnModel.AA_ALPHABET;
	for (dResidue = 0; dResidue <= numOfIterations; dResidue++) {
	    // go though the database sequence for maximised alignment
	    rho1 = 0.0f;
	    for (r = 0; r < ORonnModel.AA_ALPHABET; r++) {
		// go through the query sequence for one alignment
		rho1 += RonnConstraint.Blosum62[seqAA[sResidue + r]][dbAARow[dResidue
			+ r]];
	    }
	    if (rho1 > maxScore) {
		maxScore = rho1;
		maxIdx = dResidue;
	    }
	}
	for (r = 0; r < ORonnModel.AA_ALPHABET; r++) {
	    rho0 += RonnConstraint.Blosum62[dbAARow[maxIdx + r]][dbAARow[maxIdx
		    + r]];
	}
	return new float[] { rho0, maxScore };
    }

    public ORonnModel(final String sequence, final Model model,
	    final float disorder) throws NumberFormatException {
	this.disorder_weight = disorder;
	this.model = model;
	query = sequence.toCharArray();
	seqAA = new short[query.length];
	assert model != null;
	assert model.numOfDBAAseq > 0;
	for (int sResidue = 0; sResidue < sequence.length(); sResidue++) {
	    seqAA[sResidue] = RonnConstraint.INDEX[query[sResidue] - 'A'];
	    if ((seqAA[sResidue] < 0) || (seqAA[sResidue] > 19)) {
		System.err.printf("seqAA[sResidue]=%d(%c)\n", seqAA[sResidue],
			query[sResidue]);
		System.exit(1);
	    }
	}
    }

}
