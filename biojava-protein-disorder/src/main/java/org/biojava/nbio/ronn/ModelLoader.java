/* 
 *  @(#)ModelLoader.java	1.0 June 2010
 * 
 *  Copyright (c) 2010 Peter Troshin
 *  
 *  JRONN version: 3.1     
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
 
package org.biojava.nbio.ronn;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Class that loads data from the model files into {@link Model} objects
 * 
 * @author Peter Troshin
 * @version 1.0
 * @since 3.0.2

 */
public final class ModelLoader {

	private static final Logger logger = LoggerFactory.getLogger(ModelLoader.class);

    /**
     * Represents a Threshold
     * 
     */
    public final static class Threshold {

	final float mu0;
	final float mu1;
	final float sigma0;
	final float sigma1;

	public Threshold(final int modelNum) {
	    final float[] values = RonnConstraint.Threshold
		    .getTreshold(modelNum);
	    mu0 = values[0];
	    mu1 = values[1];
	    sigma0 = values[2];
	    sigma1 = values[3];
	}

    }

    /**
     * Represent a RONN model
     * 
     */
    public static class Model {

	/**
	 * Stores encoded sequences from the model similar to seqAA
	 * 
	 * 190 is a maximum length of the sequence in the model
	 */
	final short[][] dbAA;// = new short[RonnConstraint.maxD][190];
	/**
	 * This array contains the length of each sequence from the model file
	 * The length of this array vary with the number of sequences in the
	 * mode
	 */
	final short[] Length;// = new int[RonnConstraint.maxD];
	/**
	 * Holds the values from the second column of model file
	 */
	final float[] W;// = new float[RonnConstraint.maxD];

	int numOfDBAAseq;
	int modelNum = -1;

	public Model(final int modelNum, final int numberofSequence) {
	    this.modelNum = modelNum;
	    numOfDBAAseq = numberofSequence;
	    dbAA = new short[numberofSequence][190];
	    Length = new short[numberofSequence];
	    W = new float[numberofSequence];
	}

	@Override
	public int hashCode() {
	    final int prime = 31;
	    int result = 1;
	    result = prime * result + Arrays.hashCode(Length);
	    result = prime * result + Arrays.hashCode(W);
	    result = prime * result + Arrays.hashCode(dbAA);
	    result = prime * result + modelNum;
	    result = prime * result + numOfDBAAseq;
	    return result;
	}

	@Override
	public boolean equals(final Object obj) {
	    if (this == obj) {
		return true;
	    }
	    if (obj == null) {
		return false;
	    }
	    if (getClass() != obj.getClass()) {
		return false;
	    }
	    final Model other = (Model) obj;
	    if (!Arrays.equals(Length, other.Length)) {
		return false;
	    }
	    if (!Arrays.equals(W, other.W)) {
		return false;
	    }
	    if (!Arrays.equals(dbAA, other.dbAA)) {
		return false;
	    }
	    if (modelNum != other.modelNum) {
		return false;
	    }
	    if (numOfDBAAseq != other.numOfDBAAseq) {
		return false;
	    }
	    return true;
	}

	@Override
	public String toString() {
	    return "Model [modelNum=" + modelNum + ", numOfDBAAseq="
		    + numOfDBAAseq + "]";
	}

    }

    private static final Map<Integer, Model> models = new HashMap<Integer, Model>();

    public Model getModel(final int modelNum) {
	return ModelLoader.models.get(modelNum);
    }

    void loadModels() throws NumberFormatException, IOException {

	for (int i = 0; i < 10; i++) {
	    final BufferedReader bfr = new BufferedReader(
		    new InputStreamReader(ModelLoader.class.getResourceAsStream(
			    "model" + i + ".rec"),
			    "ISO-8859-1"));
	    String line = null;
	    final Scanner scan = new Scanner(bfr);
	    scan.useDelimiter(System.getProperty("line.separator"));
	    final int numberOfSeqs = scan.nextInt();
	    final Model model = new Model(i, numberOfSeqs);
	    // ignore this one, its always 19 defined in RonnConstrain
	    scan.nextInt();
	    for (int j = 0; j < numberOfSeqs; j++) {
		line = scan.next();
		final char[] dbseq = line.trim().toCharArray();
		assert dbseq.length < Short.MAX_VALUE;
		model.Length[j] = (short) dbseq.length;
		for (int dResidue = 0; dResidue < dbseq.length; dResidue++) {
		    model.dbAA[j][dResidue] = RonnConstraint.INDEX[dbseq[dResidue] - 'A'];
		    assert !((model.dbAA[j][dResidue] < 0) || (model.dbAA[j][dResidue] > 19));
		}
		model.W[j] = scan.nextFloat();
	    }
	    ModelLoader.models.put(model.modelNum, model);
	    bfr.close();
	    scan.close();
	}
    }

    public static void main(final String[] args) throws NumberFormatException,
	    IOException {
	final ModelLoader loader = new ModelLoader();
	loader.loadModels();
	logger.info("{}", ModelLoader.models.get(0));
	logger.info("{}", ModelLoader.models.get(9));
	logger.info("{}", ModelLoader.models.size());
    }
}
