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

import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.forester.evoinference.distance.NeighborJoining;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * TreeConstructor uses the forester tree library to build phylogenetic trees.
 *
 * @author Scooter Willis
 * @author Aleix Lafita
 * 
 */
public class TreeConstructor<C extends AbstractSequence<D>, D extends Compound>
		extends Thread {

	private static final Logger logger = LoggerFactory
			.getLogger(TreeConstructor.class);

	private TreeType treeType;
	private ScoreMatrixType scoreType;
	private List<NJTreeProgressListener> progessListenerVector = new ArrayList<NJTreeProgressListener>();
	private MultipleSequenceAlignment<C, D> multipleSequenceAlignment = new MultipleSequenceAlignment<C, D>();

	public TreeConstructor(
			MultipleSequenceAlignment<C, D> multipleSequenceAlignment,
			TreeType _treeType,
			ScoreMatrixType _scoreType,
			NJTreeProgressListener _treeProgessListener) {

		treeType = _treeType;
		scoreType = _scoreType;
		addProgessListener(_treeProgessListener);
		this.multipleSequenceAlignment = multipleSequenceAlignment;

	}

	public TreeConstructor(BasicSymmetricalDistanceMatrix _matrix,
			TreeType _treeType,
			TreeConstructionAlgorithm _treeConstructionAlgorithm,
			NJTreeProgressListener _treeProgessListener) {
		matrix = _matrix;
		copyDistanceMatrix = CheckTreeAccuracy.copyMatrix(matrix);
		treeType = _treeType;
		scoreType = _treeConstructionAlgorithm;
		addProgessListener(_treeProgessListener);

	}

	public void outputPhylipDistances(String fileName) throws Exception {
		DistanceMatrix distances = getDistanceMatrix();
		if (distances == null) {
			throw new Exception("distance matrix has not been calculated. "
					+ "Requires process() method to be called first");
		}
		FileOutputStream fo = new FileOutputStream(fileName);
		PrintStream dos = new PrintStream(fo);
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(5);
		df.setMinimumFractionDigits(5);
		for (int row = 0; row < distances.getSize(); row++) {
			dos.print(distances.getIdentifier(row));
			for (int col = 0; col < distances.getSize(); col++) {
				dos.print(" " + df.format(distances.getValue(col, row)));
			}
			dos.println();
		}
		dos.close();
		fo.close();
	}

	public DistanceMatrix getDistanceMatrix() {
		return copyDistanceMatrix;
	}

	public void cancel() {
		// if (njtree != null) {
		// njtree.cancel();
		// }
	}

	boolean verbose = false;
	Phylogeny p = null;
	BasicSymmetricalDistanceMatrix matrix = null;
	DistanceMatrix copyDistanceMatrix = null;

	public void process() throws Exception {

		if (matrix == null) {
			double[][] distances = DistanceCalculator.calculateDistanceMatrix(
					this, multipleSequenceAlignment, scoreType);
			matrix = new BasicSymmetricalDistanceMatrix(
					multipleSequenceAlignment.getSize());
			for (int i = 0; i < matrix.getSize(); i++) {
				matrix.setIdentifier(i, multipleSequenceAlignment
						.getAlignedSequence(i + 1).getAccession().getID());
			}
			for (int col = 0; col < matrix.getSize(); col++) {
				for (int row = 0; row < matrix.getSize(); row++) {
					matrix.setValue(col, row, distances[col][row]);

				}
			}
			copyDistanceMatrix = CheckTreeAccuracy.copyMatrix(matrix);
		}

		final List<Phylogeny> ps = new ArrayList<Phylogeny>();
		final NeighborJoining nj = NeighborJoining.createInstance(verbose);

		ps.add(nj.execute(matrix));
		p = ps.get(0);

	}

	// public void getTreeAccuracy(){
	// CheckTreeAccuracy checkTreeAccuracy = new CheckTreeAccuracy();
	// checkTreeAccuracy.process(p,distanceMatrix );
	// }
	@Override
	public void run() {
		try {
			process();
		} catch (Exception e) {
			logger.error("Exception: ", e);
		}
	}

	public String getNewickString(boolean simpleNewick,
			boolean writeDistanceToParent) throws Exception {
		final PhylogenyWriter w = new PhylogenyWriter();
		StringBuffer newickString = w.toNewHampshire(p, simpleNewick,
				writeDistanceToParent);
		return newickString.toString();
	}

	public void addProgessListener(NJTreeProgressListener treeProgessListener) {
		if (treeProgessListener != null) {
			progessListenerVector.add(treeProgessListener);
		}
	}

	public void removeProgessListener(NJTreeProgressListener treeProgessListener) {
		if (treeProgessListener != null) {
			progessListenerVector.remove(treeProgessListener);
		}
	}

	public void broadcastComplete() {
		for (NJTreeProgressListener treeProgressListener : progessListenerVector) {
			treeProgressListener.complete(this);
		}
	}

	public void updateProgress(String state, int percentage) {
		for (NJTreeProgressListener treeProgressListener : progessListenerVector) {
			treeProgressListener.progress(this, state, percentage);
		}
	}

	public void updateProgress(String state, int currentCount, int totalCount) {
		for (NJTreeProgressListener treeProgressListener : progessListenerVector) {
			treeProgressListener
					.progress(this, state, currentCount, totalCount);
		}
	}

}