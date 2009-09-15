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
 * @(#)TrainRegression.java      0.1 00/01/15
 *
 * By Matthew Pocock <mrp@sanger.ac.uk>
 */

package org.biojava.stats.svm.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import org.biojava.stats.svm.PolynomialKernel;
import org.biojava.stats.svm.SMORegressionTrainer;
import org.biojava.stats.svm.SVMRegressionModel;
import org.biojava.stats.svm.SparseVector;

/**
 * @author Ewan Birney
 * @author Matthew Pocock
 */
public class TrainRegression {
  public static void main(String[] args) throws Throwable {
    if (args.length < 2) {
	    throw new Exception("usage: stats.svm.tools.TrainRegression <train_examples> <model_file>");
    }
    String trainFile = args[0];

    List examples = new ArrayList();
    BufferedReader r = new BufferedReader(new FileReader(trainFile));
    String line;

    while ((line = r.readLine()) != null) {
	    if (line.length() == 0 || line.startsWith("#")) {
        continue;
      }
	    examples.add(SVM_Light.parseExample(line));
    }
    r.close();
  
  	SVMRegressionModel model = new SVMRegressionModel(examples.size());
    double[] target = new double[examples.size()];
    for (int i = 0; i < examples.size(); ++i) {
	    SVM_Light.LabelledVector ex = (SVM_Light.LabelledVector) examples.get(i);
	    model.addVector(ex.getVector());
	    target[i] = ex.getLabel();
    }

    PolynomialKernel k = new PolynomialKernel();
    k.setNestedKernel(SparseVector.kernel);
    k.setOrder(2);
    model.setKernel(k);
    System.out.println("Calculating kernel " + k);
    model.calcKernel();
    SMORegressionTrainer trainer = new SMORegressionTrainer();
    trainer.setEpsilon(0.00000000001);
    trainer.setC(1000);
    System.out.println("\nTraining");
    trainer.trainModel(model, target);
    System.out.println("\nDone");

    for (int i=0; i < model.size(); ++i) {
      System.err.println("y=" + target[i] + "\tf(x)=" + model.internalClassify(i)
                         + "    (" + model.getAlpha(i) + ",\t"
                         + model.getAlphaStar(i) + ")" + "\t" + (model.internalClassify(i) - model.getThreshold()));
    }
    System.err.println("b=" + model.getThreshold());

//    SVM_Light.writeModelFile(model, modelFile);
  }
}
