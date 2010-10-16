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
 * @(#)Train.java      0.1 00/01/15
 *
 * By Thomas Down <td2@sanger.ac.uk>
 */

package org.biojava.stats.svm.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Iterator;

import org.biojava.stats.svm.CachingKernel;
import org.biojava.stats.svm.DiagonalCachingKernel;
import org.biojava.stats.svm.NormalizingKernel;
import org.biojava.stats.svm.PolynomialKernel;
import org.biojava.stats.svm.SMOTrainer;
import org.biojava.stats.svm.SVMClassifierModel;
import org.biojava.stats.svm.SVMTarget;
import org.biojava.stats.svm.SimpleSVMTarget;
import org.biojava.stats.svm.SparseVector;
import org.biojava.stats.svm.TrainingEvent;
import org.biojava.stats.svm.TrainingListener;

/**
 * @author Ewan Birney
 * @author Matthew Pocock
 */
public class Train {
  public static void main(String[] args) throws Throwable {
    if (args.length != 2) {
	    throw new Exception("usage: stats.svm.tools.Classify <train_examples> <model_file>");
    }
    String trainFile = args[0];
    String modelFile = args[1];

    BufferedReader r = new BufferedReader(new FileReader(trainFile));
    String line;

    SVMTarget target = new SimpleSVMTarget();
    while ((line = r.readLine()) != null) {
      if (line.length() == 0 || line.startsWith("#")) {
        continue;
      }
	    SVM_Light.LabelledVector ex = SVM_Light.parseExample(line);
      target.addItemTarget(ex.getVector(), ex.getLabel());
    }
    r.close();

    PolynomialKernel pK = new PolynomialKernel();
    pK.setOrder(2.0);
    pK.setNestedKernel(SparseVector.kernel);
    DiagonalCachingKernel gcK = new DiagonalCachingKernel();
    gcK.setNestedKernel(pK);
    NormalizingKernel nK = new NormalizingKernel();
    nK.setNestedKernel(gcK);
    CachingKernel cK = new CachingKernel();
    cK.setNestedKernel(nK);

    SMOTrainer trainer = new SMOTrainer();
    trainer.setEpsilon(1.0e-9);
    trainer.setC(1000);
    TrainingListener tl = new TrainingListener() {
	    public void trainingCycleComplete(TrainingEvent e) {
        System.out.print('.');
	    }
	    public void trainingComplete(TrainingEvent e) {
        System.out.println("");
	    }
    };
    
    System.out.println("Training");
	  SVMClassifierModel model = trainer.trainModel(target, cK, tl);
    System.out.println("Done");

    for(Iterator i = target.items().iterator(); i.hasNext(); ) {
      Object item = i.next();
	    System.out.println(target.getTarget(item) + "\t" +
                         model.classify(item) + "\t(" +
                         model.getAlpha(item) + ")");
    }

    SVM_Light.writeModelFile(model, modelFile);
  }
}
