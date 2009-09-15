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
 * @(#)Classify.java      0.1 00/01/15
 *
 * By Thomas Down <td2@sanger.ac.uk>
 */

package org.biojava.stats.svm.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.biojava.stats.svm.SVMClassifierModel;

/**
 * @author Ewan Birney;
 * @author Matthew Pocock
 */
public class Classify {
    public static void main(String[] args) throws Throwable {
	if (args.length < 3) {
	    throw new Exception("usage: stats.svm.tools.Classify <model> <test_examples> <results_log>");
	}

	String modelName = args[0];
	String examplesName = args[1];
	String resultsName = args[2];

	SVMClassifierModel model = SVM_Light.readModelFile(modelName);

	BufferedReader r = new BufferedReader(new FileReader(examplesName));
	PrintWriter w = new PrintWriter(new FileWriter(resultsName));
	
	String line;
	int right = 0, wrong = 0;

	while ((line = r.readLine()) != null) {
	    if (line.length() == 0 || line.startsWith("#"))
		continue;
	    SVM_Light.LabelledVector ex = SVM_Light.parseExample(line);
	    double result = model.classify(ex.getVector());
	    w.println(result);

	    if (sign(result) == sign(ex.getLabel()))
		right++;
	    else
		wrong++;
	}

	System.out.println("" + ((double) right / (right + wrong))*100 + "% correct");

	r.close();
	w.close();
    }

    public static int sign(double d) {
	if (d < 0)
	    return -1;
	else if (d == 0)
	    return 0;
	return 1;
    }
}
