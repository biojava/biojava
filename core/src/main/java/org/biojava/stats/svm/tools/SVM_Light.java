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
 * @(#)SVM_Light.java      0.1 00/01/15
 *
 * By Thomas Down <td2@sanger.ac.uk>
 */

package org.biojava.stats.svm.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.StringTokenizer;

import org.biojava.stats.svm.CachingKernel;
import org.biojava.stats.svm.PolynomialKernel;
import org.biojava.stats.svm.RadialBaseKernel;
import org.biojava.stats.svm.SVMClassifierModel;
import org.biojava.stats.svm.SVMKernel;
import org.biojava.stats.svm.SimpleSVMClassifierModel;
import org.biojava.stats.svm.SparseVector;

/**
 * @author Thomas Down
 * @author Greg Cox
 */

public class SVM_Light {
    public static class LabelledVector {
	private SparseVector v;
	private double label;
	private String comment = null;

	public LabelledVector(SparseVector v, double label) {
	    this.v = v;
	    this.label = label;
	}

	public LabelledVector(SparseVector v, double label, String comment) {
	    this.v = v;
	    this.label = label;
	    this.comment = comment;
	}

	public SparseVector getVector() {
	    return v;
	}

	public double getLabel() {
	    return label;
	}

	public String getComment() {
	    return comment;
	}
    }

    public static LabelledVector parseExample(String ex)
	throws NumberFormatException
    {
	String comment = null;
	int hashPos = ex.indexOf('#');
	if (hashPos >= 0) {
	    comment = ex.substring(hashPos + 1);
	    ex = ex.substring(0, hashPos);
	}

	StringTokenizer toke = new StringTokenizer(ex);
	double label = Double.parseDouble(toke.nextToken());

  	int size = toke.countTokens();
	SparseVector v = new SparseVector(size);
	while (toke.hasMoreTokens()) {
	    String dim = toke.nextToken();
	    int cut = dim.indexOf(':');
	    if (cut < 0) {
		throw new NumberFormatException("Bad dimension "+dim);
	    }
	    int dnum = Integer.parseInt(dim.substring(0, cut));
	    double dval = Double.parseDouble(dim.substring(cut + 1));
	    v.put(dnum, dval);
	}

	return new LabelledVector(v, label, comment);
    }

    public static String vectorToString(SparseVector v) {
	StringBuffer sb = new StringBuffer();
	boolean first = true;

	for (int i = 0; i < v.size(); ++i) {
	    double x = v.getValueAtIndex(i);

	    if (first) {
		first = false;
	    } else {
		sb.append(' ');
	    }

	    sb.append(v.getDimAtIndex(i));
	    sb.append(':');
	    sb.append(x);
	}
	return sb.substring(0);
    }

    public static SVMClassifierModel readModelFile(String fileName)
	throws IOException
    {
	BufferedReader r = new BufferedReader(new FileReader(fileName));
	r.readLine(); // format
	String kType = firstToken(r.readLine());
	String dParam = firstToken(r.readLine());
	String gParam = firstToken(r.readLine());
	r.readLine(); // s
	r.readLine(); // r
	r.readLine(); // u
	r.readLine(); // numSV
	String threshString = firstToken(r.readLine());

	SVMKernel kernel = null;
	try {
	    switch (Integer.parseInt(kType)) {
	    case 0:
		kernel = SparseVector.kernel;
		break;
	    case 1:
		int order = Integer.parseInt(dParam);
		PolynomialKernel k = new PolynomialKernel();
		k.setOrder(order);
		k.setNestedKernel(SparseVector.kernel);
		kernel = k;
		break;
  	    case 2:
		RadialBaseKernel rbk = new RadialBaseKernel();
		double width = Double.parseDouble(gParam);
		rbk.setWidth(width);
		rbk.setNestedKernel(SparseVector.kernel);
		kernel = rbk;
		break;
  	    default:
		throw new IOException("Couldn't create kernel");
	    }

	    SimpleSVMClassifierModel model = new SimpleSVMClassifierModel(kernel);
	    model.setThreshold(Double.parseDouble(threshString));
	    String line;
	    while ((line = r.readLine()) != null) {
		LabelledVector ex = parseExample(line);
		model.addItemAlpha(ex.getVector(), ex.getLabel());
	    }
	    r.close();

	    return model;
	} catch (NumberFormatException ex) {
	    throw new IOException("Couldn't parse model file");
	}
    }

  public static void writeModelFile(SVMClassifierModel model, String fileName)
  throws IOException {
    SVMKernel k = model.getKernel();

    int kType = 0;
    int d = 3;
    double g = 1;
    double s = 1;
    double r = 1;
    String u = "empty";

    while(k instanceof CachingKernel) {
      k = ((CachingKernel) k).getNestedKernel();
    }
    
    if (k == SparseVector.kernel) {
      kType = 0;
    } else if (k instanceof PolynomialKernel) {
      kType = 1;
      d = (int) ((PolynomialKernel) k).getOrder();
    } else if (k instanceof RadialBaseKernel) {
      kType = 2;
      g = ((RadialBaseKernel) k).getWidth();
    } else {
      throw new IOException("Can't write SVM_Light file with kernel type " + k.getClass().toString());
    }


    PrintWriter pw = new PrintWriter(new FileWriter(fileName));
    pw.println("SVM-light Version V3.01");
    pw.println("" + kType + " # kernel type");
    pw.println("" + d + " # kernel parameter -d");
    pw.println("" + g + " # kernel parameter -g");
    pw.println("" + s + " # kernel parameter -s");
    pw.println("" + r + " # kernel parameter -r");
    pw.println(u + " # kernel parameter -u");

    int numSV = 0;
    for(Iterator i = model.items().iterator(); i.hasNext(); ) {
      Object item = i.next();
      if (model.getAlpha(item) != 0) {
        numSV++;
      }
    }

    pw.println("" + numSV + " # number of support vectors");
    pw.println("" + model.getThreshold() + " # threshold b");

    for(Iterator i = model.items().iterator(); i.hasNext(); ) {
      Object item = i.next();
      if (model.getAlpha(item) == 0) {
        continue;
      }
      pw.print(model.getAlpha(i));

      SparseVector v = (SparseVector) item;
      for (int j = 0; j <= v.maxIndex(); ++j) {
        double x = v.get(j);
        if (x != 0.0)
        pw.print(" " + j + ":" + x);
      }
      pw.println("");
    }

    pw.close();
  }

  public static String firstToken(String s) {
    int ndx = s.indexOf(" ");
    if (ndx < 0) {
	    return s;
    } else {
      return s.substring(0, ndx);
    }
  }
}
