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

package org.biojava.stats.svm;

/**
 * Train a regression support vector machine using the Sequential Minimal
 * Optimization algorithm.  See "A Tutorial on Support Vector Regression"
 * by Smola and Scholkopf.
 *
 * <p>
 * WARNING: This doesn't work right now -- use and fix (both
 * at own risk...)
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */

public class SMORegressionTrainer {
  private double C = 1000;
  private double epsilon = 0.000001;
    
  // Working variables for the trainer: protected by the
  // synchronization on trainModel.

  private SVMRegressionModel model;
  private double[] target;
  private double[] E;

  public void setC(double C) {
    this.C = C;
  }

  public void setEpsilon(double epsilon) {
    this.epsilon = epsilon;
  }

  private boolean takeStep(int i1, int i2) {
    //System.out.print("+");

    if (i1 == i2) {
	    return false;
    }
    
    double y1 = target[i1];
    double y2 = target[i2];
    double alpha1 = model.getAlpha(i1);
    double alpha2 = model.getAlpha(i2);
    double alpha1star = model.getAlphaStar(i1);
    double alpha2star = model.getAlphaStar(i2);
    double phi1 = getError(i1);
    double phi2 = getError(i2);

    System.out.println("y1=" + y1 + "\ty2=" + y2);
    System.out.println("alpha1=" + alpha1 + "\talpha2=" + alpha2);
    System.out.println("alpha1star=" + alpha1star + "\talpha2star=" + alpha2star);
    System.out.println("phi1=" + phi1 + "\tphi2=" + phi2);
    
    double k11 = model.getKernelValue(i1, i1);
    double k12 = model.getKernelValue(i1, i2);
    double k22 = model.getKernelValue(i2, i2);
    double eta = k11 + k22 - 2.0 * k12; // from improvement
    // double eta = 2.0 * k12 - k11 - k22; // from tutorial, but always gives negative eta

    System.out.println("k11=" + k11 + "\tk12=" + k12 + "\tk22=" + k22);
    System.out.println("eta=" + eta);
    
    boolean case1, case2, case3, case4, finnished, changed;
    double deltaphi = phi2 - phi1;
    case1 = case2 = case3 = case4 = finnished = changed = false;

    System.out.println("deltaphi=" + deltaphi);
    
    double L, H;
    double a1 = 0.0, a2 = 0.0;
    double gamma = alpha1 - alpha1star + alpha2 - alpha2star;
    System.out.println("gamma=" + gamma);
    System.out.println("epsilon=" + epsilon);
    
    if(eta <= 0) {
      System.out.println("Negative eta");
      return false;
    }
    while (!finnished) {
      if(!case1 &&
      (alpha1 > 0.0 || (alpha1star == 0.0 && deltaphi > 0.0)) &&
      (alpha2 > 0.0 || (alpha2star == 0.0 && deltaphi < 0.0))
      ) {
        L = Math.max(0.0, gamma - C);
        H = Math.min(C, gamma);
        System.out.println("L=" + L + "\tH=" + H);
        if (L < H) {
          a2 = alpha2 - deltaphi / eta;
          System.out.println("Ideal a2 = " + a2);
          a2 = Math.min(a2, H);
          a2 = Math.max(L, a2);
          a1 = alpha1 - (a2 - alpha2);
          System.out.println("a1=" + a1 + ", a2=" + a2);
          // updatae alpah1, alpha2 if change greater than some eps
          if(Math.abs(a2 - alpha2) > epsilon * (a2 + alpha2 + 1.0 + epsilon)) {
            model.setAlpha(i1, a1);
            model.setAlpha(i2, a2);
            alpha1 = a1;
            alpha2 = a2;
            System.out.println("case1 worked");
            changed = true;
          } else {
            System.out.println("case1: change too small: " + (a2 - alpha2));
          }
        } else {
          System.out.println("case1: L > H");
          finnished = true;
        }
        case1 = true;
      } else if (!case2 &&
      (alpha1     > 0.0 || (alpha1star == 0.0 && deltaphi > 2.0 * epsilon)) &&
      (alpha2star > 0.0 || (alpha2     == 0.0 && deltaphi > 2.0 * epsilon))
      ) {
        L = Math.max(0.0, gamma);
        H = Math.min(C, C + gamma);
        System.out.println("L=" + L + "\tH=" + H);
        if(L < H) {
          a2 = alpha2star + (deltaphi - 2.0 * epsilon) / eta;
          System.out.println("Ideal a2 = " + a2);
          a2 = Math.min(a2, H);
          a2 = Math.max(L, a2);
          a1 = alpha1 + (a2 - alpha2star);
          System.out.println("a1=" + a1 + ", a2=" + a2);
          // updatae alpah1, alpha2star if change greater than some eps
          if(Math.abs(a2 - alpha2star) > epsilon * (a2 + alpha2star + 1.0 + epsilon)) {
            model.setAlpha(i1, a1);
            model.setAlphaStar(i2, a2);
            alpha1     = a1;
            alpha2star = a2;
            System.out.println("case2 worked");
            changed = true;
          } else {
            System.out.println("case2: change too small: " + (a2 - alpha2star));
          }
        } else {
          System.out.println("case2: L > H");
          finnished = true;
        }
        case2 = true;
      } else if (!case3 &&
      (alpha1star > 0.0 || (alpha1     == 0.0 && deltaphi < 2.0 * epsilon)) &&
      (alpha2     > 0.0 || (alpha2star == 0.0 && deltaphi < 2.0 * epsilon))
      ) {
        L = Math.max(0.0, -gamma);
        H = Math.min(C, -gamma + C);
        System.out.println("L=" + L + "\tH=" + H);
        if(L < H) {
          a2 = alpha2 - (deltaphi + 2.0 * epsilon) / eta; // according to improvement
          //a2 = alpha2 - (deltaphi - 2.0 * epsilon) / eta; // according to tutorial
          System.out.println("Ideal a2 = " + a2);
          a2 = Math.min(a2, H);
          a2 = Math.max(L, a2);
          a1 = alpha1star + (a2 - alpha2); // according to improvement
          //a1 = alpha1star - (a2 - alpha2); // according to tutorial
          System.out.println("a1=" + a1 + ", a2=" + a2);
          // update alpha1star, alpha2 if change is greater than some eps
          if(Math.abs(a2 - alpha2) > epsilon * (a2 + alpha2 + 1.0 + epsilon)) {
            model.setAlphaStar(i1, a1);
            model.setAlpha(i2, a2);
            alpha1star = a1;
            alpha2     = a2;
            System.out.println("case3 worked");
            changed = true;
          } else {
            System.out.println("case3: change too small: " + (a2 - alpha2));
          }
        } else {
          System.out.println("case3: L > H");
          finnished = true;
        }
        case3 = true;
      } else if(!case4 &&
      (alpha1star > 0.0 || (alpha1 == 0 && deltaphi < 0.0)) &&
      (alpha2star > 0.0 || (alpha2 == 0 && deltaphi > 0.0))
      ) {
        L = Math.max(0.0, -gamma - C);
        H = Math.min(C, -gamma);
        System.out.println("L=" + L + "\tH=" + H);
        if(L < H) {
          a2 = alpha2star + deltaphi/eta;
          System.out.println("Ideal a2 = " + a2);
          a2 = Math.min(a2, H);
          a2 = Math.max(L, a2);
          a1 = alpha1star - (a2 - alpha2star);
          System.out.println("a1=" + a1 + ", a2=" + a2);
          // update alpha1star, alpha2star if change is larger than some eps
          if(Math.abs(a2 - alpha2star) > epsilon * (a2 + alpha2star + 1.0 + epsilon)) {
            model.setAlphaStar(i1, a1);
            model.setAlphaStar(i2, a2);
            alpha1star = a1;
            alpha2star = a2;
            System.out.println("case4 worked");
            changed = true;
          } else {
            System.out.println("case4: change too small: " + (a2 - alpha2star));
          }
        } else {
          System.out.println("case4: L > H");
          finnished = true;
        }
        case4 = true;
      } else {
        finnished = true;
      }
      System.out.println("!!!Errors: " + getError(i1) + "    " + getError(i2));

      // update deltaphi
      /*

      deltaphi = phi2 - phi1 + eta * ( (alpha1 - alpha1star)
                                     + (alpha1old - alpha1starold) );

      */
      deltaphi = getError(i2) - getError(i1);
      System.out.println("deltaphi=" + deltaphi);
    }
    

    // Calculate new threshold
  	double b;
    double bOld = model.getThreshold();
    System.out.println("b was " + bOld);

    /*

    if(!isBound(alpha1)) {
      b = phi1 + y1*(alpha1 - alpha1old)*k11 + y2*(alpha2 - alpha2old)*k12 + bOld;
    } else if(!isBound(alpha2)) {
      b = phi2 + y1*(alpha1 - alpha1old)*k12 + y2*(alpha2 - alpha2old)*k22 + bOld;
    } else if(!isBound(alpha1star)) {
      b = phi1 + y1*(alpha1star - alpha1starold)*k11 + y2*(alpha2star - alpha2starold)*k12 + bOld;
    } else if(!isBound(alpha2star)) {
      b = phi2 + y1*(alpha1star - alpha1starold)*k12 + y2*(alpha2star - alpha2starold)*k22 + bOld;
    } else {
      // no suitable alpha found to infer b - take middle of allowed interval
      double b1 = phi1 + y1*(alpha1 - alpha1old)*k11 + y2*(alpha2 - alpha2old)*k12 + bOld;
      double b2 = phi2 + y1*(alpha1 - alpha1old)*k12 + y2*(alpha2 - alpha2old)*k22 + bOld;
      b = (b1 + b2) / 2.0;
    }
    model.setThreshold(b);
    System.out.println("b is " + b);

    */

    // double bOld = model.getThreshold();
    if (!isBound(alpha1)) {
	b = y1 - model.internalClassify(i1) + bOld - epsilon;
    } else if (!isBound(alpha1star)) {
	b = y1 - model.internalClassify(i1) + bOld - epsilon;
    } else {
	b = y1 - model.internalClassify(i1) + bOld;
    }
    model.setThreshold(b);   
    
    // Update error cache
/*    E[i1] = 0;
    E[i2] = 0;
  	for (int l = 0; l < E.length; ++l) {
  	  if (l==i1 || l==i2) {
  		  continue;
      }
  	  if (!(isBound(model.getAlpha(l)))) {
        E[l] += y1*(a1-alpha1)*model.getKernelValue(i1, l) + y2*(a2-alpha2)*model.getKernelValue(i2, l) + bOld - b;
  	  }
  	}*/

    if(changed) {
      System.out.println("Successfuly changed things");
      return true;
    } else {
      System.out.println("Nothing changed");
      return false;
    }
  }
  
  private int examineExample(int i2) {
    double alpha2 = model.getAlpha(i2);
    double alpha2star = model.getAlphaStar(i2);
    double phi2 = getError(i2);

    /*

    if (Math.abs(phi2) < epsilon)
	return 0;

    */

    if ( (+phi2 < epsilon && alpha2star < C)   ||
         (+phi2 < epsilon && alpha2star > 0.0) ||
         (-phi2 > epsilon && alpha2 < C)       ||
         (-phi2 > epsilon && alpha2 > 0.0))
     {

	 /*

	    int secondChoice = -1;
	    double step = 0.0;
	    for (int l = 0; l < model.size(); ++l) {
        if (!isBound(model.getAlpha(l))) {
          double thisStep = Math.abs(getError(l) - phi2);
          if (thisStep > step) {
            step = thisStep;
            secondChoice = l;
          }
        }
	    }

	    if (secondChoice >= 0) {
        System.out.println("Using secondChoice heuristic for " + secondChoice + ", " + i2);
        if (takeStep(secondChoice, i2)) {
          return 1;
        }
	    }

	    int randomStart = (int) Math.floor(Math.random() * model.size());
	    for (int l = 0; l < model.size(); ++l) {
        int i1 = (l + randomStart) % model.size();
        if (!isBound(model.getAlpha(i1))) {
          System.out.println("Using unbound huristic for " + i1 + ", " + i2);
  		    if (takeStep(i1, i2)) {
            return 1;
          }
  	    }
      }
      
	 */
      
	 int randomStart = (int) Math.floor(Math.random() * model.size());

	    // The second pass should look at ALL alphas, but
	    // we've already checked the non-bound ones.
	    for (int l = 0; l < model.size(); ++l) {
        int i1 = (l + randomStart) % model.size();
        // if (isBound(model.getAlpha(i1))) {
          System.out.println("Using bound huristic for " + i1 + ", " + i2);
          if (takeStep(i1, i2)) {
            return 1;
          }
	  // }
	    }
    }
	  return 0;
  }

  private boolean isBound(double alpha) {
    return (alpha <= 0 || alpha >= C);
  }

  private double getError(int i) {
//    if (E[i] == Double.NEGATIVE_INFINITY || 
//  	    isBound(model.getAlpha(i))) {
	    E[i] = model.internalClassify(i) - target[i];
      System.out.println("Calculated error: " + E[i]);
//    }
    return E[i];
  }

  public synchronized void trainModel(SVMRegressionModel m, double[] t) {
    model = m;
    target = t;

    E = new double[model.size()];
    for (int i = 0; i < t.length; ++i) {
	    E[i] = Double.NEGATIVE_INFINITY;
    }

    int numChanged = 0;
    boolean examineAll = true;
    //int SigFig = -100;

    while (numChanged > 0 || examineAll /*|| SigFig < 3*/) {
	System.out.print(".");

	numChanged = 0;
	if (examineAll) {
	    System.out.println("Running full iteration");
	    for (int i = 0; i < model.size(); ++i) {
		numChanged += examineExample(i);
	    } 
	} else {
	    System.out.println("Running non-bounds iteration");
	    for (int i = 0; i < model.size(); ++i) {
		double alpha = model.getAlpha(i);
		if (!isBound(alpha)) {
		    numChanged += examineExample(i);
		}
	    }
	}
	
	if (examineAll) {
	    examineAll = false;
	} else {
	    examineAll = (numChanged == 0);
	}
    }

    E = null;
  }
}
