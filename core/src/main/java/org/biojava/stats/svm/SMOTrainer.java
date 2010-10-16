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

import java.util.Iterator;
import java.util.Set;

/**
 * Train a support vector machine using the Sequential Minimal
 * Optimization algorithm.  See Kernel Methods book.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class SMOTrainer {
  private double _C = 1000;
  private double _epsilon = 0.000001;

  public void setC(double C) {
    this._C = C;
  }
    
  public double getC() {
    return _C;
  }

  public void setEpsilon(double epsilon) {
    this._epsilon = epsilon;
  }
    
  public double getEpsilon() {
    return _epsilon;
  }

  private boolean takeStep(SMOTrainingContext trainingContext, int i1, int i2) {
    // //System.out.print("+");

    if (i1 == i2) {
	    return false;
    }
    
    double y1 = trainingContext.getTarget(i1);
    double y2 = trainingContext.getTarget(i2);
    double alpha1 = trainingContext.getAlpha(i1);
    double alpha2 = trainingContext.getAlpha(i2);
    double E1 = trainingContext.getError(i1);
    double E2 = trainingContext.getError(i2);
    double s = y1 * y2;
    double C = trainingContext.getC();
    double epsilon = trainingContext.getEpsilon();
	
  	double L, H;
    if (y2 != y1) /* preferred (s<0) */ {
	    // targets in opposite directions
	    L = Math.max(0, alpha2 - alpha1);
	    H = Math.min(C, C + alpha2 - alpha1);
    } else {
	    // Equal targets.
	    L = Math.max(0, alpha1 + alpha2 - C);
	    H = Math.min(C, alpha1 + alpha2);
    }
    if (L == H) {
      ////System.out.print("h");
      return false;
    }

    double k11 = trainingContext.getKernelValue(i1, i1);
    double k12 = trainingContext.getKernelValue(i1, i2);
    double k22 = trainingContext.getKernelValue(i2, i2);
    double eta = 2 * k12 - k11 - k22;
	
  	double a1 = 0, a2 = 0;
    if (eta > 0 && eta < epsilon) {
      eta = 0.0;
    }
  
  	if (eta < 0) {
	    a2 = alpha2 - y2 * (E1 - E2) / eta;
	    if (a2 < L) {
        a2 = L;
      } else if (a2 > H) {
        a2 = H;
      }
    } else {
      //System.out.println("Positive eta!");

      /*

      double gamma = alpha1 + s*alpha2;
      double v1 = model.classify(model.getVector(i1)) + model.getThreshold() - y1*alpha1*k11 - y2*alpha2*k12;
      double v2 = model.classify(model.getVector(i2)) + model.getThreshold() - y1*alpha1*k12 - y2*alpha2*k22;

      double Lobj = gamma - s * L + L - 0.5*k11*Math.pow(gamma - s*L,2) - 0.5*k22*Math.pow(L,2) - s*k12*(gamma-s*L)*L-y1*(gamma-s*L) - y1*(gamma - s*L)*v1 - y2*L*v2; 
      double Hobj = gamma - s * H + H - 0.5*k11*Math.pow(gamma - s*H,2) - 0.5*k22*Math.pow(H,2) - s*k12*(gamma-s*H)*H-y1*(gamma-s*H) - y1*(gamma - s*H)*v1 - y2*H*v2;
      if (Lobj > Hobj+epsilon) {
        a2 = L;
      } else if (Lobj < Hobj-epsilon) {
        a2 = H;
      } else {
        a2 = alpha2;
      }
      */
      ////System.out.print("+");
      return false;
    }
	
  	a1 = alpha1 + s*(alpha2 - a2);
    if (Math.abs(a1 - alpha1) < epsilon * (a1 + alpha1+1 +epsilon)) {
      //    //System.out.print("s");
      return false;
    }

    // Calculate new threshold
	
  	double b;
    double bOLD = trainingContext.getThreshold();

  	if (0 < a1 && a1 < C) {
      // use "b1 formula"
	    // //System.out.println("b1");
  	  b = E1 + y1*(a1 - alpha1)*k11 + y2*(a2 - alpha2)*k12 + bOLD;
  	} else if (0 < a2 && a2 < C) {
  	  // use "b2 formula"
  	  b = E2 + y1*(a1 - alpha1)*k12 + y2*(a2 - alpha2)*k22 + bOLD;
	    // //System.out.println("b2");
  	} else {
	    // Both are at bounds -- use `half way' method.
	    double b1, b2;
	    b1 = E1 + y1*(a1 - alpha1)*k11 + y2*(a2 - alpha2)*k12 + bOLD;
	    b2 = E2 + y1*(a1 - alpha1)*k12 + y2*(a2 - alpha2)*k22 + bOLD;
	    // //System.out.println("hybrid");
	    b = (b1 + b2) / 2.0;
    }
    trainingContext.setThreshold(b);
    trainingContext.setAlpha(i1, a1);
    trainingContext.setAlpha(i2, a2);

    // Update error cache

    trainingContext.resetError(i1);
    trainingContext.resetError(i2);

    for (int l = 0; l < trainingContext.size(); ++l) {
      if (l==i1 || l==i2) {
        continue;
      }
      if (!trainingContext.isBound(
        trainingContext.getAlpha(l)
      )) {
        trainingContext.updateError(
          l,
          y1*(a1-alpha1)*trainingContext.getKernelValue(i1, l) +
          y2*(a2-alpha2)*trainingContext.getKernelValue(i2, l) +
          bOLD - b
        );
      }
    }

    return true;
  }

  private int examineExample(SMOTrainingContext trainingContext, int i2) {
    double y2 = trainingContext.getTarget(i2);
    double alpha2 = trainingContext.getAlpha(i2);
    double E2 = trainingContext.getError(i2);
    double r2 = E2 * y2;
    double epsilon = trainingContext.getEpsilon();
    double C = trainingContext.getC();

    //System.out.println("r2 = " + r2);
    //System.out.println("alpha2 = " + alpha2);
    //System.out.println("epsilon = " + epsilon);
    //System.out.println("C = " + C);
    if ((r2 < -epsilon && alpha2 < C) || (r2 > epsilon && alpha2 > 0)) {
	    int secondChoice = -1;
	    double step = 0.0;
      //System.out.println("First choice heuristic");
	    for (int l = 0; l < trainingContext.size(); ++l) {
        if (!trainingContext.isBound(
          trainingContext.getAlpha(l)
        )) {
          double thisStep = Math.abs(trainingContext.getError(l) - E2);
          if (thisStep > step) {
            step = thisStep;
            secondChoice = l;
          }
        }
	    }

	    if (secondChoice >= 0) {
        if (takeStep(trainingContext, secondChoice, i2)) {
          return 1;
        }
	    }

      //System.out.println("Unbound");
  	  int randomStart = (int) Math.floor(Math.random() * trainingContext.size());
  	  for (int l = 0; l < trainingContext.size(); ++l) {
        int i1 = (l + randomStart) % trainingContext.size();
        if (!trainingContext.isBound(
          trainingContext.getAlpha(i1)
        )) {
  		    if (takeStep(trainingContext, i1, i2)) {
            return 1;
          }
        }
  	  }
	    // The second pass should look at ALL alphas, but
	    // we've already checked the non-bound ones.
      //System.out.println("Bound");
	    for (int l = 0; l < trainingContext.size(); l++) {
        int i1 = (l + randomStart) % trainingContext.size();
        if (trainingContext.isBound(
          trainingContext.getAlpha(i1)
        )) {
          if (takeStep(trainingContext, i1, i2)) {
            return 1;
          }
        }
	    }
    } else {
      //System.out.print("Nothing to optimize");
    }
    return 0;
  }

  public SVMClassifierModel
  trainModel(SVMTarget target, SVMKernel kernel, TrainingListener l) {
    SMOTrainingContext trainingContext = new SMOTrainingContext(target, kernel, l);

    int numChanged = 0;
    boolean examineAll = true;
    
    while (numChanged > 0 || examineAll) {
      numChanged = 0;
	    if (examineAll) {
        //System.out.println("Running full iteration");
        for(int i = 0; i < trainingContext.size(); i++) {
          //System.out.println("Item " + i);
          numChanged += examineExample(trainingContext, i);
        } 
	    } else {
        //System.out.println("Running non-bounds iteration");
        for(int i = 0; i < trainingContext.size(); i++) {
          double alpha = trainingContext.getAlpha(i);
          if (!trainingContext.isBound(alpha)) {
            numChanged += examineExample(trainingContext, i);
          }
        }
      }
      if (examineAll) {
        examineAll = false;
      } else {
        examineAll = (numChanged == 0);
      }
	    
  	  trainingContext.trainingCycleCompleted();
    }
    trainingContext.trainingCompleted();
    
    return trainingContext.getModel();
  }


  final class SMOTrainingContext implements TrainingContext {
    private double C;
    private double epsilon;
    private TrainingListener listener;
    private int cycle = 0;
    private TrainingEvent ourEvent;
    private SVMTarget target;
    private SVMClassifierModel model;

    private Object [] items;
    private double [] alphas;
    private double [] targets;
    private double [] E;

    private boolean isBound(double alpha) {
      return (alpha <= 0 || alpha >= getC());
    }
    
    public int size() {
      return items.length;
    }
    
    public int getCurrentCycle() {
      return cycle;
    }
    
    public void trainingCycleCompleted() {
      cycle++;
      if(listener != null) {
        listener.trainingCycleComplete(ourEvent);
      }
    }

    public void trainingCompleted() {
      for(int i = 0; i < size(); i++) {
        if(getAlpha(i) == 0) {
          model.removeItem(getItem(i));
        }
      }
      if (listener != null) {
        listener.trainingComplete(ourEvent);
      }
    }

    public Object getItem(int i) {
      return items[i];
    }
    
    public double getAlpha(int i) {
      return alphas[i];
    }
    
    public void setAlpha(int i, double a) {
      alphas[i] = a;
      model.setAlpha(getItem(i), getAlpha(i) * getTarget(i));
    }
    
    public double getTarget(int i) {
      return targets[i];
    }
    
    public double getC() {
      return C;
    }
    
    public double getEpsilon() {
      return epsilon;
    }
    
    public SVMTarget getTarget() {
      return target;
    }
    
    public SVMClassifierModel getModel() {
      return model;
    }
    
    public void setThreshold(double t) {
      model.setThreshold(t);
    }
    
    public double getThreshold() {
      return model.getThreshold();
    }
    
    public double getError(int i) {
      double alpha = getAlpha(i);
      if (isBound(alpha)) {
        return E[i] = getModel().classify(getItem(i)) - getTarget(i);
      }
      return E[i];
    }
   
    public void updateError(int i, double delta) {
      E[i] += delta;
    }

    public void resetError(int i) {
      E[i] = getModel().classify(getItem(i)) - getTarget(i);
    }

    public double getKernelValue(int i1, int i2) {
      return getModel().getKernel().evaluate(getItem(i1), getItem(i2));
    }
    
    public SMOTrainingContext(SVMTarget target, SVMKernel kernel, TrainingListener l) {
      C = SMOTrainer.this.getC();
      epsilon = SMOTrainer.this.getEpsilon();
      model = new SimpleSVMClassifierModel(kernel, target);
      model.setThreshold(0.0);
      listener = l;
      ourEvent = new TrainingEvent(this);
      cycle = 0;
      Set itemSet = target.items();
      int size = itemSet.size();
      items = new Object[size];
      alphas = new double[size];
      targets = new double[size];
      E = new double[size];
      Iterator itemI = itemSet.iterator();
      for (int i = 0; itemI.hasNext(); ++i) {
        Object item = itemI.next();
        items[i] = item;
        targets[i] = target.getTarget(item);
        alphas[i] = model.getAlpha(item) / targets[i];
        E[i] = - targets[i];
      }
    }
  }
}
