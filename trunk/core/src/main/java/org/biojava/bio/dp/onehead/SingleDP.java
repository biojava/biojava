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


package org.biojava.bio.dp.onehead;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.BackPointer;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPMatrix;
import org.biojava.bio.dp.DotState;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.MagicalState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleStatePath;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.DoubleAlphabet;
import org.biojava.bio.symbol.GappedSymbolList;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleGappedSymbolList;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * An implementation of DP that aligns a single sequence against a single model.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Samiul Hasan
 * @author Lukas Kall
 */
public class SingleDP extends DP implements Serializable {
  protected final HashMap emissionsProb;
  protected final HashMap emissionsOdds;
  protected final HashMap emissionsNull;

  public SingleDP(MarkovModel model)
  throws IllegalSymbolException, IllegalTransitionException, BioException {
    super(model);
    emissionsProb = new HashMap();
    emissionsOdds = new HashMap();
    emissionsNull = new HashMap();
  }

  public void update() {
    // System.out.println("Updating emissions as underlying model has changed!");
    super.update();
    // workaround for bug in vm
    if(emissionsProb != null) {
      emissionsProb.clear();
    }
    if(emissionsOdds != null) {
      emissionsOdds.clear();
    }
    if(emissionsNull != null) {
      emissionsNull.clear();
    }
  }
  
    /**
     * This method is public for the benefit of training algorithms,
     * and in the future we should look at a better way of exposing
     * the emissions cache.
     */

  public double [] getEmission(Symbol sym, ScoreType scoreType)
  throws IllegalSymbolException {
    Map emissions;
    if(scoreType == ScoreType.PROBABILITY) {
      emissions = emissionsProb;
    } else if(scoreType == ScoreType.ODDS) {
      emissions = emissionsOdds;
    } else if(scoreType == ScoreType.NULL_MODEL) {
      emissions = emissionsNull;
    } else {
      throw new BioError("Unknown ScoreType object: " + scoreType);
    }
    double [] em = (double []) emissions.get(sym);
    if(em == null) {
      int dsi = getDotStatesIndex();
      em = new double[dsi];
      State [] states = getStates();
      if(sym == AlphabetManager.getGapSymbol()) {
        em[0] = 0;
      } else {
        em[0] = Double.NEGATIVE_INFINITY;
      }
      for(int i = 1; i < dsi; i++) {
        EmissionState es = (EmissionState) states[i];
        Distribution dis = es.getDistribution();
        em[i] = Math.log(scoreType.calculateScore(dis, sym));
      }
      emissions.put(sym, em);
      /*System.out.println("Emissions for " + sym);
      for(int i = 0; i < em.length; i++) {
        System.out.println("\t" + states[i] + "\t-> " + em[i]);
      }*/
    }
    return em;
  }
  
  public double forward(SymbolList [] seq, ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalSymbolException {
    if(seq.length != 1) {
      throw new IllegalArgumentException("seq must be 1 long, not " + seq.length);
    }
    lockModel();
    DPCursor dpCursor = new SmallCursor(
      getStates(),
      seq[0],
      seq[0].iterator()
    );
    double score = forward(dpCursor, scoreType);
    unlockModel();

    return score;
  }
  
  public double backward(SymbolList [] seq, ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalSymbolException {
    if(seq.length != 1) {
      throw new IllegalArgumentException("seq must be 1 long, not " + seq.length);
    }
    lockModel();
    DPCursor dpCursor = new SmallCursor(
      getStates(),
      seq[0],
      new ReverseIterator(seq[0])
    );
    double score = backward(dpCursor, scoreType);
    unlockModel();
    
    return score;
  }

  public DPMatrix forwardMatrix(SymbolList [] seq, ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalSymbolException {
    if(seq.length != 1) {
      throw new IllegalArgumentException("seq must be 1 long, not " + seq.length);
    }
    
    lockModel();
    SingleDPMatrix matrix = new SingleDPMatrix(this, seq[0]);
    DPCursor dpCursor = new MatrixCursor(matrix, seq[0].iterator(), +1);
    matrix.setScore(forward(dpCursor, scoreType));
    unlockModel();
    
    return matrix;
  }
  
  public DPMatrix backwardMatrix(SymbolList [] seq, ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalSymbolException {
    if(seq.length != 1) {
      throw new IllegalArgumentException("seq must be 1 long, not " + seq.length);
    }
    
    lockModel();
    SingleDPMatrix matrix = new SingleDPMatrix(this, seq[0]);
    DPCursor dpCursor = new MatrixCursor(matrix, new ReverseIterator(seq[0]), -1);
    matrix.setScore(backward(dpCursor, scoreType));
    unlockModel();
    
    return matrix;
  }
  
  public DPMatrix forwardMatrix(SymbolList [] seq, DPMatrix matrix, ScoreType scoreType)
  throws IllegalArgumentException, IllegalSymbolException,
  IllegalAlphabetException, IllegalSymbolException {
    if(seq.length != 1) {
      throw new IllegalArgumentException("seq must be 1 long, not " + seq.length);
    }
    
    lockModel();
    SingleDPMatrix sm = (SingleDPMatrix) matrix;
    DPCursor dpCursor = new MatrixCursor(sm, seq[0].iterator(), +1);
    sm.setScore(forward(dpCursor, scoreType));
    unlockModel();
    
    return sm;
  }

  public DPMatrix backwardMatrix(SymbolList [] seq, DPMatrix matrix, ScoreType scoreType)
  throws IllegalArgumentException, IllegalSymbolException,
  IllegalAlphabetException, IllegalSymbolException {
    if(seq.length != 1) {
      throw new IllegalArgumentException("seq must be 1 long, not " + seq.length);
    }
    
    lockModel();
    SingleDPMatrix sm = (SingleDPMatrix) matrix;
    DPCursor dpCursor = new MatrixCursor(sm, new ReverseIterator(seq[0]), -1);
    sm.setScore(backward(dpCursor, scoreType));
    unlockModel();
    
    return sm;
  }

  protected double forward(DPCursor dpCursor, ScoreType scoreType)
  throws IllegalSymbolException {
    forward_initialize(dpCursor, scoreType);
    forward_recurse(dpCursor, scoreType);
    return forward_termination(dpCursor, scoreType);
  }

  protected double backward(DPCursor dpCursor, ScoreType scoreType)
  throws IllegalSymbolException {
    backward_initialize(dpCursor, scoreType);
    backward_recurse(dpCursor, scoreType);
    return backward_termination(dpCursor, scoreType);
  }

  protected void forward_initialize(DPCursor dpCursor, ScoreType scoreType)
    throws IllegalSymbolException {
    double [] v = dpCursor.currentCol();
    State [] states = getStates();
    
    for (int l = 0; l < getDotStatesIndex(); l++) {
      if(states[l] == getModel().magicalState()) {
        //prob 1
        v[l] = 0.0;
      } else {
        //prob 0
        v[l] = Double.NEGATIVE_INFINITY;
      }
    }
    
    int [][] transitions = getForwardTransitions();
    double [][] transitionScore = getForwardTransitionScores(scoreType);
    double [] currentCol = dpCursor.currentCol();
    //l over dots
    for (int l = getDotStatesIndex(); l < states.length; l++) {
      double score = 0.0;
      int [] tr = transitions[l];
      double [] trs = transitionScore[l];
        
      int ci = 0;
      while(
        ci < tr.length  && (
        currentCol[tr[ci]] == Double.NEGATIVE_INFINITY ||
        currentCol[tr[ci]] == Double.NaN ||
        currentCol[tr[ci]] == Double.POSITIVE_INFINITY
        )
      ) {
        ci++;
      }
      double constant = (ci < tr.length) ? currentCol[tr[ci]] : 0.0;
        
      for(int kc = 0; kc < tr.length; kc++) {
        int k = tr[kc];

        if(
          currentCol[k] != Double.NEGATIVE_INFINITY &&
          currentCol[k] != Double.NaN &&
          currentCol[k] != Double.POSITIVE_INFINITY
        ) {
          double t = trs[kc];
          score += Math.exp(t + (currentCol[k] - constant));
        } else {
        }
      }
      currentCol[l] = Math.log(score) + constant;
    }
  }

  protected void backward_initialize(DPCursor dpCursor, ScoreType scoreType)
    throws IllegalSymbolException {
    double [] v = dpCursor.currentCol();
    State [] states = getStates();

    for (int l = 0; l < states.length; l++) {
      if(states[l] == getModel().magicalState()) {
        v[l] = 0.0;
      } else {
        v[l] = Double.NEGATIVE_INFINITY;
      }
    }
  }

  private void forward_recurse(DPCursor dpCursor, ScoreType scoreType)
    throws IllegalSymbolException {
    State [] states = getStates();
    int [][] transitions = getForwardTransitions();
    double [][] transitionScore = getForwardTransitionScores(scoreType);

    // int _index = 0;
    while (dpCursor.canAdvance()) {
      // _index++;
      // System.out.println("\n*** Index=" + _index + " ***");
      dpCursor.advance();
      Symbol sym = dpCursor.currentRes();
      double [] emissions = getEmission(sym, scoreType);
//      System.out.println("Consuming " + sym.getName());
      double [] currentCol = dpCursor.currentCol();
      double [] lastCol = dpCursor.lastCol();
      for (int l = 0; l < getDotStatesIndex(); l++) { //any -> emission
        double weight = emissions[l];
        if (weight == Double.NEGATIVE_INFINITY) {
          // System.out.println("*");
          currentCol[l] = Double.NEGATIVE_INFINITY;
        } else {
          double score = 0.0;
          int [] tr = transitions[l];
          double [] trs = transitionScore[l];
          // System.out.println("l=" + states[l].getName());
          int ci = 0;
          while (
            ci < tr.length &&
            (lastCol[tr[ci]] == Double.NEGATIVE_INFINITY
            || lastCol[tr[ci]] == Double.NaN
            || lastCol[tr[ci]] == Double.POSITIVE_INFINITY)
          ) {
            ci++;
          }
          double constant = (ci < tr.length) ? lastCol[tr[ci]] : 0.0;

          for (int kc = 0; kc < tr.length; kc++) {
            int k = tr[kc];
            // System.out.println("k=" + states[k].getName());
            if (lastCol[k] != Double.NEGATIVE_INFINITY) {
              double t = trs[kc];
              
              if(states[l]== getModel().magicalState())  {
		  // System.out.print("magic " + "lastCol[k]=" + lastCol[k] + " , ");
		  // System.out.println("t=" + t);
              }
              
              score += Math.exp(t + (lastCol[k] - constant));
            } else {
              // System.out.println("-");
            }
          }
          // new_l = emission_l(sym) * sum_k(transition(k, l) * old_k)
          currentCol[l] = (weight + Math.log(score)) + constant;
          
          // System.out.println("currentCol[" + states[l].getName() + "]=" + currentCol[l]);

          if(states[l] == getModel().magicalState())  {
	      // System.out.print("magic\t");
	      //System.out.print("Weight " + weight + "\t");
              // System.out.print(", score " + score + " = " + Math.log(score) + "\t");
              // System.out.println(", constant " + constant);             
          }
        }
      }
      for(int l = getDotStatesIndex(); l < states.length; l++) { // all dot states from emissions
        double score = 0.0;
        int [] tr = transitions[l];
        double [] trs = transitionScore[l];
        
        int ci = 0;
        while(
          ci < tr.length  && (
          currentCol[tr[ci]] == Double.NEGATIVE_INFINITY
          || currentCol[tr[ci]] == Double.NaN
          || currentCol[tr[ci]] == Double.POSITIVE_INFINITY)
        ) {
          ci++;
        }
        double constant = (ci < tr.length) ? currentCol[tr[ci]] : 0.0;
        //System.out.println("constant: " + constant);
        //System.out.println("read from state: " + ((ci < tr.length) ? states[tr[ci]].getName() : "none"));
        for(int kc = 0; kc < tr.length; kc++) {
          int k = tr[kc];

          if(currentCol[k] != Double.NEGATIVE_INFINITY 
          && currentCol[k] !=Double.NaN
          && currentCol[k] != Double.POSITIVE_INFINITY) {
            double t = trs[kc];
            score += Math.exp(t + (currentCol[k] - constant));
          } else {
          }
        }
        currentCol[l] = Math.log(score) + constant;
        //System.out.println("currentCol[" + states[l].getName() + "]=" + currentCol[l]);
      }
    }
  }

  protected void backward_recurse(DPCursor dpCursor, ScoreType scoreType)
    throws IllegalSymbolException {
    State [] states = getStates();
    int stateCount = states.length;
    int [][] transitions = getBackwardTransitions();
    double [][] transitionScore = getBackwardTransitionScores(scoreType);
    double [] prevScores = new double[getDotStatesIndex()];

    while (dpCursor.canAdvance()) {
      dpCursor.advance();
      Symbol sym = dpCursor.lastRes();
      double [] emissions = getEmission(sym, scoreType);
      double [] currentCol = dpCursor.currentCol();
      double [] lastCol = dpCursor.lastCol();
      for(int k = getDotStatesIndex() - 1; k >= 0; k--) {
        prevScores[k] = emissions[k];
      }
      
//System.out.println(sym.getName());
      for (int k = stateCount-1; k >= 0; k--) {
//System.out.println("State " + k + " of " + stateCount + ", " + transitions.length);
//System.out.println(states[k].getName());
        int [] tr = transitions[k];
        double [] trs = transitionScore[k];
        double score = 0.0;
        int ci = 0;
        while (
          ci < tr.length &&
          lastCol[tr[ci]] == Double.NEGATIVE_INFINITY
        ) {
          ci++;
        }
        double constant = (ci < tr.length) ? lastCol[tr[ci]] : 0.0;
//System.out.println("Chosen constant: " + constant);
        for (int lc = tr.length-1; lc >= 0; lc--) { // any->emission
          int l = tr[lc];
          if(l >= getDotStatesIndex()) {
            continue;
          }
//System.out.println(states[k].getName() + " -> " + states[l].getName());
          double weight = prevScores[l];
//System.out.println("weight = " + weight);
          if (
            lastCol[l] != Double.NEGATIVE_INFINITY &&
            weight != Double.NEGATIVE_INFINITY
          ) {
            double t = trs[lc];
            score += Math.exp(t + weight + (lastCol[l] - constant));
          }
        }
//System.out.println("Score = " + score);
        for(int lc = tr.length-1; lc >= 0; lc--) { // any->dot
          int l = tr[lc];
          if(l < getDotStatesIndex() || l <= k) {
            break;
          }
          /*System.out.println(
            "Processing dot-state transition " +
            states[k].getName() + " -> " + states[l].getName()
          );*/
          if(currentCol[l] != Double.NEGATIVE_INFINITY) {
            score += Math.exp(trs[lc] + (currentCol[l] - constant));
          }
        }
//System.out.println("Score = " + score);
        currentCol[k] = Math.log(score) + constant;
//System.out.println("currentCol = " + currentCol[k]);
      }
    }
  }

  private double forward_termination(DPCursor dpCursor, ScoreType scoreType)
    throws IllegalSymbolException {
    double [] scores = dpCursor.currentCol();
    State [] states = getStates();

    int l = 0;
    while (states[l] != getModel().magicalState())
      l++;

    return scores[l];
  }

  protected double backward_termination(DPCursor dpCursor, ScoreType scoreType)
    throws IllegalSymbolException {
    double [] scores = dpCursor.currentCol();
    State [] states = getStates();

    int l = 0;
    while (states[l] != getModel().magicalState())
      l++;

    return scores[l];
  }
  
  public StatePath viterbi(SymbolList [] symList, ScoreType scoreType)
  throws IllegalSymbolException {
    SymbolList r = symList[0];
    DPCursor dpCursor = new SmallCursor(getStates(), r, r.iterator());
    return viterbi(dpCursor, scoreType);
  }

  private StatePath viterbi(DPCursor dpCursor, ScoreType scoreType)
  throws IllegalSymbolException {
    lockModel();
    
    State [] states = getStates();

    int [][] transitions = getForwardTransitions();
    double [][] transitionScore = getForwardTransitionScores(scoreType);
    int stateCount = states.length;

    BackPointer [] oldPointers = new BackPointer[stateCount];
    BackPointer [] newPointers = new BackPointer[stateCount];

    // initialize
    {
      double [] vc = dpCursor.currentCol();
      double [] vl = dpCursor.lastCol();
      for (int l = 0; l < getDotStatesIndex(); l++) {
        if(states[l] == getModel().magicalState()) {
          //System.out.println("Initializing start state to 0.0");
          vc[l] = vl[l] = 0.0;
          oldPointers[l] = newPointers[l] = new BackPointer(states[l]);
        } else {
          vc[l] = vl[l] = Double.NEGATIVE_INFINITY;
        }
      }
      for (int l = getDotStatesIndex(); l < stateCount; l++) {
        int [] tr = transitions[l];
        double [] trs = transitionScore[l];
        double transProb = Double.NEGATIVE_INFINITY;
        double trans = Double.NEGATIVE_INFINITY;
        int prev = -1;
        for (int kc = 0; kc < tr.length; kc++) {
          int k = tr[kc];
          double t = trs[kc];
          double s = vc[k];
          double p = t + s;
          if (p > transProb) {
            transProb = p;
            prev = k;
            trans = t;
          }
        }
        if(prev != -1) {
          vc[l] = vl[l] = transProb;
          oldPointers[l] = newPointers[l] = new BackPointer(
            states[l],
            newPointers[prev],
            trans
          );
        } else {
          vc [l] = vl[l] = Double.NEGATIVE_INFINITY;
          oldPointers[l] = newPointers[l] = null;
        }
      }          
    }

    // viterbi
    while (dpCursor.canAdvance()) { // symbol i
      dpCursor.advance();
      Symbol sym = dpCursor.currentRes();
      double [] emissions = getEmission(sym, scoreType);
      //System.out.println(sym.getName());
      double [] currentCol = dpCursor.currentCol();
      double [] lastCol = dpCursor.lastCol();
      for (int l = 0; l < states.length; l++) { // don't move from magical state
        double emission;
        if(l < getDotStatesIndex()) {
          emission = emissions[l];
        } else {
          emission = 0.0;
        }
        int [] tr = transitions[l];
        //System.out.println("Considering " + tr.length + " alternatives");
        double [] trs = transitionScore[l];
        if (emission == Double.NEGATIVE_INFINITY) {
          //System.out.println(states[l].getName() + ": impossible emission");
          currentCol[l] = Double.NEGATIVE_INFINITY;
          newPointers[l] = null;
        } else {
          double transProb = Double.NEGATIVE_INFINITY;
          double trans = Double.NEGATIVE_INFINITY;
          int prev = -1;
          for (int kc = 0; kc < tr.length; kc++) {
            int k = tr[kc];
            double t = trs[kc];
            double s = (l < getDotStatesIndex()) ? lastCol[k] : currentCol[k];
            double p = t + s;
            /*System.out.println("Looking at scores from " + states[k].getName());
            System.out.println("Old = " + lastCol[k]);
            System.out.println("New = " + currentCol[k]);
            System.out.println(
              "Considering " + states[k].getName() + " -> " +
              states[l].getName() + ", " + t + " + " + s + " = " + p
            );*/
            if (p > transProb) {
              transProb = p;
              prev = k;
              trans = t;
            }
          }
          if(prev != -1) {
            currentCol[l] = transProb + emission;
            /*System.out.println(
              states[prev].getName() + "->" + states[l].getName() + ", " +
              (trans + emission)
            );*/
            newPointers[l] = new BackPointer(
              states[l],
              (l < getDotStatesIndex()) ? oldPointers[prev] : newPointers[prev],
              trans + emission
            );
            /*System.out.println("Succesfully completed " + states[l].getName());
            System.out.println("Old = " + lastCol[l]);
            System.out.println("New = " + currentCol[l]);*/
          } else {
            //System.out.println(states[l].getName() + ": Nowhere to come from");
            currentCol[l] = Double.NEGATIVE_INFINITY;
            newPointers[l] = null;
          }
        }
      }
      
      BackPointer [] bp = newPointers;
      newPointers = oldPointers;
      oldPointers = bp;
    }

    // find max in last row
    BackPointer best = null;
    double bestScore = Double.NaN;
    for (int l = 0; l < stateCount; l++) {
      if (states[l] == getModel().magicalState()) {
        best = oldPointers[l].back;
        bestScore = dpCursor.currentCol()[l];
        break;
      }
    }

    int len = 0;
    BackPointer b2 = best;
    int dotC = 0;
    int emC = 0;
    // trace back ruit to check out size of path
    while(b2.back != b2) {
      len++;
      if(b2.state instanceof EmissionState) {
        emC++;
      } else {
        dotC++;
      }
      b2 = b2.back;
    };

    Map aMap = new HashMap();
    aMap.put(dpCursor.symList(), dpCursor.symList());
    Alignment ali = new SimpleAlignment(aMap);
    GappedSymbolList symView = new SimpleGappedSymbolList(ali);
    double [] scores = new double[len];
    List stateList = new ArrayList(len);
    for (int j = 0; j < len; j++) {
      stateList.add(null);
    }

    b2 = best;
    int ri = dpCursor.symList().length()+1;
    int lc = len;
    int gaps = 0;
    while(b2.back != b2) {
      lc--;
      //System.out.println("At " + lc + " state=" + b2.state.getName() + ", score=" + b2.score + ", back=" + b2.back);
      if(b2.state instanceof MagicalState) {
        b2 = b2.back;
        continue;
      }
      stateList.set(lc, b2.state);
      if(b2.state instanceof DotState) {
        symView.addGapInSource(ri);
        gaps++;
      } else {
        ri--;
      }
      scores[lc] = b2.score;
      b2 = b2.back;
    }

    /*System.out.println("Counted " + emC + " emissions and " + dotC + " dots");
    System.out.println("Counted backpointers. Alignment of length " + len);
    System.out.println("Counted states " + stateList.size());
    System.out.println("Input list had length " + dpCursor.symList().length());
    System.out.println("Added gaps: " + gaps);
    System.out.println("Gapped view has length " + symView.length());*/

    unlockModel();
    return new SimpleStatePath(
      bestScore,
      symView,
      new SimpleSymbolList(getModel().stateAlphabet(), stateList),
      DoubleAlphabet.fromArray(scores)
    );
  }
}
