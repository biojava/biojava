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


package org.biojava.bio.dp;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.DoubleAlphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * Objects that can perform dymamic programming operations upon sequences with
 * HMMs.
 * </p>
 *
 * <p>
 * The three main DP operations are Forwards, Backwards and Viterbi. Forwards
 * and Backwards calculate the probability of the sequences having been made in
 * any way by the model. Viterbi finds the most supported way that the sequence
 * could have been made.
 * </p>
 *
 * <p>
 * Each of the functions can return the dynamic-programming matrix containing
 * the intermediate results. This may be useful for model training, or for
 * visualisation.
 * </p>
 *
 * <p>
 * Each of the funcitons can be calculated using the model probabilities, the
 * null-model probabilities or the odds (ratio between the two). For Forwards
 * and Backwards, the odds calculations produce numbers with questionable basis
 * in reality. For Viterbi with odds, you will recieve the path through the
 * model that is most different from the null model, and supported by the
 * probabilities.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public abstract class DP {
  private static List NO_ADVANCE = new ArrayList();

  private int[] getNoAdvance() {
    int heads = getModel().advance().length;
    int[] no_advance = (int[]) NO_ADVANCE.get(heads);

    if (no_advance == null) {
      no_advance = new int[heads];
      for (int i = 0; i < heads; i++) {
        no_advance[i] = 0;
      }

      NO_ADVANCE.add(heads, no_advance);
    }

    return no_advance;
  }

  /**
   * Scores the SymbolList from symbol start to symbol (start+columns) with a
   * weight matrix.
   *
   * @param matrix  the weight matrix used to evaluate the sequences
   * @param symList the SymbolList to assess
   * @param start   the index of the first symbol in the window to evaluate
   * @return  the log probability or likelyhood of this weight matrix
   *          having generated symbols start to (start + columns) of symList
   */
  public static double scoreWeightMatrix(
          WeightMatrix matrix, SymbolList symList, int start)
          throws IllegalSymbolException {
    double score = 0;
    int cols = matrix.columns();

    for (int c = 0; c < cols; c++) {
      score += Math.log(
              matrix.getColumn(c).getWeight(symList.symbolAt(c + start)));
    }

    return score;
  }

  /**
   * Scores the SymbolList from symbol start to symbol (start+columns) with a
   * weight matrix using a particular ScoreType.
   *
   * <p>
   * This method allows you to use score types such as ScoreType.ODDS. The other
   * scoreWeightMatrix methods gives a result similar or identical to
   * ScoreType.PROBABILITY.
   * </p>
   *
   * @param matrix  the weight matrix used to evaluate the sequences
   * @param symList the SymbolList to assess
   * @param scoreType the score type to apply
   * @param start   the index of the first symbol in the window to evaluate
   * @return  the sum of log scores of this weight matrix
   *          having generated symbols start to (start + columns) of symList
   * @since 1.4
   */
  public static double scoreWeightMatrix(
          WeightMatrix matrix,
          SymbolList symList,
          ScoreType scoreType,
          int start)
          throws IllegalSymbolException {
    double score = 0;
    int cols = matrix.columns();

    for (int c = 0; c < cols; c++) {
      score += Math.log(scoreType.calculateScore(
              matrix.getColumn(c), symList.symbolAt(c + start)));
    }

    return score;
  }
  public static MarkovModel flatView(MarkovModel model)
          throws IllegalAlphabetException, IllegalSymbolException {
    for (Iterator i = model.stateAlphabet().iterator(); i.hasNext();) {
      State s = (State) i.next();
      if (
              !(s instanceof DotState) &&
              !(s instanceof EmissionState)
      ) {
        return new FlatModel(model);
      }
    }

    return model;
  }

  public State[] stateList(MarkovModel mm)
          throws IllegalSymbolException, IllegalTransitionException,
          BioException {
    FiniteAlphabet alpha = mm.stateAlphabet();

    List emissionStates = new ArrayList();
    HMMOrderByTransition comp = new HMMOrderByTransition(mm);
    List dotStates = new LinkedList();
    for (Iterator addStates = alpha.iterator(); addStates.hasNext();) {
      Object state = addStates.next();
      if (state instanceof MagicalState) {
        emissionStates.add(0, state);
      } else if (state instanceof EmissionState) {
        emissionStates.add(state);
      } else {
        ListIterator checkOld = dotStates.listIterator();
        int insertPos = -1;
        while (checkOld.hasNext() && insertPos == -1) {
          Object oldState = checkOld.next();
          if (comp.compare(state, oldState) == HMMOrderByTransition.LESS_THAN) {
            insertPos = checkOld.nextIndex() - 1;
          }
        }
        if (insertPos >= 0) {
          dotStates.add(insertPos, state);
        } else {
          dotStates.add(state);
        }
      }
    }
    Collections.sort(emissionStates, new Comparator() {
      public int compare(Object o1, Object o2) {
        State s = (State) o1;
        State t = (State) o2;

        // sort by advance
        int[] sa;
        if (s instanceof EmissionState) {
          sa = ((EmissionState) s).getAdvance();
        } else {
          sa = getNoAdvance();
        }

        int[] ta;
        if (t instanceof EmissionState) {
          ta = ((EmissionState) t).getAdvance();
        } else {
          ta = getNoAdvance();
        }

        for (int i = 0; i < sa.length; i++) {
          if (sa[i] > ta[i]) {
            return -1;
          } else if (sa[i] < ta[i]) {
            return +1;
          }
        }

        // give up - sort by name
        return s.getName().compareTo(t.getName());
      }
    });
    State[] sl = new State[emissionStates.size() + dotStates.size()];
    int i = 0;
    for (Iterator si = emissionStates.iterator(); si.hasNext();) {
      EmissionState ex = (EmissionState) si.next();
      int[] ad = ex.getAdvance();
      if (ad.length != mm.advance().length) {
        throw new BioException(
                "State " + ex.getName() + " advances " + ad.length + " heads. " +
                " however, the model " + mm.stateAlphabet().getName() +
                " advances " + mm.advance().length + " heads."
        );
      }
      for (int adi = 0; ad != null && adi < ad.length; adi++) {
        if (ad[adi] != 0) {
          ad = null;
        }
      }
      if (ad != null) {
        throw new Error(
                "State " + ex.getName() + " has advance " + ad
        );
      }
      sl[i++] = ex;
    }
    for (Iterator si = dotStates.iterator(); si.hasNext();) {
      sl[i++] = (State) si.next();
    }
    return sl;
  }

  /**
   * Returns a matrix for the specified States describing all
   * valid Transitions between those States.
   * <p>
   * The matrix is 2-dimensional.  The primary array has an element
   * corresponding to every State in the states argument.  That
   * element is itself an array the elements of which identify 
   * the States that can reach that State.  The source States 
   * are identified by their index within the states [] array.
   * @param model MarkovModel to be analysed.
   * @param states The States for which the transition matrix is to be determined.
   */
  public static int[][] forwardTransitions(
          MarkovModel model,
          State[] states
          ) throws IllegalSymbolException {
    int stateCount = states.length;
    int[][] transitions = new int[stateCount][];

    for (int i = 0; i < stateCount; i++) {
      int[] tmp = new int[stateCount];
      int len = 0;
      FiniteAlphabet trans = model.transitionsTo(states[i]);
      for (int j = 0; j < stateCount; j++) {
        if (trans.contains(states[j])) {
          tmp[len++] = j;
        }
      }
      int[] tmp2 = new int[len];
      for (int j = 0; j < len; j++) {
        tmp2[j] = tmp[j];
      }
      transitions[i] = tmp2;
    }

    return transitions;
  }

  /**
   * Compute the log(score) of all transitions
   * between the specified States.  The layout
   * of the array is identical to that of the transitions
   * array.
   * <p>
   * Note that all parameters <b>MUST</b> be
   * consistent with each other!!!!
   * <p>
   * @param model The model for which the data is to be computed.
   * @param states The States within that model for which scores are required.
   * @param transitions The transition matrix obtained by calling forwardTransitions() with the above argument values.
   * @param scoreType The type of score to be evaluated.
   */
  public static double[][] forwardTransitionScores(
          MarkovModel model,
          State[] states,
          int[][] transitions,
          ScoreType scoreType
          ) {
    // System.out.println("forwardTransitionScores");
    int stateCount = states.length;
    double[][] scores = new double[stateCount][];

    for (int i = 0; i < stateCount; i++) {
      State is = states[i];
      scores[i] = new double[transitions[i].length];
      for (int j = 0; j < scores[i].length; j++) {
        try {
          scores[i][j] = Math.log(scoreType.calculateScore(
                  model.getWeights(states[transitions[i][j]]),
                  is
          ));
          /*System.out.println(
            states[transitions[i][j]] + "\t-> " +
            is.getName() + "\t = " +
            scores[i][j] + "\t(" +
            scoreType.calculateScore(
              model.getWeights(states[transitions[i][j]]),
              is
            )
          );*/
        } catch (IllegalSymbolException ite) {
          throw new BioError(
                  "Transition listed in transitions array has dissapeared.",
                  ite);
        }
      }
    }

    return scores;
  }

  public static int[][] backwardTransitions(
          MarkovModel model,
          State[] states
          ) throws IllegalSymbolException {
    int stateCount = states.length;
    int[][] transitions = new int[stateCount][];

    for (int i = 0; i < stateCount; i++) {
      int[] tmp = new int[stateCount];
      int len = 0;
      FiniteAlphabet trans = model.transitionsFrom(states[i]);
      for (int j = 0; j < stateCount; j++) {
        if (trans.contains(states[j])) {
          tmp[len++] = j;
        }
      }
      int[] tmp2 = new int[len];
      for (int j = 0; j < len; j++) {
        tmp2[j] = tmp[j];
      }
      transitions[i] = tmp2;
    }

    return transitions;
  }

  public static double[][] backwardTransitionScores(MarkovModel model,
                                                    State[] states,
                                                    int[][] transitions,
                                                    ScoreType scoreType
                                                    ) {
    int stateCount = states.length;
    double[][] scores = new double[stateCount][];

    for (int i = 0; i < stateCount; i++) {
      State is = states[i];
      scores[i] = new double[transitions[i].length];
      for (int j = 0; j < scores[i].length; j++) {
        try {
          scores[i][j] = Math.log(scoreType.calculateScore(
                  model.getWeights(is),
                  states[transitions[i][j]]
          ));
        } catch (IllegalSymbolException ite) {
          throw new BioError(
                  "Transition listed in transitions array has dissapeared",
                  ite);
        }
      }
    }

    return scores;
  }

  private MarkovModel model;
  private State[] states;
  private int[][] forwardTransitions;
  private int[][] backwardTransitions;
  private int dotStatesIndex;
  private int lockCount = 0;

  public int getDotStatesIndex() {
    return dotStatesIndex;
  }

  public MarkovModel getModel() {
    return model;
  }

  public State[] getStates() {
    return states;
  }

  public int[][] getForwardTransitions() {
    return forwardTransitions;
  }

  private Map forwardTransitionScores;
  private Map backwardTransitionScores;

  public double[][] getForwardTransitionScores(ScoreType scoreType) {
    double[][] ts = (double[][]) forwardTransitionScores.get(scoreType);
    if (ts == null) {
      forwardTransitionScores.put(scoreType, ts = forwardTransitionScores(
              getModel(), getStates(), forwardTransitions, scoreType
      ));
    }
    return ts;
  }

  public int[][] getBackwardTransitions() {
    return backwardTransitions;
  }

  public double[][] getBackwardTransitionScores(ScoreType scoreType) {
    double[][] ts = (double[][]) backwardTransitionScores.get(scoreType);
    if (ts == null) {
      backwardTransitionScores.put(scoreType, ts = backwardTransitionScores(
              getModel(), getStates(), backwardTransitions, scoreType
      ));
    }
    return ts;
  }

  public void lockModel() {
    if (lockCount++ == 0) {
      getModel().addChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
    }
  }

  public void unlockModel() {
    if (--lockCount == 0) {
      getModel().removeChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
    }
  }

  public void update() {
    try {
      if(this.states == null) {
        this.states = stateList(model);
        this.forwardTransitions = forwardTransitions(model, states);
        this.backwardTransitions = backwardTransitions(model, states);

        // Find first dot state
        int i;
        for (i = 0; i < states.length; ++i) {
          if (!(states[i] instanceof EmissionState)) {
            break;
          }
        }
        dotStatesIndex = i;
      }

      this.forwardTransitionScores.clear();
      this.backwardTransitionScores.clear();
    } catch (Exception e) {
      throw new BioError("Something is seriously wrong with the DP code", e);
    }
  }

  public DP(MarkovModel model){
    this.setModel(model);
  }
  
  /**
   * This method will result in a DP with no model. Use the setModel() method
   * to set the model before use.
   */
  public DP(){}
  
  public void setModel(MarkovModel model){
    this.model = model;
    this.forwardTransitionScores = new HashMap();
    this.backwardTransitionScores = new HashMap();
    this.update();

    model.addChangeListener(UPDATER, ChangeType.UNKNOWN);
  }

  public abstract double forward(SymbolList[] symList, ScoreType scoreType)
          throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;

  public abstract double backward(SymbolList[] symList, ScoreType scoreType)
          throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;

  public abstract DPMatrix forwardMatrix(SymbolList[] symList, ScoreType scoreType)
          throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;

  public abstract DPMatrix backwardMatrix(SymbolList[] symList, ScoreType scoreType)
          throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException;

  public abstract DPMatrix forwardMatrix(SymbolList[] symList, DPMatrix matrix, ScoreType scoreType)
          throws IllegalArgumentException, IllegalSymbolException,
          IllegalAlphabetException, IllegalTransitionException;

  public abstract DPMatrix backwardMatrix(SymbolList[] symList, DPMatrix matrix, ScoreType scoreType)
          throws IllegalArgumentException, IllegalSymbolException,
          IllegalAlphabetException, IllegalTransitionException;

  public abstract StatePath viterbi(SymbolList[] symList, ScoreType scoreType)
          throws IllegalSymbolException, IllegalArgumentException, IllegalAlphabetException, IllegalTransitionException;

  public DPMatrix forwardsBackwards(SymbolList[] symList, ScoreType scoreType)
          throws BioException {
    try {
      System.out.println("Making backward matrix");
      final DPMatrix bMatrix = backwardMatrix(symList, scoreType);
      System.out.println("Making forward matrix");
      final DPMatrix fMatrix = forwardMatrix(symList, scoreType);

      System.out.println("Making forward/backward matrix");
      return new DPMatrix() {
        public double getCell(int[] index) {
          return fMatrix.getCell(index) + bMatrix.getCell(index);
        }

        public double getScore() {
          return fMatrix.getScore();
        }

        public MarkovModel model() {
          return fMatrix.model();
        }

        public SymbolList[] symList() {
          return fMatrix.symList();
        }

        public State[] states() {
          return fMatrix.states();
        }
      };
    } catch (Exception e) {
      throw new BioException("Couldn't build forwards-backwards matrix", e);
    }
  }

  /**
   * <p>
   * Generates an alignment from a model.
   * </p>
   *
   * <p>
   * If the length is set to -1 then the model length will be sampled
   * using the model's transition to the end state. If the length is
   * fixed using length, then the transitions to the end state are implicitly
   * invoked.
   * </p>
   *
   * @param length  the length of the sequence to generate
   * @return  a StatePath generated at random
   */
  public StatePath generate(int length)
          throws IllegalSymbolException, BioException {
    List tokenList = new ArrayList();
    List stateList = new ArrayList();
    List scoreList = new ArrayList();
    double totScore = 0.0;
    double symScore = 0.0;
    int i = length;
    State oldState;
    Symbol token;

    oldState = (State) model.getWeights(model.magicalState()).sampleSymbol();
    symScore += model.getWeights(model.magicalState()).getWeight(oldState);

    DoubleAlphabet dAlpha = DoubleAlphabet.getInstance();
    if (oldState instanceof EmissionState) {
      EmissionState eState = (EmissionState) oldState;
      token = eState.getDistribution().sampleSymbol();
      symScore += eState.getDistribution().getWeight(token);
      stateList.add(oldState);
      tokenList.add(token);
      scoreList.add(dAlpha.getSymbol(symScore));
      totScore += symScore;
      symScore = 0.0;
      i--;
    }

    while (i != 0) {
      State newState = null;
      Distribution dist = model.getWeights(oldState);
      do {
        newState = (State) dist.sampleSymbol();
      } while (newState == model.magicalState() && i > 0);
      try {
        symScore += dist.getWeight(newState);
      } catch (IllegalSymbolException ise) {
        throw new BioError(
                "Transition returned from sampleTransition is invalid",
                ise);
      }

      if (newState == model.magicalState()) {
        break;
      }

      if (newState instanceof EmissionState) {
        EmissionState eState = (EmissionState) newState;
        token = eState.getDistribution().sampleSymbol();
        symScore += eState.getDistribution().getWeight(token);
        stateList.add(newState);
        tokenList.add(token);
        scoreList.add(dAlpha.getSymbol(symScore));
        totScore += symScore;
        symScore = 0.0;
        i--;
      }
      oldState = newState;
    }

    SymbolList tokens = new SimpleSymbolList(model.emissionAlphabet(), tokenList);
    SymbolList states = new SimpleSymbolList(model.stateAlphabet(), stateList);
    SymbolList scores = new SimpleSymbolList(dAlpha, scoreList);

    return new SimpleStatePath(
            totScore,
            tokens,
            states,
            scores
    );
  }

  public static class ReverseIterator implements Iterator, Serializable {
    private SymbolList sym;
    private int index;

    public ReverseIterator(SymbolList sym) {
      this.sym = sym;
      index = sym.length();
    }

    public boolean hasNext() {
      return index > 0;
    }

    public Object next() {
      return sym.symbolAt(index--);
    }

    public void remove() throws UnsupportedOperationException {
      throw new UnsupportedOperationException("This itterator can not cause modifications");
    }
  }

  private final ChangeListener UPDATER = new ChangeListener() {
    public void preChange(ChangeEvent ce)
            throws ChangeVetoException {
    }

    public void postChange(ChangeEvent ce) {
      if (ce.getType().isMatchingType(MarkovModel.ARCHITECTURE)) {
        System.out.println("architecture alterred");
        states = null;
      }

      if (
              (ce.getType().isMatchingType(MarkovModel.ARCHITECTURE)) ||
              (ce.getType().isMatchingType(MarkovModel.PARAMETER))
      ) {
        update();
      }
    }
  };

  private static class HMMOrderByTransition {
    public final static Object GREATER_THAN = new Object();
    public final static Object LESS_THAN = new Object();
    public final static Object EQUAL = new Object();
    public final static Object DISJOINT = new Object();

    private MarkovModel mm;

    private HMMOrderByTransition(MarkovModel mm) {
      this.mm = mm;
    }

    public Object compare(Object o1, Object o2)
            throws IllegalTransitionException, IllegalSymbolException {
      if (o1 == o2) {
        return EQUAL;
      }
      State s1 = (State) o1;
      State s2 = (State) o2;

      if (transitionsTo(s1, s2)) {
        return LESS_THAN;
      }
      if (transitionsTo(s2, s1)) {
        return GREATER_THAN;
      }

      return DISJOINT;
    }

    private boolean transitionsTo(State from, State to)
            throws IllegalTransitionException, IllegalSymbolException {
      Set checkedSet = new HashSet();
      Set workingSet = new HashSet();
      for (
              Iterator i = mm.transitionsFrom(from).iterator();
              i.hasNext();
              ) {
        workingSet.add(i.next());
      }

      while (workingSet.size() > 0) {
        Set newWorkingSet = new HashSet();
        for (Iterator i = workingSet.iterator(); i.hasNext();) {
          State s = (State) i.next();
          if (s instanceof EmissionState) {
            continue;
          }
          if (s == from) {
            throw new IllegalTransitionException(
                    from, from, "Loop in dot states."
            );
          }
          if (s == to) {
            return true;
          }
          for (Iterator j = mm.transitionsFrom(s).iterator(); j.hasNext();) {
            State s2 = (State) j.next();
            if (!workingSet.contains(s2) && !checkedSet.contains(s2)) {
              newWorkingSet.add(s2);
            }
          }
          checkedSet.add(s);
        }
        workingSet = newWorkingSet;
      }
      return false;
    }
  }
}

