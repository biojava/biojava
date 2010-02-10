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
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Wraps a weight matrix up so that it appears to be a very simple HMM.
 *
 * @author Matthew Pocock
 */

public class WMAsMM
  extends
    AbstractChangeable
  implements
    MarkovModel,
    Serializable
{
  private static final int [] advance = {1};

  private final WeightMatrix wm;
  private final FiniteAlphabet stateAlpha;
  private final MagicalState magicalState;
  private final EmissionState [] states;

  private final Map transFrom;
  private final Map transTo;
  private final Map transWeights;

  private final transient ChangeForwarder distForwarder;

  public int[] advance() {
    return new int[] { 1 }; // fixme: this should be cleverer:x
  }

  public Alphabet emissionAlphabet() {
    return wm.getAlphabet();
  }

  public FiniteAlphabet stateAlphabet() {
    return stateAlpha;
  }

  public int heads() {
    return 1;
  }

  public MagicalState magicalState() {
    return magicalState;
  }

  public Distribution getWeights(State source)
  throws IllegalSymbolException {
    stateAlpha.validate(source);
    return (Distribution) transWeights.get(source);
  }

  public void setWeights(State source, Distribution dist)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "Can't replace distribution in immutable model"
    );
  }

  public FiniteAlphabet transitionsFrom(State from)
  throws IllegalSymbolException {
    Alphabet sAlpha = stateAlphabet();
    sAlpha.validate(from);

    return (FiniteAlphabet) transFrom.get(from);
  }

  public FiniteAlphabet transitionsTo(State to)
  throws IllegalSymbolException {
    Alphabet sAlpha = stateAlphabet();
    sAlpha.validate(to);

    return (FiniteAlphabet) transTo.get(to);
  }

  public void registerWithTrainer(ModelTrainer modelTrainer)
  throws BioException {
/*    for(Iterator i = stateAlphabet().iterator(); i.hasNext(); ) {
      EmissionState s = (EmissionState) i.next();
      s.registerWithTrainer(modelTrainer);
    }*/
  }

  public void createTransition(State from, State to)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "destroyTransition not supported by " + getClass());
  }

  public void destroyTransition(State from, State to)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "destroyTransition not supported by " + getClass());
  }

  public void addState(State toAdd)
  throws IllegalSymbolException, ChangeVetoException {
    if(stateAlphabet().contains(toAdd)) {
      throw new IllegalSymbolException(
        toAdd,
        "Can't add a state to a model that already contains it"
      );
    }

    throw new ChangeVetoException("addState not supported by " + getClass());
  }

  public void removeState(State toAdd)
  throws IllegalSymbolException, ChangeVetoException {
    stateAlphabet().validate(toAdd);

    throw new ChangeVetoException("removeState not supported by " + getClass());
  }

  public boolean containsTransition(State from, State to)
  throws IllegalSymbolException {
    Alphabet sAlpha = stateAlphabet();
    sAlpha.validate(from);
    sAlpha.validate(to);

    return transitionsFrom(from).contains(to);
  }

  protected int index(State s) {
    for(int i = 0; i < states.length; i++) {
      if(s == states[i]) {
        return i;
      }
    }

    return -1;
  }

  public WMAsMM(WeightMatrix wm) throws IllegalSymbolException {
    try {
      ChangeSupport changeSupport = getChangeSupport(ChangeType.UNKNOWN);
      distForwarder = new ChangeForwarder.Retyper(this, changeSupport, MarkovModel.PARAMETER);
      transFrom = new HashMap();
      transTo = new HashMap();
      transWeights = new HashMap();
      this.wm = wm;
      this.magicalState = MagicalState.getMagicalState(wm.getAlphabet(), 1);
      SimpleAlphabet sa = new SimpleAlphabet();
      sa.addSymbol(magicalState);
      this.stateAlpha = sa;
      this.states = new EmissionState[wm.columns()];
      for(int i = 0; i <= wm.columns(); i++) {
        if(i < wm.columns()) {
          sa.addSymbol(
            this.states[i] = new SimpleEmissionState(
              i + "",
              Annotation.EMPTY_ANNOTATION,
              WMAsMM.advance,
              wm.getColumn(i)
            )
          );
          wm.getColumn(i).addChangeListener(distForwarder);
        }
        State prev = (i == 0) ? magicalState : states[i-1];
        State current = (i == wm.columns()) ? magicalState : states[i];
        FiniteAlphabet fa = (FiniteAlphabet) prev.getMatches();
        transFrom.put(prev, current.getMatches());
        transTo.put(current, fa);
        Distribution dist = new UniformDistribution(fa);
        transWeights.put(prev, dist);
      }
      sa.setName("Weight Matrix columns");
    } catch (ChangeVetoException cve) {
      throw new BioError(

        "Assertion Failure: Should be able to manipulate my state alphabet.", cve
      );
    } catch (IllegalSymbolException ise) {
      throw new BioError(

        "Assertion Failure: Should be able to manipulate my state alphabet.", ise
      );
    }
  }
}
