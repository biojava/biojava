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
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Samiul Hasan
 * @author Keith James
 */
public class SimpleMarkovModel
        extends
        AbstractChangeable
        implements
        MarkovModel,
        Serializable {
  public static final long serialVersionUID = -3043028839927615753l;
  private final Alphabet emissionAlpha;
  private final FiniteAlphabet stateAlpha;
  private final MagicalState magicalState;

  private final Map transFrom;
  private final Map transTo;
  private final Map transWeights;
  private transient ChangeForwarder distForwarder;

  {
    transFrom = new HashMap();
    transTo = new HashMap();
    transWeights = new HashMap();
  }

  protected ChangeSupport getChangeSupport(ChangeType ct) {
    ChangeSupport changeSupport = super.getChangeSupport(ct);

    if (
            (MarkovModel.PARAMETER.isMatchingType(ct) ||
            ct.isMatchingType(MarkovModel.PARAMETER)) &&
            distForwarder == null
    ) {
      distForwarder = new ChangeForwarder.Retyper(
              this, changeSupport, MarkovModel.PARAMETER);
      for (Iterator si = stateAlpha.iterator(); si.hasNext();) {
        State s = (State) si.next();
        if (s instanceof EmissionState) {
          EmissionState es = (EmissionState) s;
          es.addChangeListener(distForwarder, EmissionState.DISTRIBUTION);
          es.addChangeListener(ChangeListener.ALWAYS_VETO, EmissionState.ADVANCE);
        }
      }
    }

    return changeSupport;
  }

  public Alphabet emissionAlphabet() {
    return emissionAlpha;
  }

  public FiniteAlphabet stateAlphabet() {
    return stateAlpha;
  }

  public int heads() {
    return magicalState().getAdvance().length;
  }

  public MagicalState magicalState() {
    return magicalState;
  }

  public int[] advance() {
    int[] advances = new int[heads()];

    for (int i = 0; i < advances.length; i++) {
      advances[i] = 0;
    }

    for (Iterator si = stateAlphabet().iterator(); si.hasNext();) {
      State s = (State) si.next();
      if (s instanceof EmissionState) {
        int[] adv = ((EmissionState) s).getAdvance();
        for (int i = 0; i < advances.length; i++) {
          advances[i] = Math.max(advances[i], adv[i]);
        }
      }
    }

    return advances;
  }

  public Distribution getWeights(State source)
          throws IllegalSymbolException {
    stateAlphabet().validate(source);

    Distribution dist = (Distribution) transWeights.get(source);
    if (dist == null) {
      throw new BioError(
              "Model does contain " + source.getName() +
              " but the associated transition distribution is missing."
      );
    }
    return dist;
  }

  /**
   * Use this methods to customize the transition probabilities.
   * <p>
   * By default, the distribution P(destination | source) is a totally free
   * distribution. This allows the different probabilities to vary. If you
   * wish to change this behaviour (for example, to make one set of transition
   * probabilities equal to another), then use this method to replace the
   * Distribution with one of your own.
   *
   * @param source source State
   * @param dist  the new Distribution over the transition probabilites from source
   * @throws IllegalSymbolException if source is not a member of this model
   * @throws IllegalAlphabetException if dist is not a distribution over the
   *         states returned by model.transitionsFrom(source)
   * @throws ChangeVetoException if a listener vetoed this change
   */
  public void setWeights(State source, Distribution dist)
          throws IllegalSymbolException, IllegalAlphabetException, ChangeVetoException {
    FiniteAlphabet ta = transitionsFrom(source);
    if (!dist.getAlphabet().equals(ta)) {
      throw new IllegalAlphabetException(
              "Can't set distribution from state " + source.getName() +
              " as the distribution alphabet is not the alphabet of transitions: " +
              ta.getName() + " and " + dist.getAlphabet().getName()
      );
    }

    if (!hasListeners()) {
      transWeights.put(source, dist);
    } else {
      ChangeSupport cs = getChangeSupport(MarkovModel.PARAMETER);
      synchronized (cs) {
        ChangeEvent ce = new ChangeEvent(this,
                                         MarkovModel.PARAMETER,
                                         dist,
                                         transWeights.get(source));
        cs.firePreChangeEvent(ce);
        transWeights.put(source, dist);
        cs.firePostChangeEvent(ce);
      }
    }
  }

  public void createTransition(State from, State to)
          throws IllegalSymbolException, ChangeVetoException {
    stateAlphabet().validate(from);
    stateAlphabet().validate(to);

    ChangeEvent ce = new ChangeEvent(
            this, MarkovModel.ARCHITECTURE,
            new Object[]{from, to},
            null
    );

    FiniteAlphabet f = transitionsFrom(from);
    FiniteAlphabet t = transitionsTo(to);

    if (f.contains(to)) {
      throw new ChangeVetoException(
              ce,
              "Transition already exists: " + from.getName() + " -> " + to.getName()
      );
    }

    if (!hasListeners()) {
      f.addSymbol(to);
      t.addSymbol(from);
    } else {
      ChangeSupport changeSupport = getChangeSupport(MarkovModel.ARCHITECTURE);
      synchronized (changeSupport) {
        changeSupport.firePreChangeEvent(ce);

        f.addSymbol(to);
        t.addSymbol(from);

        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  public void destroyTransition(State from, State to)
          throws IllegalSymbolException, ChangeVetoException {
    stateAlphabet().validate(from);
    stateAlphabet().validate(to);

    FiniteAlphabet f = transitionsFrom(from);

    ChangeEvent ce = new ChangeEvent(
            this, MarkovModel.ARCHITECTURE,
            null,
            new Object[]{from, to}
    );

    if (!f.contains(to)) {
      throw new ChangeVetoException(
              ce,
              "Transition does not exists: " + from.getName() + " -> " + to.getName()
      );
    }

    Distribution dist = getWeights(from);
    double w = dist.getWeight(to);
    if (w != 0.0 && !Double.isNaN(w)) {
      throw new ChangeVetoException(
              ce,
              "Can't remove transition as its weight is not zero or NaN: " +
              from.getName() + " -> " + to.getName() + " = " + w
      );
    }

    if (!hasListeners()) {
      transitionsFrom(from).removeSymbol(to);
      transitionsTo(to).removeSymbol(from);
    } else {
      ChangeSupport changeSupport = getChangeSupport(MarkovModel.ARCHITECTURE);
      synchronized (changeSupport) {
        changeSupport.firePreChangeEvent(ce);

        transitionsFrom(from).removeSymbol(to);
        transitionsTo(to).removeSymbol(from);

        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  public boolean containsTransition(State from, State to)
          throws IllegalSymbolException {
    stateAlphabet().validate(to);
    return transitionsFrom(from).contains(to);
  }

  public FiniteAlphabet transitionsFrom(State from)
          throws IllegalSymbolException {
    stateAlphabet().validate(from);

    FiniteAlphabet s = (FiniteAlphabet) transFrom.get(from);
    if (s == null) {
      throw new BioError(
              "State " + from.getName() +
              " is known in states " +
              stateAlphabet().getName() +
              " but is not listed in the transFrom table"
      );
    }
    return s;
  }

  public FiniteAlphabet transitionsTo(State to)
          throws IllegalSymbolException {
    stateAlphabet().validate(to);

    FiniteAlphabet s = (FiniteAlphabet) transTo.get(to);
    if (s == null) {
      throw new BioError(
              "State " + to +
              " is known in states " +
              stateAlphabet().getName() +
              " but is not listed in the transTo table"
      );
    }
    return s;
  }

  public void addState(State toAdd)
          throws IllegalSymbolException, ChangeVetoException {
    if (toAdd instanceof MagicalState && toAdd != magicalState) {
      throw new IllegalSymbolException("Can not add a MagicalState");
    }

    if (stateAlphabet().contains(toAdd)) {
      throw new IllegalSymbolException("We already contain " + toAdd.getName());
    }

    if (toAdd instanceof EmissionState) {
      int esh = ((EmissionState) toAdd).getAdvance().length;
      if (esh != heads()) {
        throw new IllegalSymbolException(
                "This model " + stateAlphabet().getName() +
                " has " + heads() + " heads, but the state " +
                toAdd.getName() + " has " + esh + " heads"
        );
      }
    }

    if (toAdd instanceof ModelInState) {
      int esh = ((ModelInState) toAdd).getModel().advance().length;
      if (esh != heads()) {
        throw new IllegalSymbolException(
                "This model " + stateAlphabet().getName() +
                " has " + heads() + " heads, but the model-in-state " +
                toAdd.getName() + " has " + esh + " heads"
        );
      }
    }

    if (!hasListeners()) {
      doAddState(toAdd);
    } else {
      ChangeSupport changeSupport = getChangeSupport(MarkovModel.ARCHITECTURE);
      synchronized (changeSupport) {
        ChangeEvent ce = new ChangeEvent(
                this, MarkovModel.ARCHITECTURE,
                toAdd,
                null
        );
        changeSupport.firePreChangeEvent(ce);
        if (toAdd instanceof EmissionState) {
          EmissionState eState = (EmissionState) toAdd;
          eState.addChangeListener(distForwarder, EmissionState.DISTRIBUTION);
          eState.addChangeListener(ChangeListener.ALWAYS_VETO, EmissionState.ADVANCE);
        }
        doAddState(toAdd);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  private void doAddState(State toAdd)
          throws IllegalSymbolException, ChangeVetoException {
    stateAlphabet().addSymbol(toAdd);
    FiniteAlphabet fa = new SimpleAlphabet("Transitions from " + toAdd.getName());
    transFrom.put(toAdd, fa);
    transTo.put(toAdd, new SimpleAlphabet("Transitions to " + toAdd.getName()));
    transWeights.put(toAdd, new SimpleDistribution(fa));
    stateAlphabet().addSymbol(toAdd);
  }

  public void removeState(State toGo)
          throws IllegalSymbolException, IllegalTransitionException, ChangeVetoException {
    stateAlphabet().validate(toGo);
    if (toGo instanceof MagicalState) {
      throw new IllegalSymbolException("You can not remove the MagicalState");
    }
    FiniteAlphabet t;
    if ((t = transitionsFrom(toGo)).size() != 0) {
      throw new IllegalTransitionException(
              toGo, (State) t.iterator().next(),
              "You can not remove a state untill all transitions to and from it " +
              "have been destroyed"
      );
    }

    if ((t = transitionsTo(toGo)).size() != 0) {
      throw new IllegalTransitionException(
              (State) t.iterator().next(), toGo,
              "You can not remove a state untill all transitions to and from it " +
              "have been destroyed"
      );
    }

    if (!hasListeners()) {
      doRemoveState(toGo);
    } else {
      ChangeSupport changeSupport = getChangeSupport(MarkovModel.ARCHITECTURE);
      synchronized (changeSupport) {
        ChangeEvent ce = new ChangeEvent(
                this, MarkovModel.ARCHITECTURE,
                null,
                toGo
        );
        changeSupport.firePreChangeEvent(ce);
        if (toGo instanceof EmissionState) {
          EmissionState eState = (EmissionState) toGo;
          eState.removeChangeListener(distForwarder, EmissionState.DISTRIBUTION);
          eState.removeChangeListener(ChangeListener.ALWAYS_VETO, EmissionState.ADVANCE);
        }
        doRemoveState(toGo);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  private void doRemoveState(State toGo)
          throws IllegalSymbolException, ChangeVetoException
  {
    stateAlphabet().removeSymbol(toGo);
    transFrom.remove(toGo);
    transTo.remove(toGo);
  }

  public SimpleMarkovModel(int heads, Alphabet emissionAlpha, String name) {
    this(heads, emissionAlpha);
    ((SimpleAlphabet) stateAlpha).setName(name);
  }

  /**
   * @deprecated
   */
  public SimpleMarkovModel(int heads, Alphabet emissionAlpha) {
    this.emissionAlpha = emissionAlpha;
    this.stateAlpha = new SimpleAlphabet();
    this.magicalState = MagicalState.getMagicalState(emissionAlpha, heads);

    try {
      addState(magicalState);
    } catch (IllegalSymbolException ise) {
      throw new BioError("Assertion failure: Couldn't add magical state", ise);
    } catch (ChangeVetoException cve) {
      throw new BioError("Assertion failure: Couldn't add magical state", cve);
    }
  }

  public String toString() {
    return this.getClass().getName() + ":" + stateAlpha.getName();
  }
}
