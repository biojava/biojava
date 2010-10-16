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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.TranslatedDistribution;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleReversibleTranslationTable;
import org.biojava.bio.symbol.SingletonAlphabet;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.SingletonList;

/**
 * A model that guarantees to only contain emission states and dot states.
 * <p>
 * A flat model is essentialy a view onto a more complicated model that makes it
 * appear to only contain emission states and dot states. States that emit models
 * and other exotica are expressed purely in terms of these two state types.
 * <p>
 * You can train the resulting flat model, and the underlying models will be altered.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
class FlatModel
  extends
    AbstractChangeable
  implements
    MarkovModel,
    Serializable
{
  private final MarkovModel source;
  private final MarkovModel delegate;

  protected void addAState(State ourState)
  throws IllegalSymbolException {
    try {
      delegate.addState(ourState);
    } catch (ChangeVetoException cve) {
      throw new BioError("This model should be ours with no listeners", cve);
    }
  }

  public int[] advance() {
    return delegate.advance();
  }

  public Alphabet emissionAlphabet() {
    return delegate.emissionAlphabet();
  }

  public FiniteAlphabet stateAlphabet() {
    return delegate.stateAlphabet();
  }

  public int heads() {
    return delegate.heads();
  }

  public MagicalState magicalState() {
    return delegate.magicalState();
  }

  public Distribution getWeights(State source)
  throws IllegalSymbolException {
    return delegate.getWeights(source);
  }

  public boolean containsTransition(State from, State to)
  throws IllegalSymbolException {
    return delegate.containsTransition(from, to);
  }

  public FiniteAlphabet transitionsFrom(State from)
  throws IllegalSymbolException {
    return delegate.transitionsFrom(from);
  }

  public FiniteAlphabet transitionsTo(State to)
  throws IllegalSymbolException {
    return delegate.transitionsTo(to);
  }

  public void setWeights(State source, Distribution dist)
  throws ChangeVetoException {
    throw new ChangeVetoException("Can't set weights in immutable view");
  }

  public FlatModel(MarkovModel model)
  throws IllegalSymbolException, IllegalAlphabetException {
    this.source = model;
    this.delegate = new SimpleMarkovModel(
      source.advance().length,
      source.emissionAlphabet(),
      "flat"
    );

    // add all the states
    //System.out.println("Adding states");
    Map toM = new HashMap();
    Map inModel = new HashMap();
    Map misStart = new HashMap();
    Map misEnd = new HashMap();
    Map modelStart = new HashMap();
    Map modelEnd = new HashMap();

    for(Iterator i = model.stateAlphabet().iterator(); i.hasNext(); ) {
      State s = (State) i.next();
      if(s instanceof DotState) { // simple dot state in model
        DotStateWrapper dsw = new DotStateWrapper(s);
        addAState(dsw);
        inModel.put(s, model);
        toM.put(s, dsw);
        //System.out.println("Added dot state " + dsw.getName());
      } else if(s instanceof EmissionState) {  // simple emission state in model
        if(s instanceof MagicalState) {
          modelStart.put(model, model.magicalState());
          modelEnd.put(model, model.magicalState());
        } else {
          EmissionWrapper esw =
            new EmissionWrapper((EmissionState) s);
          addAState(esw);
          inModel.put(s, model);
          toM.put(s, esw);
          //System.out.println("Added emission state " + esw.getName());
        }
      } else if(s instanceof ModelInState) { // complex model inside state
        //System.out.println("Adding a model-in-state");
        ModelInState mis = (ModelInState) s;
        MarkovModel flatM = DP.flatView(mis.getModel());

        DotStateWrapper start = new DotStateWrapper(mis, "start");
        DotStateWrapper end = new DotStateWrapper(mis, "end");
        addAState(start);
        addAState(end);
        inModel.put(mis, model);
        modelStart.put(flatM, start);
        modelEnd.put(flatM, end);
        misStart.put(mis, start);
        misEnd.put(mis, end);
        //System.out.println("Added " + start.getName() + " and " + end.getName());

        for(Iterator j = flatM.stateAlphabet().iterator(); j.hasNext(); ) {
          State t = (State) j.next();
          if(t instanceof DotState) {
            DotStateWrapper dsw = new DotStateWrapper(t);
            addAState(dsw);
            inModel.put(t, flatM);
            toM.put(t, dsw);
            toM.put(((Wrapper) t).getWrapped(), dsw);
            //System.out.println("Added wrapped dot state " + dsw.getName());
          } else if(t instanceof EmissionState) {
            if(t instanceof MagicalState) {
              continue;
            }
            EmissionWrapper esw =
              new EmissionWrapper((EmissionState) t);
            addAState(esw);
            inModel.put(t, flatM);
            toM.put(t, esw);
            //toM.put(((Wrapper) t).getUnprojectedFeatures(), esw);
            //System.out.println("Added wrapped emission state " + esw.getName());
          } else { // unknown eventuality
            throw new IllegalSymbolException(s, "Don't know how to handle state: " + s.getName());
          }
        }
      } else { // unknown eventuality
        throw new IllegalSymbolException(s, "Don't know how to handle state: " + s.getName());
      }
    }

    // wire
    for(Iterator i = delegate.stateAlphabet().iterator(); i.hasNext(); ) {
      State s = (State) i.next();

      State sOrig;
      MarkovModel sModel;

      //System.out.println("Processing transitions from " + s.getName());

      // find underlying state and model for s
      if(s instanceof MagicalState) { // from magic
        sOrig = s;
        sModel = model;
      } else { // from not Magic
        Wrapper swrapper = (Wrapper) s;
        State swrapped = swrapper.getWrapped();
        MarkovModel subModel = (MarkovModel) inModel.get(swrapped);

        if(subModel != model) { // subModel -> *
          sOrig = swrapped;
          sModel = subModel;
        } else if(swrapper instanceof ModelInState) { // mis -> ?
          if(swrapper == modelStart.get(subModel)) { // mis -> subModel
            sModel = ((ModelInState) swrapped).getModel();
            sOrig = sModel.magicalState();
          } else { // mis -> model
            sModel = model;
            sOrig = swrapped;
          }
        } else { // normal
          sModel = model;
          sOrig = s;
        }
      }

      //
      // FIXME -- Matthew broked this...
      //

      TranslatedDistribution dist = null;
      // TranslatedDistribution dist = TranslatedDistribution.getDistribution(
      //  delegate.transitionsFrom(s),
      //  sModel.getWeights(sOrig)
      // );
      SimpleReversibleTranslationTable table =
        (SimpleReversibleTranslationTable) dist.getTable();
      try {
        delegate.setWeights(s, dist);
      } catch (ChangeVetoException cve) {
        throw new BioError("Couldn't edit delegate model", cve);
      }
      table.setTranslation(s, sOrig);

      // find all reachable states from s
      for(
        Iterator j = sModel.transitionsFrom(sOrig).iterator();
        j.hasNext();
      ) {
        State tOrig = (State) j.next();
        State t;
        if(tOrig instanceof MagicalState) { // * -> magic
          if(sModel == model) { // outer -> magic
            t = tOrig;
          } else { // subModel -> magic
            t = (State) modelEnd.get(sModel);
          }
        } else { // * -> normal
          t = (State) toM.get(sOrig);
        }
        table.setTranslation(t, tOrig);
      }
    }
  }

  public void createTransition(State from, State to)
  throws IllegalSymbolException, UnsupportedOperationException {
    Alphabet a = stateAlphabet();
    a.validate(from);
    a.validate(to);
    throw new UnsupportedOperationException("createTransition not supported by FlatModel");
  }

  public void destroyTransition(State from, State to)
  throws IllegalSymbolException, UnsupportedOperationException {
    Alphabet a = stateAlphabet();
    a.validate(from);
    a.validate(to);
    throw new UnsupportedOperationException("destroyTransition not supported by FlatModel");
  }

  public void addState(State toAdd)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException("addState not supported by FlatModel");
  }

  public void removeState(State toGo)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException("removeState not supported by FlatModel");
  }

  public void registerWithTrainer(ModelTrainer modelTrainer) {
    modelTrainer.registerModel(delegate);
  }

  private static class Wrapper
    extends
      AbstractChangeable
    implements
      State,
      Serializable
  {
    private final State wrapped;
    private final String extra;
    private final Alphabet matches;

    protected transient ChangeSupport changeSupport = null;

    public String getName() {
      return wrapped.getName() + "-" + extra;
    }

    public Annotation getAnnotation() {
      return wrapped.getAnnotation();
    }

    public State getWrapped() {
      return wrapped;
    }

    public Alphabet getMatches() {
      return matches;
    }

    public Set getBases() {
      return Collections.singleton(this);
    }

    public List getSymbols() {
      return new SingletonList(this);
    }

    public Wrapper(State wrapped, String extra) {
      if(wrapped == null) {
        throw new NullPointerException("Can't wrap null");
      }
      this.wrapped = wrapped;
      this.extra = extra;
      this.matches = new SingletonAlphabet(this);
    }
  }

  private static class DotStateWrapper
  extends Wrapper implements DotState {
    public DotStateWrapper(State wrapped, String extra) {
      super(wrapped, extra);
    }

    public DotStateWrapper(State wrapped)
    throws NullPointerException {
      this(wrapped, "f");
    }
  }

  private static class EmissionWrapper
  extends Wrapper implements EmissionState {
    private EmissionState getWrappedES() {
      return (EmissionState) getWrapped();
    }

    public void setAdvance(int [] advance)
    throws ChangeVetoException {
      getWrappedES().setAdvance(advance);
    }

    public int [] getAdvance() {
      return getWrappedES().getAdvance();
    }

    public Distribution getDistribution() {
      return getWrappedES().getDistribution();
    }

    public void setDistribution(Distribution dis)
    throws ChangeVetoException {
      getWrappedES().setDistribution(dis);
    }

    public void registerWithTrainer(ModelTrainer trainer) {}

    public EmissionWrapper(EmissionState wrapped) {
      this(wrapped, "-f");
    }

    public EmissionWrapper(EmissionState wrapped, String extra) {
      super(wrapped, extra);
    }
  }
}
