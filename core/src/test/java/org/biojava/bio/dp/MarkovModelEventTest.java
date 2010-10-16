package org.biojava.bio.dp;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Test the event forwarding mess in MarkovModel & impls.
 *
 * @author Matthew Pocock
 */
public class MarkovModelEventTest
extends TestCase {
  public void testAddRemoveState() {
    try {
      // set up a model
      SimpleMarkovModel smm = new SimpleMarkovModel(1, DNATools.getDNA(), "add/remove test");

      EventCounter archC = new EventCounter("Architecture counter");
      EventCounter paraC = new EventCounter("Parameter counter");

      smm.addChangeListener(archC, MarkovModel.ARCHITECTURE);
      smm.addChangeListener(paraC, MarkovModel.PARAMETER);


      // make some silly states
      DotState ds = new SimpleDotState("d1");
      EmissionState es = new SimpleEmissionState(
              "e1",
              Annotation.EMPTY_ANNOTATION,
              new int[]{1},
              DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA()));


      // add ds & check events
      smm.addState(ds);
      assertEquals("All architecture events allowed: " + archC, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("One architecture event: " + archC, 1, archC.getPostCounts());
      assertEquals("No parameter events: " + paraC, 0, paraC.getPostCounts());

      archC.zeroCounts();
      paraC.zeroCounts();

      // add es & check events
      smm.addState(es);
      assertEquals("All architecture events allowed: " + archC, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("One architecture event: " + archC, 1, archC.getPostCounts());
      assertEquals("No parameter events: " + paraC, 0, paraC.getPostCounts());

      archC.zeroCounts();
      paraC.zeroCounts();

      // create transition
      smm.createTransition(ds, es);
      assertEquals("All architecture events allowed: " + archC, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("One architecture event: " + archC, 1, archC.getPostCounts());
      assertEquals("No parameter events: " + paraC, 0, paraC.getPostCounts());

      archC.zeroCounts();
      paraC.zeroCounts();

      // delete transition
      smm.destroyTransition(ds, es);
      assertEquals("All architecture events allowed: " + archC, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("One architecture event: " + archC, 1, archC.getPostCounts());
      assertEquals("No parameter events: " + paraC, 0, paraC.getPostCounts());

      archC.zeroCounts();
      paraC.zeroCounts();

      // remove ds
      smm.removeState(ds);
      smm.removeState(es);
    } catch (IllegalAlphabetException iae) {
      throw (AssertionError) new AssertionError("Could not create distributions.").initCause(iae);
    } catch (ChangeVetoException cve) {
      throw (AssertionError) new AssertionError("Could not modify model.").initCause(cve);
    } catch (IllegalSymbolException ise) {
      throw (AssertionError) new AssertionError("Could not modify model.").initCause(ise);
    } catch (IllegalTransitionException ite) {
      throw (AssertionError) new AssertionError("Could not modify model.").initCause(ite);
    }
  }

  public void testSetWeights() {
    try {
      // set up a model
      SimpleMarkovModel smm = new SimpleMarkovModel(1, DNATools.getDNA(), "add/remove test");

      EventCounter everything = new EventCounter("Everything");
      EventCounter archC = new EventCounter("Architecture counter");
      EventCounter paraC = new EventCounter("Parameter counter");

      smm.addChangeListener(everything, ChangeType.UNKNOWN);
      smm.addChangeListener(archC, MarkovModel.ARCHITECTURE);
      smm.addChangeListener(paraC, MarkovModel.PARAMETER);

      // make some silly states
      DotState ds = new SimpleDotState("d1");
      EmissionState es = new SimpleEmissionState(
              "e1",
              Annotation.EMPTY_ANNOTATION,
              new int[]{1},
              DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA()));

      smm.addState(ds);
      smm.addState(es);
      smm.createTransition(ds, es);

      everything.zeroCounts();
      archC.zeroCounts();
      paraC.zeroCounts();

      // set the transition weight distribution
      smm.setWeights(ds, new UniformDistribution((FiniteAlphabet) smm.getWeights(ds).getAlphabet()));
      assertEquals("All architecture events allowed: " + archC + "\n\t" + everything, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC + "\n\t" + everything, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("No architecture event: " + archC + "\n\t" + everything, 0, archC.getPostCounts());
      assertEquals("One parameter events: " + paraC + "\n\t" + everything, 1, paraC.getPostCounts());

      everything.zeroCounts();
      archC.zeroCounts();
      paraC.zeroCounts();

      // set the distribution associated with es
      Distribution dist = DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
      es.setDistribution(dist);
      assertEquals("All architecture events allowed: " + archC + "\n\t" + everything, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC + "\n\t" + everything, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("No architecture event: " + archC + "\n\t" + everything, 0, archC.getPostCounts());
      assertEquals("One parameter events: " + paraC + "\n\t" + everything, 1, paraC.getPostCounts());

      everything.zeroCounts();
      archC.zeroCounts();
      paraC.zeroCounts();

      // set two individual weights in dist
      dist.setWeight(DNATools.a(), 0.9);
      dist.setWeight(DNATools.t(), 0.1);
      assertEquals("All architecture events allowed: " + archC + "\n\t" + everything, archC.getPreCounts(), archC.getPostCounts());
      assertEquals("All parameter events allowed: " + paraC + "\n\t" + everything, paraC.getPreCounts(), paraC.getPostCounts());
      assertEquals("No architecture event: " + archC + "\n\t" + everything, 0, archC.getPostCounts());
      assertEquals("Two parameter events: " + paraC + "\n\t" + everything, 2, paraC.getPostCounts());
    } catch (IllegalAlphabetException iae) {
      throw (AssertionError) new AssertionError("Could not create distributions.").initCause(iae);
    } catch (ChangeVetoException cve) {
      throw (AssertionError) new AssertionError("Could not modify model.").initCause(cve);
    } catch (IllegalSymbolException ise) {
      throw (AssertionError) new AssertionError("Could not modify model.").initCause(ise);
    }
  }


}
