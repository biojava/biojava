package org.biojava3.core.sequence.template;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Map.Entry;

import org.biojava3.core.sequence.compound.NucleotideCompound;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;

/**
 *
 * @author Andy Yates
 * @param <C> Type of compound this set will contain but must extend
 * NucleotideCompound
 */
public abstract class AbstractNucleotideCompoundSet<C extends NucleotideCompound>
  extends AbstractCompoundSet<C> {

  protected void addNucleotideCompound(String base, String complement, String... equivalents) {

    String[] upperEquivalents = new String[equivalents.length];
    String[] lowerEquivalents = new String[equivalents.length];
    for(int i=0; i<equivalents.length; i++) {
      upperEquivalents[i] = equivalents[i].toUpperCase();
      lowerEquivalents[i] = equivalents[i].toLowerCase();
    }

    C upper = newNucleotideCompound(base.toUpperCase(), complement.toUpperCase(), upperEquivalents);
    C lower = newNucleotideCompound(base.toLowerCase(), complement.toLowerCase(), lowerEquivalents);

    List<C> equivalentCompounds = new ArrayList<C>();

    for(int i=0; i<equivalents.length; i++) {
      equivalentCompounds.add(getCompoundForString(upperEquivalents[i]));
      equivalentCompounds.add(getCompoundForString(lowerEquivalents[i]));
    }

    addCompound(upper, lower, equivalentCompounds);
  }

  protected abstract C newNucleotideCompound(String base, String complement, String... equivalents);

  /**
   * Loops through all known nucelotides and attempts to find which are
   * equivalent to each other. Also takes into account lower casing
   * nucleotides as well as upper-cased ones.
   */
  @SuppressWarnings("unchecked")
  protected void calculateIndirectAmbiguities() {
    SetMultimap<NucleotideCompound, NucleotideCompound> equivalentsMap =
      HashMultimap.create();

    List<NucleotideCompound> ambiguousCompounds = new ArrayList<NucleotideCompound>();
    for(NucleotideCompound compound: getAllCompounds()) {
      if (!compound.isAmbiguous()) {
        continue;
      }
      ambiguousCompounds.add(compound);
    }

    for(NucleotideCompound sourceCompound: ambiguousCompounds) {
      Set<NucleotideCompound> compoundConstituents = sourceCompound.getConsituents();
      for(NucleotideCompound targetCompound: ambiguousCompounds) {
        Set<NucleotideCompound> targetConstituents = targetCompound.getConsituents();
        if(targetConstituents.containsAll(compoundConstituents)) {
          NucleotideCompound lcSourceCompound = toLowerCase(sourceCompound);
          NucleotideCompound lcTargetCompound = toLowerCase(targetCompound);

          equivalentsMap.put(sourceCompound, targetCompound);
          equivalentsMap.put(sourceCompound, lcTargetCompound);

          equivalentsMap.put(targetCompound, sourceCompound);
          equivalentsMap.put(lcTargetCompound, sourceCompound);

          equivalentsMap.put(lcSourceCompound, targetCompound);
          equivalentsMap.put(lcSourceCompound, lcTargetCompound);
        }
      }
    }

    //And once it's all done start adding them to the equivalents map
    for(Entry<NucleotideCompound, Collection<NucleotideCompound>> entry: equivalentsMap.asMap().entrySet()) {
      for(NucleotideCompound target: entry.getValue()) {
        //Just to keep it happy
        addEquivalent((C)entry.getKey(), (C)target);
        addEquivalent((C)target, (C)entry.getKey());
      }
    }
  }

  private NucleotideCompound toLowerCase(NucleotideCompound compound) {
    return getCompoundForString(compound.getBase().toLowerCase());
  }

  /**
   * Calculates the best symbol for a collection of compounds. For example
   * if you gave this method a AC it will return a M which is the ambiguity
   * symbol for these compounds.
   *
   * @param compounds Compounds to calculate ambiguity for
   * @return The ambiguity symbol which represents this set of nucleotides best
   */
  public NucleotideCompound getAmbiguity(NucleotideCompound... compounds) {
    Set<NucleotideCompound> settedCompounds = new HashSet<NucleotideCompound>();
    for(NucleotideCompound compound: compounds) {
      for(NucleotideCompound subCompound: compound.getConsituents()) {
        settedCompounds.add(getCompoundForString(subCompound.getBase().toUpperCase()));
      }
    }
    for(NucleotideCompound compound: getAllCompounds()) {
      if(compound.getConsituents().equals(settedCompounds)) {
        return compound;
      }
    }
    return null;
  }

}
