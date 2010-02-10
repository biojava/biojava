package org.biojava3.core.sequence.template;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 *
 * @author Andy Yates
 * @param <C> Type of compound this set will contain but must extend
 * NucleotideCompound
 */
public abstract class AbstractNucleotideCompoundSet<C extends NucleotideCompound>
  extends AbstractCompoundSet<C> {

  protected void addNucleotideCompound(String base, String complement, String... equivalents) {
    C upper = newNucleotideCompound(base.toUpperCase(), complement.toUpperCase());
    C lower = newNucleotideCompound(base.toLowerCase(), complement.toLowerCase());

    List<C> equivalentCompounds = new ArrayList<C>();

    for(int i=0; i<equivalents.length; i++) {
      equivalentCompounds.add(getCompoundForString(equivalents[i]));
      equivalentCompounds.add(getCompoundForString(equivalents[i].toLowerCase()));
    }

    addCompound(upper, lower, equivalentCompounds);
  }

  protected abstract C newNucleotideCompound(String base, String complement);

}
