package org.biojava3.core.sequence.template;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    Map<NucleotideCompound, List<NucleotideCompound>> equivalentsMap = new HashMap<NucleotideCompound, List<NucleotideCompound>>();

    List<NucleotideCompound> ambiguousCompounds = new ArrayList<NucleotideCompound>();
    for(NucleotideCompound compound: getAllCompounds()) {
      if (!compound.isAmbiguous()) {
        continue;
      }
      ambiguousCompounds.add(compound);
    }

    for(NucleotideCompound sourceCompound: ambiguousCompounds) {
      Set<NucleotideCompound> compoundConstituents = sourceCompound.getConstituents();
      for(NucleotideCompound targetCompound: ambiguousCompounds) {
        Set<NucleotideCompound> targetConstituents = targetCompound.getConstituents();
        if(targetConstituents.containsAll(compoundConstituents)) {
          NucleotideCompound lcSourceCompound = toLowerCase(sourceCompound);
          NucleotideCompound lcTargetCompound = toLowerCase(targetCompound);

        //equivalentsMap.put(sourceCompound, targetCompound);
    	//      equivalentsMap.put(sourceCompound, lcTargetCompound);
	        
          
          checkAdd(equivalentsMap, sourceCompound, targetCompound);
          checkAdd(equivalentsMap, sourceCompound, lcTargetCompound);
        
          checkAdd(equivalentsMap,targetCompound,sourceCompound);
          checkAdd(equivalentsMap, lcTargetCompound, sourceCompound);
          
          checkAdd(equivalentsMap, lcSourceCompound, targetCompound);
          checkAdd(equivalentsMap, lcSourceCompound, lcTargetCompound);
                    
        }
      }
    }

    //And once it's all done start adding them to the equivalents map
    
    for ( NucleotideCompound key: equivalentsMap.keySet()){
    	List<NucleotideCompound> vals = equivalentsMap.get(key);
    	for (NucleotideCompound value: vals){
    		addEquivalent((C)key,(C)value);
    		addEquivalent((C)value,(C)key);
    	}
    }
  }

  private void checkAdd(
		Map<NucleotideCompound, List<NucleotideCompound>> equivalentsMap,
		NucleotideCompound key,
		NucleotideCompound value) {

	  
      List<NucleotideCompound> listS = equivalentsMap.get(key);
      if ( listS == null){
    	  listS = new ArrayList<NucleotideCompound>();
    	  equivalentsMap.put(key, listS);
      }
      listS.add(value);
      
	
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
      for(NucleotideCompound subCompound: compound.getConstituents()) {
        settedCompounds.add(getCompoundForString(subCompound.getBase().toUpperCase()));
      }
    }
    for(NucleotideCompound compound: getAllCompounds()) {
      if(compound.getConstituents().equals(settedCompounds)) {
        return compound;
      }
    }
    return null;
  }

    /**
     * NucleotideCompounds can always complement
     */
    @Override
    public boolean isComplementable() {
        return true;
    }
}
