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

package org.biojava.bio.seq.impl;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.AbstractSymbolList;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Support class for applications which need to patch together sections
 * of sequence into a single SymbolList.  This class isn't intended
 * for direct use in user code -- instead, it is a helper for people
 * implementing the full BioJava assembly model.  See SimpleAssembly
 * for an example.
 *
 * @author Thomas Down
 * @author Greg Cox
 * @author Matthew Pocock
 * @author David Huen (support for overlapped features)
 * @since 1.1
 */

public class NewAssembledSymbolList extends AbstractSymbolList {
    private boolean autoLength = true;
    private int length = 0;

    private SortedMap components; // maps between locations and the actual components
    private List componentList; // an ordered list of component locations in assembly space

    private final Symbol noninformativeSymbol = DNATools.n();


    private class TranslatedSymbolList extends AbstractSymbolList
    {
        ComponentFeature cf;
        int translation;
        int length;
        SymbolList underlyingSymList;

        private TranslatedSymbolList(ComponentFeature cf)
//            throws IllegalAlphabetException
        {
//            if (cf.getComponentSequence().getAlphabet() != DNATools.getDNA()) 
//                throw new IllegalAlphabetException("TranslatedSymbolList only accepts DNA.");

            this.cf = cf;

            if (cf.getStrand() == StrandedFeature.POSITIVE) {
                translation = cf.getLocation().getMin() - cf.getComponentLocation().getMin();
            }
            else {
                translation = cf.getLocation().getMax() + cf.getComponentLocation().getMin();
            }

            // cache these for speed
            underlyingSymList = cf.getComponentSequence();
            length = underlyingSymList.length();
        }

        /**
         * this is the only method that needs implementation. (need to reverse complement?)
         */
        public Symbol symbolAt(int i)
        {
            // compute index
            int idx;

            try {
                if (cf.getStrand() == StrandedFeature.POSITIVE) {
                    idx = i - translation;
                    if ((idx<1) || (idx>length)) return noninformativeSymbol;
//                    System.out.println("##########+++ " + i + " " + idx);
                    return underlyingSymList.symbolAt(idx);                
                }
                else {
                    idx = translation - i;
                    if ((idx<1) || (idx>length)) return noninformativeSymbol;
//                    System.out.println("##########--- " + i + " " + idx);
                    return DNATools.complement(underlyingSymList.symbolAt(idx));
                }
            }
            catch (IllegalSymbolException ise) {
                return noninformativeSymbol;
            }
        }


        /**
         * should I bother with this at all?
         */
        public Alphabet getAlphabet() { return DNATools.getDNA(); }

        /**
         * this one has to be faked to see much longer
         * than it is so the translated locations appear valid.
         */
        public int length() { return cf.getLocation().getMax(); }

        public ComponentFeature getFeature() { return cf; }
    }


    {
        components = new TreeMap(Location.naturalOrder);
        componentList = new ArrayList();
    }

    public void setLength(int len) {
	autoLength = false;
	length = len;
    }

    public void putComponent(ComponentFeature f) {
       // when receiving a new entry there are
       // two passes.

       // in pass one, we remove any redundant entries from the
       // existing components map and modify the new entry to deal
       // with existing entries.

       // in pass two, we recreate componentList from the comonents
       // key.

       // components has as key the assembly sequence ranges
       // that each component maps.  The value is the component
       // SymbolList itself.

       // componentList stores the ranges as Locations to accelerate
       // searches for the current component.

       // PASS 1
       // ======
       // when entering a new component, the following checks
       // are needed:-
       // 1) if the new range is completely contained by
       //    a current range, discard it and its component.
       // 2) if the new range completely contains a current
       //    range, discard the current range and its component.
       // 3) if the new range overlaps an existing range, adjust
       //    the new range to remove the overlap.


       // find the first range  that contains/right of the
       // left limit of the current range
       Location newLoc = new RangeLocation(f.getLocation().getMin(), f.getLocation().getMax());
       int startIdx = idxRightOfPoint(f.getLocation().getMin());

//       System.out.println("testing " + f.getLocation() + " result is " + startIdx);

       if (startIdx == NOTFOUND) {
           // either there are no entries yet or that all of 
           // them are to left of this point.
           // in this case, we just insert the new entry
//           System.out.println("putting " + f.getLocation() + " " + f.getLocation() + " " + f.getStrand() + " " + f.getComponentSequence().seqString());
           components.put(f.getLocation(), new TranslatedSymbolList(f));
           recreateList();
           return;
       }
       else {
           // we are going to reconstruct the mappings
           // according to above rules

           for (int i=startIdx;
               i < componentList.size(); // check for end of componentList
               i++) {

               Location thisLoc = (Location) componentList.get(i);
//               System.out.println("comparing " + newLoc + " against " + thisLoc);

               // check if there are any other relevant entries
               if (newLoc.getMin() > thisLoc.getMax()) break;

               if (newLoc.contains(thisLoc)) {
                   // eliminate this entry from both componentList
                   // and from the components
                   // old mapping can also be discarded
                   components.remove(thisLoc);
               }
               else if (thisLoc.contains(newLoc)) {
                   // throw away new entry
                   return;
               }
               else if (newLoc.overlaps(thisLoc)) {
                   // resolve the overlap
                   if (newLoc.getMin() < thisLoc.getMin()) {
                       // newLoc on left: correct right limit
                       newLoc = new RangeLocation(newLoc.getMin(), thisLoc.getMin() - 1);
                   }
                   else {
                       // right overlap
                       newLoc = new RangeLocation(thisLoc.getMax() + 1, newLoc.getMax());
                   }
               }
               else {
                   // no overlap of any kind at all
                   // this can only mean we can slot in the entry here
//                   System.out.println("putting " + newLoc + " " + f.getLocation() + " " + f.getStrand() + " " + f.getComponentSequence().seqString());
                   components.put(newLoc, new TranslatedSymbolList(f));
                   recreateList();
                   return;
               }
           }

           // you can get here if you have "consumed" the last location
           // or overlapped it.  It only remains to put in the new
           // entry.  In this case, it is either going to be last
           // (because you have consumed the last, or second last
           // (because you have overlapped the last).
//           System.out.println("putting " + newLoc + " " + f.getLocation() + " " + f.getStrand() + " " + f.getComponentSequence().seqString());
           components.put(newLoc, new TranslatedSymbolList(f));
           recreateList();
           return;
       }

    }

    private void recreateList()
    {
        componentList.clear();
        componentList.addAll(components.keySet());
    }

    public void removeComponent(ComponentFeature f) {
        // search thru current components to
        // remove the specified component
        for (int i=0; i< componentList.size(); i++) {
            Location loc = (Location) componentList.get(i);
            TranslatedSymbolList sl = (TranslatedSymbolList) components.get(loc);

            if (sl.getFeature() != f) continue;

            // have a match, remove from both componentList
            // and component
            // loc is now the mapping range of the to-be-removed feature
            componentList.remove(i);

            // at this point the former left and right neighbours of the
            // removed feature are i-1 and i.

            // check bounds of both neighbouring locations
            // and extend if possible to cover the void
            // left by the removal
            // extension would not occur if the neighbours
            // are not abutting to begin with
            int minExtensionFromRight = loc.getMin();

            // the way in which the locations are assembled
            // mean that at most, the left and right neighbours can be extended
            // to some point within the removed range.
            if (i !=0 ) {
                // there was a left neighbour, extend it if poss
                Location leftLoc = (Location) componentList.get(i-1);
                Location leftCFLoc = (Location) components.get(leftLoc);

                if (leftCFLoc.getMax() > loc.getMin()) {
                    // there is scope for extension
                    // compute the right limit
                    int newRightLimitOnLeftRange = Math.min(leftCFLoc.getMax(), loc.getMax());

                    Location newLeftLoc = new RangeLocation(
                        leftLoc.getMin(),
                        newRightLimitOnLeftRange
                        );

                    // update all tables
                    componentList.remove(i-1);
                    componentList.add(i-1, newLeftLoc);
                    Object value = components.remove(leftLoc);
                    components.put(newLeftLoc, value);

                    minExtensionFromRight = newRightLimitOnLeftRange + 1;
                }
            }

            if (i != componentList.size()) { // we use i because that's the index of the right neighbour after deletion
                // there was a right neighbour
                Location rightLoc = (Location) componentList.get(i);
                Location rightCFLoc = (Location) components.get(rightLoc);

                if ((loc.getMax() > rightCFLoc.getMin()) && // did the removed range limit the right neighbour?
                    (minExtensionFromRight < rightLoc.getMin()) ) { // is there any potential extension?
                    
                    // there is potential for extension
                    int newLeftLimitOnRightRange = Math.max(minExtensionFromRight, rightCFLoc.getMin());

                    Location newRightLoc = new RangeLocation(
                        newLeftLimitOnRightRange,
                        rightLoc.getMax()
                        );

                    // update all tables
                    componentList.remove(i);
                    componentList.add(i, newRightLoc);
                    Object value = components.remove(rightLoc);
                    components.put(newRightLoc, value);

                }
            }
        }
    }

    private SymbolList getComponentSymbols(Location loc) {
        // the only kind of object in this class is a TranslatedSymbolList
        // this already handles the coordinate transformations
        // I only need to select the right one
        return (SymbolList) components.get(loc);
    }


    public Set getComponentLocationSet() {
//	return components.keySet();
        TreeSet newSet = new TreeSet(Location.naturalOrder);
        newSet.addAll(componentList);
        return newSet;
    }

    /**
     * Find the location containing p in a sorted list of non-overlapping contiguous
     * locations.
     */

    private Location lastLocation = Location.empty;

    // searches for first range that contains the point.
    // returns null on failure
    private Location locationOfPoint(int p) {
	if (lastLocation.contains(p)) {
	    return lastLocation;
	}

	int first = 0;
	int last = componentList.size() - 1;

	while (first <= last) {
	    int check = (first + last) / 2;
	    Location checkL = (Location) componentList.get(check);
	    if (checkL.contains(p)) {
		lastLocation = checkL;
		return checkL;
	    }

	    if (p < checkL.getMin()) {
		last = check - 1;
	    } else {
		first = check + 1;
	    }
	}

	return null;
    }

    private int NOTFOUND = -1;

    private int idxRightOfPoint(int p)
    {

        int first = 0;
        int last = componentList.size() - 1;

        int check = 0;
        Location checkL = null;
        while (first <= last) {
            check = (first + last) / 2;
            checkL = (Location) componentList.get(check);
            if (checkL.contains(p)) {
//                System.out.println("idxRightOfPoint("+p+") returning " + check + " on " + checkL);
                return check;
            }

            if (p < checkL.getMin())
                last = check - 1;
            else
                first = check + 1;
        }

        try {
            if (p < checkL.getMin()) {
                return check;
            } else {
                int nextIdx = check+1;
                if ((nextIdx < componentList.size())
                    && (p < ((Location) componentList.get(nextIdx)).getMin()) )
                    return nextIdx;
                else
                    return NOTFOUND;
            }
        }
        catch (NullPointerException npe) {
            // this deals with the start
            // when there are no entries
            return NOTFOUND;
        }
    }

    public Alphabet getAlphabet() {
	return DNATools.getDNA();
    }

    public int length() {
      if(autoLength) {
        int componentCount = componentList.size();

        if (componentCount == 0)
          // there's nothing in 'ere.
          return 0;
        else {
          Location last = (Location) componentList.get(componentCount - 1);
          return last.getMax();
        }
      } else {
        return length;
      }
    }

  public Symbol symbolAt(int pos) {
      Location l = locationOfPoint(pos);

      if (l == null) {
//          System.out.println("symbol at " + pos + " is n");
          return noninformativeSymbol;
      }

      SymbolList syms = getComponentSymbols(l);

      try {
//          System.out.println("symbol at " + pos + " is " + DNATools.getDNA().getTokenization("token").tokenizeSymbol(syms.symbolAt(pos)) 
//          + " location is " + l + " symbolList is " + syms);
      }
      catch (Exception be) { be.printStackTrace(); }
      return syms.symbolAt(pos);
  }

  public SymbolList subList(int start, int end) 
  {
    Location l = locationOfPoint(start);
    if (l != null && l.contains(end)) {
      SymbolList symbols = getComponentSymbols(l);
//      int tstart = start - l.getMin() + 1;
//      int tend = end - l.getMin() + 1;
      return symbols.subList(start, end);
    }

    // All is lost.  Fall back onto `view' subList from AbstractSymbolList

    return super.subList(start, end);
  }
/*
  public Symbol symbolAt(int pos) {
      // System.out.println(this + "  symbolAt(" + pos + ")");
    Location l = locationOfPoint(pos);
    if (l != null) {
      SymbolList syms = getComponentSymbols(l);
      return syms.symbolAt(pos - l.getMin() + 1);
    }

    return noninformativeSymbol;
  }


  public SymbolList subList(int start, int end) {
      // System.out.println(this + "  subList(" + start + ", " + end + ")");
    Location l = locationOfPoint(start);
    if (l != null && l.contains(end)) {
      SymbolList symbols = getComponentSymbols(l);
      int tstart = start - l.getMin() + 1;
      int tend = end - l.getMin() + 1;
      return symbols.subList(tstart, tend);
    }

    // All is lost.  Fall back onto `view' subList from AbstractSymbolList

    return super.subList(start, end);
  }
*/
}
