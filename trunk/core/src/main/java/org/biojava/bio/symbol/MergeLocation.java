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
 */

package org.biojava.bio.symbol;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.biojava.bio.BioException;

/**
 * Produced by LocationTools as a result of union operations.
 * It's a <code>RangeLocation</code> and can be used
 * as such but it also retains knowledge of which individual components made
 * it up. None of the methods of RangeLocation are overridden only new methods
 * have been added to get the subcomponents.</p>
 *
 *  <p>For example a union operation between the following locations
 * [1,20],[27,45],[30-70] will produce a <code>CompoundLocation</code>
 * like this: {[1,20],[27,70]}, the last two locations have been merged into a
 * <code>MergeLocation</code> which contains the two subcomponents.</p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: AgResearch</p>
 * @author Mark Schreiber
 * @author Francois Pepin
 * @version 1.0
 */

public class MergeLocation extends RangeLocation {
  List componentLocations;

  /**
   * Private constructor, use static method mergeLocations to make an
   * instance of this class
   * @param min the minimum point spanned by the location
   * @param max the maximum point spanned by the location
   * @param componentLocations the locations that make up this.
   */
  private MergeLocation(int min, int max, List componentLocations){
    super(min, max);
    this.componentLocations = componentLocations;
  }

  /**
   * Gets the component locations that make up this one
   *
   * @param recurse if true the method lists the component locations of all
   * nested <code>MergedLocation</code>s.
   *
   * @return a <code>List</code> of <code>Location</code>s.
   */
  public List getComponentList(boolean recurse){
    if(! recurse)
      return componentLocations;

    else{
      List l = new ArrayList();
      for(Iterator i = componentLocations.iterator(); i.hasNext();){
        Object o = i.next();

        //what if its decorated?
        if(o instanceof AbstractLocationDecorator){
          o = ((AbstractLocationDecorator)o).getWrapped();
        }

        if (o instanceof MergeLocation) {
          List ll = ((MergeLocation)o).getComponentList(true);
          l.addAll(ll);
        }
        else {
          l.add(o);
        }

      }

      return l;
    }
  }

  /**
   * @return A <code>ListIterator</code> over the component locations. The iterator
   * does not recurse nested MergeLocations.
   */
  public ListIterator componentLocationIterator(){
    return componentLocations.listIterator();
  }

  /**
   * Static Factory method for getting an instance of a <code>MergeLocation</code>
   * @param componentLocations the <code>Locations to Merge</code>
   * @return the merged location
   * @throws BioException if the list contains objects that are not <code>Locations</code>
   *  or if the locations don't represent a contiguous block. Use <code>
   *  LocationTools.union()</code> if you want to merge discontinuous blocks.
   */
  public static MergeLocation mergeLocations(List componentLocations)
      throws BioException{

    Collections.sort(componentLocations, Location.naturalOrder);

    int lastMax = -1;
    int lastMin = -1;

    int min = Integer.MAX_VALUE;
    int max = 0;

    for (Iterator i = componentLocations.iterator(); i.hasNext(); ) {
      Object item = i.next();

      if(!( item instanceof Location)){
        throw new BioException(
            "All members of the component locations list must be Location objects");
      }

      Location loc = (Location)item;

      //make sure blocks are contiguous
      if(lastMin != -1 && lastMax != -1){
        if(!(loc.getMin()+1 <= lastMax)){
          //blocks aren't contiguous
          throw new BioException(
            "All members of the component locations list must be contiguous");
        }
      }

      //determine what the min and max should be for the MergeLocation
      if (loc.getMin() < min) {
        min = loc.getMin();
      }
      if (loc.getMax() > max) {
        max = loc.getMax();
      }

      lastMin = loc.getMin(); lastMax = loc.getMax();
    }

    return new MergeLocation(min, max, componentLocations);
  }

  public static MergeLocation mergeLocations(Location locA,
      Location locB) throws BioException{

    int min = Math.min(locA.getMin(), locB.getMin());
    int max = Math.max(locA.getMax(), locB.getMax());

    List l = new ArrayList(2);
    l.add(locA);
    l.add(locB);

    return new MergeLocation(min, max, l);
  }
}
