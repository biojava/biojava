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

/**
 * <p>AbstractULAlignment is an abstract base class for alignments
 * where the constituent sequences have unequal lengths.</p>
 *
 * @author David Waring
 * @author Mark Schreiber (docs and some minor changes)
 */


package org.biojava.bio.alignment;

import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.AbstractSymbolList;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;


public abstract class AbstractULAlignment extends AbstractSymbolList implements UnequalLengthAlignment{

    protected Alphabet alphabet;
    /** this will return the ambiguity symbol associated with all symbols in that column
    */

    public Symbol symbolAt(int index) {
        try {
          return alphabet.getSymbol(new ColAsList(index));
        } catch (IllegalSymbolException ire) {
          throw new BioError(
            "Somehow my crossproduct alphabet is incompatible with column " + index, ire
          );
        }
    }

    public List labelsAt(int column){
        return labelsInRange(new RangeLocation(column,column));
    }

    public List labelsInRange(Location loc){
        List labels = getLabels();
        Location seqLoc;
        Object label;
        List labelList = new Vector();

        for (Iterator labelIterator = labels.iterator();labelIterator.hasNext();){
            label = labelIterator.next();
            seqLoc = locInAlignment(label);
            if (LocationTools.overlaps(loc,seqLoc)){
                labelList.add(label);
            }
        }
        return labelList;
    }

    public Iterator symbolListIterator() {
      return new Alignment.SymbolListIterator(this);
    }

    protected void debug(String s){
        System.out.println(s);
    }
    /**
    * leftMost and rightMost return labels. If there are more than one
    * that start at the same location it returns the longest, if they are
    * the same length it returns the first one it found;
    */
    public Object leftMost(){
        List labels = getLabels();
        Object leftMost = null;
        Object label;
        Location leftLoc= null;
        Location loc = null;

        for (Iterator labelIterator = labels.iterator();labelIterator.hasNext();){
            label = labelIterator.next();
            loc = locInAlignment(label);
            if (leftMost == null){
                leftMost = label;
                leftLoc= loc;
                //loc = loc;
            }
            if (loc.getMin() < leftLoc.getMin()){
                leftMost = label;
                leftLoc = loc;
            }else if (loc.getMin() == leftLoc.getMin()){
                if ((loc.getMax() - loc.getMin()) > (leftLoc.getMax() - leftLoc.getMin())){
                    leftMost = label;
                    leftLoc = loc;
                }
            }
        }
        return leftMost();
    }


    public Object rightMost(){
        List labels = getLabels();
        Object rightMost = null;
        Object label;
        Location rightLoc= null;
        Location loc= null;

        for (Iterator labelIterator = labels.iterator();labelIterator.hasNext();){
            label = labelIterator.next();
            loc = locInAlignment(label);
            if (rightMost == null){
                rightMost = label;
                rightLoc= loc;
                //loc= loc;
            }
            if (loc.getMin() > rightLoc.getMin()){
                rightMost = label;
                rightLoc = loc;
            }else if (loc.getMin() == rightLoc.getMin()){
                if ((loc.getMax() - loc.getMin()) > (rightLoc.getMax() - rightLoc.getMin())){
                    rightMost = label;
                    rightLoc = loc;
                }
            }
        }
        return rightMost();
    }

    /**
     * Retrieves a subalignment specified by the location.
     * <p><b>WARNING:</b> It is assumed that the location is contiguous. If
     * the location is non-contiguous it may be preferable to use a block iterator
     * to retrieve each sub location independently.
     * @see #subAlignment(Set labels, int min, int max)
     */
    public Alignment subAlignment(Set labels, Location loc) throws IndexOutOfBoundsException{
        return new SubULAlignment(labels,loc);
    }

    /**
      * Retreives a subAlignment
      * @param labels the labels of the <code>SymbolLists</code> to be in the Alignment
      * @param min the left most coordinate
      * @param max the right most coordinate
      * @return an Alignment
      * @throws NoSuchElementException if one of the values in <code>labels</code> is not in the parent alignment
      */
    public Alignment subAlignment(Set labels, int min, int max) throws NoSuchElementException{
       return subAlignment(labels, LocationTools.makeLocation(min, max));
    }

    public SortedSet orderedLables(Comparator comp){
        TreeSet orderedSet = new TreeSet(comp);
        orderedSet.addAll(getLabels());
        return orderedSet;
    }







   //////////////////////////////////////////////
   //////////////////////////////////////////////
   // INNER CLASSES
   //////////////////////////////////////////////
   //////////////////////////////////////////////

    private final class ColAsList extends AbstractList implements Serializable {
        private final int col;
        private List labels;

        public ColAsList(int col) {
          this.col = col;
          labels = getLabels();
        }

        protected ColAsList() {
          this.col = 0;
        }

        public Object get(int indx) {
          return symbolAt(labels.get(indx), col);
        }

        public int size() {
          return labels.size();
        }
    }

    /**
    * Orders by location left to right. If they both start at the same location (o1.getMin() == o2.getMin())
    * it the larger of the two considered to the Left
    */

    public class LeftRightLocationComparator<T> implements Comparator<T> {
       public int compare(Object o1, Object o2){
            int ret = 1;
            Location loc1,loc2;
            loc1 = locInAlignment(o1);
            loc2 = locInAlignment(o2);
            if (loc1.getMin() > loc2.getMin()){
                ret = 1;
            }else if (loc1.getMin() < loc2.getMin()){
                ret = -1;
            }else if (loc1.getMin() == loc2.getMin()){
                int s1 = (loc1.getMax() - loc1.getMin()) + 1;
                int s2 = (loc2.getMax() - loc2.getMin()) + 1;
                if (s1 == s2){
                    ret = 1;
                }else{
                    ret =  s2 - s1;
                }
            }
//            debug (" result " + ret);
            return ret;
        }
    }


    ///////////////
    // This inner class should take care of all subAlignments
    ///////////////

    public class SubULAlignment extends AbstractSymbolList implements UnequalLengthAlignment{
        private int start,end;
        private List subLabels;  // will be left null if constructed with null Set of labels
                                 // this allows labels added to underlying Alignment to be seen
                                 // unless it was constructed with a Set of labels  NOT

        protected SubULAlignment (Set labels, Location loc) throws IndexOutOfBoundsException{
            this.start = loc.getMin();
            this.end = loc.getMax();
            if (start < 1 || end > AbstractULAlignment.this.length()){
                throw new IndexOutOfBoundsException();
            }
            if (labels != null){
                subLabels = new ArrayList();
                subLabels.addAll(labels);
            }else{
                subLabels = AbstractULAlignment.this.labelsInRange(new RangeLocation(start,end));
            }

        }

        /**
        * realPosition is the position in the underlying Alignment corresponding to
        * a position in the subAlignment
        */

        private int realPosition(int pos){
            return pos + start - 1;
        }
        /////////////////////////
        // methods from Interface UnequalLengthAlignment
        ////////////////////////
        public int length(){
            return end - start + 1;
        }

        /**
        * The location of an individual SymbolList relative to overall Alignment
        */
        public Location locInAlignment(Object label){
            Location origLoc = AbstractULAlignment.this.locInAlignment(label);
            int min = origLoc.getMin() - start + 1;
            int max = origLoc.getMax() - start + 1;
            return new RangeLocation(min,max);
        }




        public Alignment subAlignment(Set labels, Location loc) throws NoSuchElementException{
            int min = realPosition(loc.getMin());
            int max = realPosition(loc.getMax());
            // if the Subalignment has a limited subLabel we want to keep labels of the sub sub aligment limited to that list too

            if (labels == null){
                labels = new TreeSet(subLabels);
            }else{
                Vector l = new Vector(labels);
                labels = new TreeSet(listIntersection(l,subLabels));
            }
//            debug (" labels now " + labels);
            return new SubULAlignment(labels,new RangeLocation(min,max));
        }
        protected List listIntersection(List s1,List s2){
           List common = new Vector(s1);

            for (Iterator i = s1.iterator();i.hasNext();){
                 Object label = i.next();
                 if ( ! s2.contains(label)){
                    common.remove(label);
                 }
            }
            return common;
        }
        public List labelsAt(int column) throws IndexOutOfBoundsException {
            return labelsInRange(new RangeLocation(column,column));
        }

        public List labelsInRange(Location loc) throws IndexOutOfBoundsException {
//            debug ("looking for labelsInRange " + loc.getMin() + "-" + loc.getMax());
            int min = realPosition(loc.getMin());
            int max = realPosition(loc.getMax());
            if (min < start || max > end){
                throw new IndexOutOfBoundsException();
            }
            return listIntersection(subLabels,AbstractULAlignment.this.labelsInRange(new RangeLocation(min,max)));
//            List aLabels = AbstractULAlignment.this.labelsInRange(new RangeLocation(min,max));
//            List sLabels = AbstractULAlignment.this.labelsInRange(new RangeLocation(min,max));
//            for (Iterator i = aLabels.iterator();i.hasNext();){
//                 Object label = i.next();
//                 if ( ! subLabels.contains(label)){
//                    sLabels.remove(label);
//                 }
//            }
//            return sLabels;

        }



        /////////////////////////
        // methods from Interface Alignment
        ////////////////////////
        public List getLabels(){
            return subLabels;
        }

        public Symbol symbolAt(Object label, int column) throws NoSuchElementException{
           return AbstractULAlignment.this.symbolAt(label,realPosition(column));
        }

        public Symbol symbolAt(int column) throws NoSuchElementException{
           return AbstractULAlignment.this.symbolAt(realPosition(column));
        }

        public SymbolList symbolListForLabel(Object label) throws NoSuchElementException{
           return AbstractULAlignment.this.symbolListForLabel(label);
        }

        public Alphabet getAlphabet(){
           return AbstractULAlignment.this.getAlphabet();
        }

        public Iterator symbolListIterator() {
          return new Alignment.SymbolListIterator(this);
        }
    }
}
