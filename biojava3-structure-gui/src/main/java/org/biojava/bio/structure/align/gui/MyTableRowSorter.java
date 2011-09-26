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
 * Created on Mar 30, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;

import java.util.Comparator;

import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

public class MyTableRowSorter extends TableRowSorter<TableModel>
{

   public MyTableRowSorter(TableModel tm)
   {
      super(tm);
   }
   
   @Override
   public Comparator<?> getComparator(int column)
   {
     
      return new MyComparator(column);
      
   }
  
}

 class MyComparator implements Comparator<String> {
  
    int column;
    public MyComparator(int column){
       this.column = column;
    }
    
  public int compare(String s1, String s2){
	  if (( column >= 2 && column <= 4)|| (column==9)){
		  return compare(Float.parseFloat(s1), Float.parseFloat(s2));
	  } else if  (  column > 4 && column < 10) {
        return compare(Integer.parseInt(s1), Integer.parseInt(s2));
     } else
        return s1.compareTo(s2);
  }
    
  public int compare(Float f1, Float f2){
     return f1.compareTo(f2);
  }
   public int compare(Integer o1, Integer o2)
   {
     return o1.compareTo(o2);
   }
}
