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

package org.biojava.bio.seq.io.agave;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;



/**

 * 

 * @author Hanning Ni    Doubletwist Inc
  * @author Greg Cox

 */

public class AGAVERelatedAnnot {

  private List element_ids ;

  private List props ;

  private String rel ;

  private String score ;



  /** construct.. **/

  public AGAVERelatedAnnot() {

      element_ids = new ArrayList(1) ;

  }





  public void setRel(String rel)

  {

      this.rel = rel  ;

  }

  public void setScore(String score)

  {

      this.score = score ;

  }

  /**  @param id is content of the tag of <element_id> **/

  public void addElementId(String id)

  {

      element_ids .add( id );

  }

  public void addProp(AGAVEProperty prop)

  {

      if( props == null ) //lazy realization

          props = new ArrayList(1) ;

      props.add(prop);

  }

  public String getScore()

  {

      return score ;

  }

  public String getRel()

  {

      return rel ;

  }

  public String[] getElementIds()

  {

      String[]  tmp = new String[ element_ids.size() ] ;

      return (String[]) element_ids.toArray( tmp ) ;

  }

  public AGAVEProperty[] getProps()

  {

      AGAVEProperty[] tmp = new AGAVEProperty[ props.size() ] ;

      return (AGAVEProperty[]) props.toArray( tmp ) ;

  }

  public String toString(String indent, String indent_unit)

  {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append(indent +  "<related_annot rel=\"" + rel + "\"" ) ;

      if( score != null )

         tmp.append(" score=\"" + score + "\"") ;

      tmp.append(">" + "\n") ;

      Iterator it = element_ids.iterator() ;

      while( it.hasNext() )

      {

          tmp.append( indent + indent_unit + "<element_id id=\"" + (String)it.next() + "\"/>" + "\n" ) ;

      }

      it = props.iterator() ;

      while( it.hasNext() )

      {

          tmp.append( ((AGAVEProperty) it.next()).toString(indent + indent_unit, indent_unit) ) ;

      }

      tmp.append(indent + "</related_annot>" + "\n") ;

      return tmp.substring(0) ;

  }

   public String toString()

  {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append( "<related_annot rel=\"" + rel + "\"" ) ;

      if( score != null )

         tmp.append(" score=\"" + score + "\"") ;

      tmp.append(">" + "\n") ;

      Iterator it = element_ids.iterator() ;

      while( it.hasNext() )

      {

          tmp.append("<element_id id=\"" + (String)it.next() + "\"/>" + "\n" ) ;

      }

      it = props.iterator() ;

      while( it.hasNext() )

      {

          tmp.append( (AGAVEProperty) it.next() ) ;

      }

      tmp.append("</related_annot>" + "\n") ;

      return tmp.substring(0) ;

  }

}