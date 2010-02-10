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

public class AGAVEMapPosition {

  private String pos  ;

  private List db_ids ;



  public AGAVEMapPosition() {

  }

  public void  setPos(String pos)

  {

      this.pos = pos ;

  }

  public String getPos()

  {

      return pos ;

  }

  public void addDbId(AGAVEDbId id)

  {

      if( db_ids == null )

          db_ids = new ArrayList(1) ;

      db_ids.add( id );

  }

  public Iterator getDbIds()

  {

      return db_ids.iterator() ;

  }

  public String toString(String indent, String indent_unit)

  {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append(indent + "<map_position pos=\"" + pos + "\">" + "\n" ) ;

      Iterator it  = db_ids.iterator() ;

      while( it.hasNext() )

      {

          tmp.append( ((AGAVEDbId) it.next()).toString(indent + indent_unit, indent_unit) ) ;

      }

      tmp.append(indent + "</map_position>") ;

      return tmp.substring(0) ;

  }

  public String toString()

  {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append("<map_position pos=\"" + pos + "\">" + "\n" ) ;

      Iterator it  = db_ids.iterator() ;

      while( it.hasNext() )

      {

          tmp.append( (AGAVEDbId) it.next() ) ;

      }

      tmp.append("</map_position>") ;

      return tmp.substring(0) ;

  }

}