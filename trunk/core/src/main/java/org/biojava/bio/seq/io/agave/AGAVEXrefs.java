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

 * xrefs

 *

 * @author Hanning Ni    Doubletwist Inc
  * @author Greg Cox

*/

public class  AGAVEXrefs{

     private List db_ids ;

     private List xrefs ;

     /** add @param id **/

     public void addDbId(AGAVEDbId id)

     {

         if( db_ids == null )

             db_ids =  new ArrayList(1) ;

         db_ids .add( id ) ;

     }

     /** add @param xref **/

     public void addXref(AGAVEXref xref)

     {

         if( xrefs == null )

             xrefs = new ArrayList(1) ;

         xrefs.add( xref ) ;

     }

     /** return a set of DbId **/

     public Iterator getDbIds()

     {

         return db_ids.iterator() ;

     }

     /** return a set of AGAVEXref **/

     public Iterator getXrefs()

     {

         return xrefs.iterator() ;

     }

     /**

      * @param indent  the leading space

      *  @param indent_unit the unit of indenting for xml format

      *

      **/

     public String toString(String indent, String indent_unit)

     {

       StringBuffer tmp = new StringBuffer();

       tmp.append(indent +  "<xrefs>" + "\n" );

       Iterator it = db_ids.iterator() ;

       while( it.hasNext() )

       {

           tmp.append( ((AGAVEDbId)it.next()).toString(indent + indent_unit, indent_unit) ) ;

       }

       it = xrefs.iterator() ;

       while( it.hasNext() )

       {

           tmp.append( ((AGAVEXref)it.next()).toString(indent + indent_unit, indent_unit) ) ;

       }

       tmp.append(indent + "</xrefs>" + "\n" );

       return tmp.substring(0) ;

     }

     /** the agave xml representation of xrefs **/

     public String toString()

     {

       StringBuffer tmp = new StringBuffer();

       tmp.append( "<xrefs>" + "\n" );

       Iterator it = db_ids.iterator() ;

       while( it.hasNext() )

       {

           tmp.append(it.next().toString() ) ;

       }

       it = xrefs.iterator() ;

       while( it.hasNext() )

       {

           tmp.append(it.next().toString() ) ;

       }

       tmp.append("</xrefs>" + "\n" );

       return tmp.substring(0) ;

     }

}

