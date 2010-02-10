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


/**

 * @author Hanning Ni    Doubletwist Inc
  * @author Greg Cox

 */

public class AGAVEDbId {

   private String id ;

   private String version;

   private String db_code ;

   public String getId(){

      return id ;

   }

   public String getVersion()

   {

      return version ;

   }

   public String getDbCode()

   {

      return db_code ;

   }

   public void setId(String id )

   {

      this.id = id ;

   }

   public void setVersion(String version)

   {

      this.version = version ;

   }

   public void setDbCode(String code)

   {

      this.db_code = code ;

   }

   public String toString(String indent, String indent_unit)

   {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append(indent + "<db_id id=\"" + id +   "\" db_code=\"" + db_code +"\" " ) ;

      if( version != null )

          tmp.append(  " version=\"" + version + "\"") ;

      tmp.append("/>" + "\n");

      return tmp.substring(0) ;



   }

   public String toString()

   {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append("<db_id id=\"" + id +   "\" db_code=\"" + db_code +"\" " ) ;

      if( version != null )

          tmp.append(  " version=\"" + version + "\"") ;

      tmp.append("/>" + "\n");

      return tmp.substring(0) ;



   }

}

