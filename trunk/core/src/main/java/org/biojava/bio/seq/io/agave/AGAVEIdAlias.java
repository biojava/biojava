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

public class AGAVEIdAlias {

  private String id, type ;



  AGAVEIdAlias(String id, String type) {

      this.id = id ;

      this.type = type;

  }

  public String getID(){

      return id ;

  }

  public String getType(){

      return type ;

  }

  public String toString(String indent, String indent_unit)

  {

     StringBuffer sb = new StringBuffer();

     sb.append(indent + "<id_alias id=\"" + id + "\"");

     if(  type != null )

        sb.append(" type=\"" + type + "\"");

     sb.append("/>");

     return sb.substring(0);

  }

  public String toString()

  {

     StringBuffer sb = new StringBuffer();

     sb.append("<id_alias id=\"" + id + "\"");

     if(  type != null )

        sb.append(" type=\"" + type + "\"");

     sb.append("/>");

     return sb.substring(0);

  }



}