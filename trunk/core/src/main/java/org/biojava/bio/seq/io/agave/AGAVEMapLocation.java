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

public class AGAVEMapLocation {

  private AGAVEMapPosition first ;

  private AGAVEMapPosition second ;

  private String map_type ;

  private String source ;

  private String units ;

  private String chromosome ;

  private String subseq_start ;

  private String orientation ;

  public AGAVEMapLocation() {

  }

  public void  addPosition(AGAVEMapPosition pos)

  {

      if( first == null )

          first = pos ;

      else

          second = pos ;

  }

  public AGAVEMapPosition getPosition1()

  {

      return first ;

  }



  public AGAVEMapPosition getPosition2()

  {

      return second;

  }



  public void setMapType(String value)

  {

      map_type = value ;

  }

  public String getMapType()

  {

     return map_type ;

  }

  public void setUnits(String value)

  {

      units = value ;

  }

  public String getUnits()

  {

      return units ;

  }

  public void setSource(String value)

  {

      source= value ;

  }

  public String getSource()

  {

      return source ;

  }

  public void setChromosome(String value)

  {

      chromosome = value ;

  }

  public String getChromosome()

  {

      return chromosome ;

  }

  public void setSubSeqStart(String value)

  {

      subseq_start = value ;

  }

  public String getSubSeqStart()

  {

      return subseq_start ;

  }

    public void setOrientation(String value)

  {

      orientation = value ;

  }

  public String getOrientation()

  {

      return orientation ;

  }

  public String toString(String indent, String indent_unit)

  {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append(indent + "<map_location ") ;

      if( map_type != null )

         tmp.append( " map_type=\"" + map_type + "\" " );

      if( source != null )

         tmp.append( " source=\"" + source + "\" " );

      if( units != null )

         tmp.append( " units=\"" + units + "\" " );

      if( orientation != null )

         tmp.append( " orientation=\"" + orientation + "\" " );

      if( chromosome != null )

         tmp.append( " chromosome=\"" + chromosome + "\" " );

      if( subseq_start != null )

         tmp.append( " subseq_start=\"" + subseq_start + "\" " );

      tmp.append( ">" + "\n") ;

      if( first != null )

         tmp.append ( first.toString(indent + indent_unit, indent_unit) ) ;

      if( second != null )

         tmp.append( second.toString(indent + indent_unit, indent_unit) ) ;

      tmp.append(indent + "</map_location>"+ "\n" ) ;

      return tmp.substring(0) ;

  }

  public String toString()

  {

      StringBuffer tmp = new StringBuffer() ;

      tmp.append("<map_location ") ;

      if( map_type != null )

         tmp.append( " map_type=\"" + map_type + "\" " );

      if( source != null )

         tmp.append( " source=\"" + source + "\" " );

      if( units != null )

         tmp.append( " units=\"" + units + "\" " );

      if( orientation != null )

         tmp.append( " orientation=\"" + orientation + "\" " );

      if( chromosome != null )

         tmp.append( " chromosome=\"" + chromosome + "\" " );

      if( subseq_start != null )

         tmp.append( " subseq_start=\"" + subseq_start + "\" " );

      tmp.append( ">" + "\n") ;

      if( first != null )

         tmp.append ( first ) ;

      if( second != null )

         tmp.append( second ) ;

      tmp.append( "</map_location>"+ "\n" ) ;

      return tmp.substring(0) ;

  }



}