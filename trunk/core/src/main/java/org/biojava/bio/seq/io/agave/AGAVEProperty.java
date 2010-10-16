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
 * @author Primary author unknown
 * @author Greg Cox
 */

package org.biojava.bio.seq.io.agave;



public class AGAVEProperty

{

   public static final String SCI_PROPERTY = "sci_property" ;

   public static final String RESULT_PROPERTY = "result_property" ;

   public static final String XREF_PROPERTY = "xref_property" ;

   public static final String QUALIFIER = "qualifier" ;



   private String prop ;

   private String data_type = "varchar";

   private String value ;

   private String property_class = "result_property" ;

   public AGAVEProperty(String property_class, String prop, String data_type, String value)

   {

       if( property_class != null && !( property_class.equals( SCI_PROPERTY) ||

                                        property_class.equals( RESULT_PROPERTY) || 

					property_class.equals( XREF_PROPERTY) ||

					property_class.equals( QUALIFIER)) )

           throw new RuntimeException("unkownn parameter for AGAVEProperty constructor :" + property_class) ;

       if( property_class != null )

           this.property_class = property_class  ;

       this.prop = prop ;

       if( data_type != null )

           this.data_type =   data_type  ;

       this.value = value ;

   }

   public String getPropClass()

   {

      return property_class ;

   }

   public String getPropType()

   {

      return prop ;

    }

    public String getDataType()

    {

       return data_type;

    }

    public String getValue()

    {

       return value ;

    }

    public String toString(String indent, String indent_unit)

    {

        StringBuffer sb = new StringBuffer() ;

        sb.append(indent + "<" + property_class ) ;

        if( property_class.equals(QUALIFIER)) {

           sb.append( " qualifier_type=\"" ) ;

        } else {

	    sb.append( " prop_type=\"" ) ;   // Fixed.  Must be a bug, since it was giving duplicated attributes THOMASD

	}

	sb.append(prop);

        sb.append("\" data_type=\"" + data_type + "\">" + value + "</" + property_class + ">" + "\n" ) ;



        return sb.substring(0);

    }

    public String toString()

    {

        StringBuffer sb = new StringBuffer() ;

        sb.append("<" + property_class ) ;

        if( property_class.equals(QUALIFIER)) {

           sb.append( " qualifier_type=\"" ) ;

        } else {

	    sb.append( " prop_type=\"" ) ;   // Fixed.  Must be a bug, since it was giving duplicated attributes THOMASD

	}

	sb.append(prop);

        sb.append("\" data_type=\"" + data_type + "\">" + value + "</" + property_class + ">" + "\n" ) ;



        return sb.substring(0);

    }

}

