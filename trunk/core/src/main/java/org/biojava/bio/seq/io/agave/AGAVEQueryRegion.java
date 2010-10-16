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

 * 

 * @author Hanning Ni    Doubletwist Inc
  * @author Greg Cox

 */

public class AGAVEQueryRegion

{

    private int start ;

    private int end ;

    private AGAVEDbId db_id ;





    public void setStart(int start)

    {

        this.start = start ;

    }



    public void setEnd(int end)

    {

        this.end = end ;

    }



    public void setDbId(AGAVEDbId id)

    {

        this.db_id = id ;

    }



    public int getStart()

    {

        return start ;

    }

    public int getEnd()

    {

        return end ;

    }

    public AGAVEDbId getDbId()

    {

        return db_id ;

    }



    public String toString(String indent, String indent_unit)

    {

        StringBuffer sb = new StringBuffer();

        sb.append( indent +"<query_region start=\"" + start + "\" end=\"" + end +"\">" + "\n");

        if( db_id != null )

          sb.append( db_id.toString(indent + indent_unit, indent_unit) ) ;

        sb.append( indent +  "</query_region>" + "\n" );

        return sb.substring(0) ;

    }



    public String toString()

    {

        return "<query_region start=\"" + start + "\" end=\"" + end +"\">" + "\n"

               + db_id == null ? "" : db_id

               + "</query_region>" + "\n" ;

    }

}

