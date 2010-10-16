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
 * xref
 * @author Hanning Ni    Doubletwist Inc
 * @author Greg Cox
*/
public class AGAVEXref {
    private AGAVEDbId id ;
    private List xrefs ;
    private String rel ;
    public void setDbId( AGAVEDbId id )
    {
        this.id = id ;
    }
    public void addProp(AGAVEProperty prop)
    {
        if( xrefs  == null )
            xrefs = new ArrayList(1) ;
        xrefs.add( prop );
    }
    public void setRel(String rel)
    {
        this.rel = rel ;
    }
    public String getRel()
    {
        return rel  ;
    }
    public AGAVEDbId getDbId()
    {
        return id ;
    }
    public Iterator getXrefProps()
    {
        return xrefs.iterator() ;
    }

    /** return the agave xml representation of this instance **/
    public String toString(String indent, String indent_unit)
    {
       StringBuffer tmp = new StringBuffer();
       tmp.append(indent +  "<xref");
       if( rel != rel )
           tmp.append(" relationship=\"" + rel + "\"") ;
       tmp.append( ">" + "\n") ;
       tmp.append( id.toString(indent + indent_unit, indent_unit) ) ;
       Iterator it = xrefs.iterator() ;
       while( it.hasNext() )
       {
           tmp.append( ((AGAVEProperty)it.next()).toString(indent + indent_unit, indent_unit) ) ;
       }
       tmp.append(indent + "</xref>" + "\n" );
       return tmp.substring(0) ;
    }
    /** return the agave xml representation of this instance **/
    public String toString()
    {
       StringBuffer tmp = new StringBuffer();
       tmp.append( "<xref");
       if( rel != rel )
           tmp.append(" relationship=\"" + rel + "\"") ;
       tmp.append( ">" + "\n") ;
       tmp.append( id ) ;
       Iterator it = xrefs.iterator() ;
       while( it.hasNext() )
       {
           tmp.append( it.next().toString() ) ;
       }
       tmp.append("</xref>" + "\n" );
       return tmp.toString() ;
    }
}
