// $Id: PhylogenyDataUtil.java,v 1.31 2009/11/01 04:18:53 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny.data;

import java.awt.Graphics;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;

import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public final class PhylogenyDataUtil {

    public static void appendClose( final Writer w, final String element_name ) throws IOException {
        w.write( "</" );
        w.write( element_name );
        w.write( ">" );
    }

    public static void appendElement( final Writer w, final String element_name, final String value )
            throws IOException {
        appendOpen( w, element_name );
        w.write( replaceIllegalXmlCharacters( value ) );
        appendClose( w, element_name );
    }

    public static void appendElement( final Writer w,
                                      final String element_name,
                                      final String value,
                                      final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        w.write( PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
        // Something like this replacement needs to be done in a more systematic manner.
        appendElement( w, element_name, value );
    }

    public static void appendElement( final Writer w,
                                      final String element_name,
                                      final String value,
                                      final String attribute_name,
                                      final String attribute_value ) throws IOException {
        appendOpen( w, element_name, attribute_name, attribute_value );
        w.write( replaceIllegalXmlCharacters( value ) );
        appendClose( w, element_name );
    }

    public static void appendElement( final Writer w,
                                      final String element_name,
                                      final String value,
                                      final String attribute_name,
                                      final String attribute_value,
                                      final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        w.write( PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
        appendOpen( w, element_name, attribute_name, attribute_value );
        w.write( replaceIllegalXmlCharacters( value ) );
        appendClose( w, element_name );
    }

    public static void appendElement( final Writer w,
                                      final String element_name,
                                      final String value,
                                      final String attribute1_name,
                                      final String attribute1_value,
                                      final String attribute2_name,
                                      final String attribute2_value,
                                      final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        w.write( PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
        appendOpen( w, element_name, attribute1_name, attribute1_value, attribute2_name, attribute2_value );
        w.write( replaceIllegalXmlCharacters( value ) );
        appendClose( w, element_name );
    }

    public static void appendElement( final Writer w,
                                      final String element_name,
                                      final String attribute1_name,
                                      final String attribute1_value,
                                      final String attribute2_name,
                                      final String attribute2_value,
                                      final String attribute3_name,
                                      final String attribute3_value,
                                      final String attribute4_name,
                                      final String attribute4_value,
                                      final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        appendOpen( w,
                    element_name,
                    attribute1_name,
                    attribute1_value,
                    attribute2_name,
                    attribute2_value,
                    attribute3_name,
                    attribute3_value,
                    attribute4_name,
                    attribute4_value );
        appendClose( w, element_name );
    }

    public static void appendElement( final Writer w,
                                      final String element_name,
                                      final String value,
                                      final String attribute1_name,
                                      final String attribute1_value,
                                      final String attribute2_name,
                                      final String attribute2_value,
                                      final String attribute3_name,
                                      final String attribute3_value,
                                      final String attribute4_name,
                                      final String attribute4_value,
                                      final String attribute5_name,
                                      final String attribute5_value,
                                      final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        w.write( PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
        appendOpen( w,
                    element_name,
                    attribute1_name,
                    attribute1_value,
                    attribute2_name,
                    attribute2_value,
                    attribute3_name,
                    attribute3_value,
                    attribute4_name,
                    attribute4_value,
                    attribute5_name,
                    attribute5_value );
        w.write( replaceIllegalXmlCharacters( value ) );
        appendClose( w, element_name );
    }

    public static void appendOpen( final Writer w, final String element_name ) throws IOException {
        w.write( "<" );
        w.write( element_name );
        w.write( ">" );
    }

    public static void appendOpen( final Writer w,
                                   final String element_name,
                                   final String attribute_name,
                                   final String attribute_value ) throws IOException {
        w.write( "<" );
        w.write( element_name );
        if ( !ForesterUtil.isEmpty( attribute_value ) ) {
            w.write( " " );
            w.write( attribute_name );
            w.write( "=\"" );
            w.write( attribute_value );
            w.write( "\"" );
        }
        w.write( ">" );
    }

    public static void appendOpen( final Writer w,
                                   final String element_name,
                                   final String attribute1_name,
                                   final String attribute1_value,
                                   final String attribute2_name,
                                   final String attribute2_value ) throws IOException {
        w.write( "<" );
        w.write( element_name );
        if ( !ForesterUtil.isEmpty( attribute1_value ) ) {
            w.write( " " );
            w.write( attribute1_name );
            w.write( "=\"" );
            w.write( attribute1_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute2_value ) ) {
            w.write( " " );
            w.write( attribute2_name );
            w.write( "=\"" );
            w.write( attribute2_value );
            w.write( "\"" );
        }
        w.write( ">" );
    }

    public static void appendOpen( final Writer w,
                                   final String element_name,
                                   final String attribute1_name,
                                   final String attribute1_value,
                                   final String attribute2_name,
                                   final String attribute2_value,
                                   final String attribute3_name,
                                   final String attribute3_value ) throws IOException {
        w.write( "<" );
        w.write( element_name );
        if ( !ForesterUtil.isEmpty( attribute1_value ) ) {
            w.write( " " );
            w.write( attribute1_name );
            w.write( "=\"" );
            w.write( attribute1_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute2_value ) ) {
            w.write( " " );
            w.write( attribute2_name );
            w.write( "=\"" );
            w.write( attribute2_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute2_value ) ) {
            w.write( " " );
            w.write( attribute3_name );
            w.write( "=\"" );
            w.write( attribute3_value );
            w.write( "\"" );
        }
        w.write( ">" );
    }

    public static void appendOpen( final Writer w,
                                   final String element_name,
                                   final String attribute1_name,
                                   final String attribute1_value,
                                   final String attribute2_name,
                                   final String attribute2_value,
                                   final String attribute3_name,
                                   final String attribute3_value,
                                   final String attribute4_name,
                                   final String attribute4_value ) throws IOException {
        w.write( "<" );
        w.write( element_name );
        if ( !ForesterUtil.isEmpty( attribute1_value ) ) {
            w.write( " " );
            w.write( attribute1_name );
            w.write( "=\"" );
            w.write( attribute1_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute2_value ) ) {
            w.write( " " );
            w.write( attribute2_name );
            w.write( "=\"" );
            w.write( attribute2_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute3_value ) ) {
            w.write( " " );
            w.write( attribute3_name );
            w.write( "=\"" );
            w.write( attribute3_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute4_value ) ) {
            w.write( " " );
            w.write( attribute4_name );
            w.write( "=\"" );
            w.write( attribute4_value );
            w.write( "\"" );
        }
        w.write( ">" );
    }

    public static void appendOpen( final Writer w,
                                   final String element_name,
                                   final String attribute1_name,
                                   final String attribute1_value,
                                   final String attribute2_name,
                                   final String attribute2_value,
                                   final String attribute3_name,
                                   final String attribute3_value,
                                   final String attribute4_name,
                                   final String attribute4_value,
                                   final String attribute5_name,
                                   final String attribute5_value ) throws IOException {
        w.write( "<" );
        w.write( element_name );
        if ( !ForesterUtil.isEmpty( attribute1_value ) ) {
            w.write( " " );
            w.write( attribute1_name );
            w.write( "=\"" );
            w.write( attribute1_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute2_value ) ) {
            w.write( " " );
            w.write( attribute2_name );
            w.write( "=\"" );
            w.write( attribute2_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute3_value ) ) {
            w.write( " " );
            w.write( attribute3_name );
            w.write( "=\"" );
            w.write( attribute3_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute4_value ) ) {
            w.write( " " );
            w.write( attribute4_name );
            w.write( "=\"" );
            w.write( attribute4_value );
            w.write( "\"" );
        }
        if ( !ForesterUtil.isEmpty( attribute5_value ) ) {
            w.write( " " );
            w.write( attribute5_name );
            w.write( "=\"" );
            w.write( attribute5_value );
            w.write( "\"" );
        }
        w.write( ">" );
    }

    /**
     * Creates a deep copy of ArrayList of PhylogenyData objects.
     * 
     * @param list
     *            an ArrayList of PhylogenyData objects
     * @return a deep copy of ArrayList list
     */
    public static ArrayList<PhylogenyData> copy( final ArrayList<PhylogenyData> list ) {
        final ArrayList<PhylogenyData> l = new ArrayList<PhylogenyData>( list.size() );
        for( int i = 0; i < list.size(); ++i ) {
            l.add( ( list.get( i ) ).copy() );
        }
        return l;
    }

    public static void drawLine( final double x1, final double y1, final double x2, final double y2, final Graphics g ) {
        g.drawLine( org.forester.util.ForesterUtil.roundToInt( x1 ),
                    org.forester.util.ForesterUtil.roundToInt( y1 ),
                    org.forester.util.ForesterUtil.roundToInt( x2 ),
                    org.forester.util.ForesterUtil.roundToInt( y2 ) );
    }

    public static void drawString( final String str, final double x, final double y, final Graphics g ) {
        g.drawString( str, org.forester.util.ForesterUtil.roundToInt( x ), org.forester.util.ForesterUtil
                .roundToInt( y ) );
    }

    public static String replaceIllegalXmlCharacters( final String value ) {
        String v = value.replaceAll( "&", "&amp;" );
        v = v.replaceAll( "<", "&lt;" );
        v = v.replaceAll( ">", "&gt;" );
        v = v.replaceAll( "'", "&apos;" );
        v = v.replaceAll( "\"", "&quot;" );
        return v;
    }
}
