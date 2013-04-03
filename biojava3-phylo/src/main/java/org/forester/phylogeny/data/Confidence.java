// $Id: Confidence.java,v 1.22 2009/12/17 02:28:13 cmzmasek Exp $
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

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.util.ForesterUtil;

public class Confidence implements PhylogenyData, Comparable<Confidence> {

    public final static double CONFIDENCE_DEFAULT_VALUE = -9999.0;
    private double             _value;
    private String             _type;

    public Confidence() {
        init();
    }

    public Confidence( final double value, final String type ) {
        setValue( value );
        setType( type );
    }

    public StringBuffer asSimpleText() {
        return new StringBuffer().append( ForesterUtil.FORMATTER_6.format( getValue() ) );
    }

    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getType() ) ) {
            sb.append( "[" );
            sb.append( getType() );
            sb.append( "] " );
        }
        sb.append( ForesterUtil.FORMATTER_6.format( getValue() ) );
        return sb;
    }

    @Override
    public int compareTo( final Confidence confidence ) {
        if ( this == confidence ) {
            return 0;
        }
        return getType().compareToIgnoreCase( confidence.getType() );
    }

    public PhylogenyData copy() {
        return new Confidence( getValue(), new String( getType() ) );
    }

    public String getType() {
        return _type;
    }

    public double getValue() {
        return _value;
    }

    public void init() {
        setValue( CONFIDENCE_DEFAULT_VALUE );
        setType( "" );
    }

    public boolean isEqual( final PhylogenyData confidence ) {
        if ( confidence == null ) {
            return false;
        }
        if ( !( confidence instanceof Confidence ) ) {
            return false;
        }
        final Confidence s = ( Confidence ) confidence;
        if ( s.getValue() != getValue() ) {
            return false;
        }
        if ( !s.getType().equals( getType() ) ) {
            return false;
        }
        return true;
    }

    public void setType( final String type ) {
        _type = type;
    }

    public void setValue( final double value ) {
        _value = value;
    }

    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( NHXtags.SUPPORT );
        sb.append( getValue() );
        return sb;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( getValue() == CONFIDENCE_DEFAULT_VALUE ) {
            return;
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendElement( writer,
                                         PhyloXmlMapping.CONFIDENCE,
                                         String.valueOf( ForesterUtil
                                                 .round( getValue(),
                                                         PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ),
                                         PhyloXmlMapping.CONFIDENCE_TYPE_ATTR,
                                         ForesterUtil.isEmpty( getType() ) ? "unknown" : getType() );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
