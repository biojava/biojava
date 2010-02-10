// $Id: Reference.java,v 1.10 2009/12/12 00:14:39 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.util.ForesterUtil;

public class Reference implements PhylogenyData {

    String _value;
    String _doi;

    public Reference( final String value ) {
        _value = value;
        _doi = "";
    }

    public Reference( final String value, final String doi ) {
        _value = value;
        _doi = doi;
    }

    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() );
    }

    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getDoi() ) ) {
            sb.append( "[doi:" );
            sb.append( getDoi() );
            sb.append( "] " );
        }
        sb.append( getValue() );
        return sb;
    }

    public PhylogenyData copy() {
        return new Reference( new String( getValue() ), new String( getDoi() ) );
    }

    public String getDoi() {
        return _doi;
    }

    public String getValue() {
        return _value;
    }

    public boolean isEqual( final PhylogenyData data ) {
        if ( ( data == null ) || ( getValue() == null ) ) {
            return false;
        }
        return ( ( Reference ) data ).getValue().equals( getValue() )
                && ( ( Reference ) data ).getDoi().equals( getDoi() );
    }

    public void setDoi( final String doi ) {
        if ( !ForesterUtil.isEmpty( doi ) && !PhyloXmlUtil.LIT_REF_DOI_PATTERN.matcher( doi ).matches() ) {
            throw new PhyloXmlDataFormatException( "illegal doi: [" + doi + "]" );
        }
        _doi = doi;
    }

    public void setValue( final String value ) {
        _value = value;
    }

    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.REFERENCE, PhyloXmlMapping.REFERENCE_DOI_ATTR, getDoi() );
        if ( !ForesterUtil.isEmpty( getValue() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.REFERENCE_DESC, getValue(), indentation );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.REFERENCE );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}