// $Id: Uri.java,v 1.23 2008/09/24 16:42:51 cmzmasek Exp $
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
import java.net.URI;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;

public class Uri implements PhylogenyData {

    final private URI    _uri;
    final private String _description;
    final private String _type;

    public Uri( final String uri_str, final String description, final String type ) {
        if ( uri_str == null ) {
            throw new IllegalArgumentException( "attempt to create Uri from null" );
        }
        _uri = URI.create( uri_str );
        _description = description;
        _type = type;
    }

    public Uri( final URI uri ) {
        if ( uri == null ) {
            throw new IllegalArgumentException( "attempt to create Uri from null URI" );
        }
        _uri = uri;
        _description = "";
        _type = "";
    }

    public Uri( final URI uri, final String description, final String type ) {
        if ( uri == null ) {
            throw new IllegalArgumentException( "attempt to create Uri from null URI" );
        }
        _uri = uri;
        _description = description;
        _type = type;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue().toString() );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        sb.append( "[" );
        sb.append( getDescription() );
        sb.append( " " );
        sb.append( getType() );
        sb.append( "] " );
        sb.append( getValue().toString() );
        return sb;
    }

    @Override
    public PhylogenyData copy() {
        return new Uri( getValue().toString(), new String( getDescription() ), new String( getType() ) );
    }

    public String getDescription() {
        return _description;
    }

    public String getType() {
        return _type;
    }

    public URI getValue() {
        return _uri;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        PhylogenyDataUtil.appendElement( writer,
                                         PhyloXmlMapping.URI,
                                         getValue().toString(),
                                         PhyloXmlMapping.TYPE_ATTR,
                                         getType(),
                                         PhyloXmlMapping.URI_DESC_ATTR,
                                         getDescription(),
                                         indentation );
    }
}
