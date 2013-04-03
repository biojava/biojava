// $Id: WebLink.java,v 1.2 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.archaeopteryx;

import java.net.URL;

class WebLink {

    private final URL    _url;
    private final String _desc;
    private final String _source_identifier;

    WebLink( final URL url, final String desc, final String source_identifier ) {
        _url = url;
        _desc = desc;
        _source_identifier = source_identifier;
    }

    String getDesc() {
        return _desc;
    }

    String getSourceIdentifier() {
        return _source_identifier;
    }

    URL getUrl() {
        return _url;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( getDesc() );
        sb.append( " [" );
        sb.append( getSourceIdentifier() );
        sb.append( "]: " );
        sb.append( getUrl().toString() );
        return sb.toString();
    }
}
