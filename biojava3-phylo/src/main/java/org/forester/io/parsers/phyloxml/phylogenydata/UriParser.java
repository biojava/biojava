// $Id: UriParser.java,v 1.3 2009/11/03 19:16:34 cmzmasek Exp $
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

package org.forester.io.parsers.phyloxml.phylogenydata;

import java.net.URI;
import java.net.URISyntaxException;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Uri;

public class UriParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new UriParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private UriParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhylogenyParserException {
        String type = "";
        String desc = "";
        URI uri = null;
        try {
            uri = new URI( element.getValueAsString() );
        }
        catch ( final URISyntaxException e ) {
            throw new PhylogenyParserException( "ill formatted Uri: " + element.getValueAsString() );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.URI_DESC_ATTR ) ) {
            desc = element.getAttribute( PhyloXmlMapping.URI_DESC_ATTR );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.TYPE_ATTR ) ) {
            type = element.getAttribute( PhyloXmlMapping.TYPE_ATTR );
        }
        return new Uri( uri, desc, type );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
