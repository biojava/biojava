// $Id: ParserBasedPhylogenyFactory.java,v 1.10 2008/08/19 21:34:21 cmzmasek Exp
// $
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

package org.forester.phylogeny.factories;

import java.io.IOException;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

public class ParserBasedPhylogenyFactory extends BasicPhylogenyFactory {

    private final static PhylogenyFactory _instance;
    static {
        try {
            _instance = new ParserBasedPhylogenyFactory();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private ParserBasedPhylogenyFactory() {
        // Private constructor.
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    public synchronized Phylogeny[] create( final Object source, final Object parser, final List<Object> parameters )
            throws IOException {
        if ( !( parser instanceof PhylogenyParser ) ) {
            throw new IllegalArgumentException( "attempt to use object of type other than PhylogenyParser as creator for ParserBasedPhylogenyFactory" );
        }
        final PhylogenyParser my_parser = ( PhylogenyParser ) parser;
        my_parser.setSource( source );
        return my_parser.parse();
    }

    public synchronized Phylogeny[] create( final Object source,
                                            final Object parser,
                                            final String schema_location,
                                            final List<Object> parameters ) throws IOException {
        if ( !( parser instanceof PhylogenyParser ) ) {
            throw new IllegalArgumentException( "attempt to use object of type other than PhylogenyParser as creator for ParserBasedPhylogenyFactory." );
        }
        if ( !( parser instanceof PhyloXmlParser ) ) {
            throw new IllegalArgumentException( "attempt to use schema location with other than phyloXML parser" );
        }
        final PhyloXmlParser xml_parser = ( PhyloXmlParser ) parser;
        if ( !ForesterUtil.isEmpty( schema_location ) ) {
            xml_parser.setValidateAgainstSchema( schema_location );
        }
        xml_parser.setSource( source );
        return xml_parser.parse();
    }

    public static PhylogenyFactory getInstance() {
        return _instance;
    }
}
